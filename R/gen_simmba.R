#' Generate simulated multi-omics data for benchmarking
#'
#' Generates simulated multi-omics datasets with specified parameters,
#' including sample size, signal-to-noise ratio, DE probabilities,
#' and response variable generation mode. It also splits the data into
#' training and testing sets.
#'
#' @param nsample Sample size
#' @param snr Signal to noise ratio (for continuous outcomes)
#' @param p.train Train-test split ratio
#' @param de.prob DE probability across all modalities (vector of length 3)
#' @param de.downProb Down-regulation probability (vector of length 3)
#' @param de.facLoc DE factor location (vector of length 3)
#' @param de.facScale DE factor scale (vector of length 3)
#' @param ygen.mode Y generation mode: "LM", "Friedman", or "Friedman2"
#' @param outcome.type Outcome type: "continuous", "binary", or "survival"
#' @param surv.hscale Multiplicative scale on hazard for survival outcomes
#'   (hazard_i = surv.hscale * exp(eta_i))
#' @param cens.lower Lower bound for uniform censoring time
#' @param cens.upper Upper bound for uniform censoring time
#' @param nrep Number of repetitions
#' @param seed Random seed
#'
#' @return List containing training and testing datasets and simulation parameters:
#'   \itemize{
#'     \item trainDat: list of length nrep, each an object like pcl with training samples
#'     \item testDat: list of length nrep, same structure for test samples
#'     \item snr, p.train, de.prob, de.downProb, de.facLoc, de.facScale, nrep, seed
#'   }
#' @export
#' @examples
#' simulated_data <- gen_simmba(nsample = 100)
gen_simmba <- function(nsample,                      # Sample size
                       snr = 1,                      # Signal to noise ratio (continuous)
                       p.train = 0.7,                # Train-test split ratio
                       de.prob = rep(0.1, 3),        # DE probability across modalities
                       de.downProb = rep(0.5, 3),    # Down-regulation probability
                       de.facLoc = rep(1, 3),        # DE factor location
                       de.facScale = rep(0.4, 3),    # DE factor scale
                       ygen.mode = c("LM", "Friedman", "Friedman2"),
                       outcome.type = c("continuous", "binary", "survival"),
                       surv.hscale = 1,              # hazard scale for survival
                       cens.lower = 1,               # censoring lower bound
                       cens.upper = 3,               # censoring upper bound
                       nrep = 100,                   # Number of repetitions
                       seed = 1234) {                # Random seed

  set.seed(seed)

  ygen.mode    <- match.arg(ygen.mode)
  outcome.type <- match.arg(outcome.type)

  # Initialize lists to store training and testing datasets
  trainDat <- testDat <- vector("list", nrep)
  names(trainDat) <- names(testDat) <- paste("Rep", 1:nrep, sep = "_")

  # Loop over the number of repetitions
  for (k in seq_len(nrep)) {

    ## 1. Generate feature data using InterSIM
    pcl <- trigger_InterSIM(n = nsample)  # user-defined wrapper
    X <- pcl$feature_table %>% t() %>% as.matrix()

    ## 2. Construct beta0 using Splatter factors (kept fixed afterwards)
    nfeature <- table(pcl$feature_metadata$featureType)

    de.facs <- vector("list", 3)
    for (i in 1:3) {
      de.facs[[i]] <- splatter:::getLNormFactors(
        n.facs   = nfeature[i],
        sel.prob = de.prob[i],
        neg.prob = de.downProb[i],
        fac.loc  = de.facLoc[i],
        fac.scale = de.facScale[i]
      )
    }

    # Convert to log fold changes (LFCs): this is your beta
    beta0 <- log2(unlist(de.facs))

    ## 3. Base linear predictor from beta (this is what you want to preserve)
    eta_lin <- as.numeric(X %*% beta0)

    ## Optional non-linear Friedman component (beta0 itself is unchanged)
    if (ygen.mode %in% c("Friedman", "Friedman2")) {

      # Choose 5 non-zero beta features for the Friedman function
      nonzero_index <- which(beta0 != 0)
      if (length(nonzero_index) < 5) {
        stop("Not enough non-zero coefficients to select 5 Friedman features.")
      }
      friedman_index <- sample(nonzero_index, 5)
      X.friedman <- X[, friedman_index, drop = FALSE]

      # Friedman function (embedded here for self-containment)
      f <- function(x) {
        10 * sin(pi * x[, 1] * x[, 2]) +
          20 * (x[, 3] - 0.5)^2 +
          10 * x[, 4] +
          5  * x[, 5]
      }

      Xbeta.friedman <- f(X.friedman)
    }

    ## 4. Final linear predictor eta, used for *all* outcome types
    if (ygen.mode == "LM") {
      eta <- eta_lin
    } else if (ygen.mode == "Friedman") {
      eta <- Xbeta.friedman
    } else if (ygen.mode == "Friedman2") {
      eta <- eta_lin + Xbeta.friedman
    }

    # Store Xbeta in metadata for reference
    pcl$sample_metadata$Xbeta <- eta

    ###################################
    # 5. Generate outcome from eta    #
    ###################################

    if (outcome.type == "continuous") {

      # Y = eta + noise, SNR controlled by snr
      sigma2 <- as.vector(var(eta) / snr)
      Y <- eta + rnorm(nsample) * sqrt(sigma2)
      pcl$sample_metadata$Y <- as.vector(Y)

    } else if (outcome.type == "binary") {

      # Bernoulli with probability p = logit^{-1}(eta)
      p <- plogis(eta)
      Ybin <- rbinom(nsample, size = 1, prob = p)
      pcl$sample_metadata$Y <- Ybin

    } else if (outcome.type == "survival") {

      # Simple survival model with exponential event times:
      # hazard_i = surv.hscale * exp(eta_i)
      h  <- as.vector(surv.hscale * exp(eta))   # rate parameter for rexp
      X0 <- rexp(nsample, rate = h)             # true event times

      # Censoring times
      C <- runif(nsample, cens.lower, cens.upper)

      # Event indicator and observed time
      delta <- ifelse(C >= X0, 1L, 0L)
      Tobs  <- ifelse(C >= X0, X0, C)           # NOTE: uses X0, not X

      pcl$sample_metadata$time   <- Tobs
      pcl$sample_metadata$status <- delta
    }

    #################################
    # 6. Train / test split         #
    #################################

    train <- test <- pcl
    tr.row <- sample.int(nsample, size = round(nsample * p.train), replace = FALSE)

    train$sample_metadata <- pcl$sample_metadata[tr.row, , drop = FALSE]
    test$sample_metadata  <- pcl$sample_metadata[-tr.row, , drop = FALSE]

    train$feature_table <- pcl$feature_table[, tr.row, drop = FALSE]
    test$feature_table  <- pcl$feature_table[, -tr.row, drop = FALSE]

    # Store in lists
    trainDat[[k]] <- train
    testDat[[k]]  <- test
  }

  # Return synthetic data and simulation parameters
  return(list(
    trainDat    = trainDat,
    testDat     = testDat,
    snr         = snr,
    p.train     = p.train,
    de.prob     = de.prob,
    de.downProb = de.downProb,
    de.facLoc   = de.facLoc,
    de.facScale = de.facScale,
    nrep        = nrep,
    seed        = seed,
    ygen.mode   = ygen.mode,
    outcome.type = outcome.type,
    surv.hscale = surv.hscale,
    cens.lower  = cens.lower,
    cens.upper  = cens.upper
  ))
}
