#' Generate simulated multi-study, multiview data for benchmarking
#'
#' Generates *multi-study* simulated multiview datasets with specified parameters,
#' including per-study sample size, signal-to-noise ratio, DE probabilities,
#' and response variable generation mode. It also splits the data into
#' training and testing sets within each study.
#'
#' The model uses a Craig-style shared/individual coefficient structure:
#' for study s,
#'   beta_s = sqrt(rho.beta) * beta_shared + sqrt(1 - rho.beta) * beta_indiv_s,
#' where beta_shared is common across studies and beta_indiv_s is study-specific.
#'
#' Special cases:
#'   - nstudy = 1, rho.beta = 1, sigma.alpha = 0, tau.snr = 0:
#'       => reduces to your current single-study generator (one study).
#'   - rho.beta = 0 (any nstudy):
#'       => studies are independent (no shared beta).
#'
#' @param nsample Number of samples per study (integer)
#' @param nstudy Number of studies (integer, default 1)
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
#' @param rho.beta Shared-vs-individual mixing in [0,1].
#'        1 = fully shared beta across studies; 0 = independent betas.
#' @param sigma.alpha SD of study-specific intercepts (0 = no shift)
#' @param tau.snr SD of log-SNR across studies (0 = same SNR)
#' @param nrep Number of repetitions
#' @param seed Random seed
#'
#' @return List containing training and testing datasets and simulation parameters:
#'   \itemize{
#'     \item trainDat: list of length nrep; each element is a list of length nstudy,
#'           with one pcl-like object per study containing training samples.
#'     \item testDat: same, for test samples.
#'     \item snr, p.train, de.prob, de.downProb, de.facLoc, de.facScale,
#'           nrep, seed, nstudy, rho.beta, sigma.alpha, tau.snr
#'   }
#' @export
gen_simmba_multistudy <- function(nsample,                      # per-study sample size
                                  nstudy = 1,                   # number of studies
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
                                  rho.beta = 1,                 # shared vs individual beta
                                  sigma.alpha = 0,              # study-specific intercept SD
                                  tau.snr = 0,                  # study-specific log-SNR SD
                                  nrep = 100,                   # Number of repetitions
                                  seed = 1234) {                # Random seed

  set.seed(seed)

  ygen.mode    <- match.arg(ygen.mode)
  outcome.type <- match.arg(outcome.type)

  # Outer lists: over repetitions
  trainDat <- testDat <- vector("list", nrep)
  names(trainDat) <- names(testDat) <- paste("Rep", 1:nrep, sep = "_")

  ######################################################################
  # SPECIAL CASE: recover original single-study behavior exactly-ish   #
  #   nstudy = 1, rho.beta = 1, sigma.alpha = 0, tau.snr = 0           #
  ######################################################################
  if (nstudy == 1 && rho.beta == 1 && sigma.alpha == 0 && tau.snr == 0) {

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

      ## 3. Base linear predictor from beta
      eta_lin <- as.numeric(X %*% beta0)

      ## Optional non-linear Friedman component
      if (ygen.mode %in% c("Friedman", "Friedman2")) {

        nonzero_index <- which(beta0 != 0)
        if (length(nonzero_index) < 5) {
          stop("Not enough non-zero coefficients to select 5 Friedman features.")
        }
        friedman_index <- sample(nonzero_index, 5)
        X.friedman <- X[, friedman_index, drop = FALSE]

        f <- function(x) {
          10 * sin(pi * x[, 1] * x[, 2]) +
            20 * (x[, 3] - 0.5)^2 +
            10 * x[, 4] +
            5  * x[, 5]
        }

        Xbeta.friedman <- f(X.friedman)
      }

      ## 4. Final linear predictor eta
      if (ygen.mode == "LM") {
        eta <- eta_lin
      } else if (ygen.mode == "Friedman") {
        eta <- Xbeta.friedman
      } else if (ygen.mode == "Friedman2") {
        eta <- eta_lin + Xbeta.friedman
      }

      pcl$sample_metadata$Xbeta <- eta

      ###################################
      # 5. Generate outcome from eta    #
      ###################################
      if (outcome.type == "continuous") {

        sigma2 <- as.vector(var(eta) / snr)
        Y <- eta + rnorm(nsample) * sqrt(sigma2)
        pcl$sample_metadata$Y <- as.vector(Y)

      } else if (outcome.type == "binary") {

        p <- plogis(eta)
        Ybin <- rbinom(nsample, size = 1, prob = p)
        pcl$sample_metadata$Y <- Ybin

      } else if (outcome.type == "survival") {

        h  <- as.vector(surv.hscale * exp(eta))
        X0 <- rexp(nsample, rate = h)

        C <- runif(nsample, cens.lower, cens.upper)

        delta <- ifelse(C >= X0, 1L, 0L)
        Tobs  <- ifelse(C >= X0, X0, C)

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

      # wrap in a list of one study to keep structure [rep][study]
      trainDat[[k]] <- list(Study_1 = train)
      testDat[[k]]  <- list(Study_1 = test)
    }

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
      cens.upper  = cens.upper,
      nstudy      = nstudy,
      rho.beta    = rho.beta,
      sigma.alpha = sigma.alpha,
      tau.snr     = tau.snr
    ))
  }

  #############################################################
  # GENERAL CASE: true multi-study with shared/individual betas
  #############################################################

  for (k in seq_len(nrep)) {

    ## --- 1. Draw a "template" InterSIM dataset to define feature structure
    pcl_template <- trigger_InterSIM(n = nsample)
    X_template   <- pcl_template$feature_table %>% t() %>% as.matrix()
    nfeature     <- table(pcl_template$feature_metadata$featureType)

    # number of features
    p <- nrow(pcl_template$feature_table)

    ## --- 2. Shared beta via Splatter (as in original)
    de.facs.shared <- vector("list", 3)
    for (i in 1:3) {
      de.facs.shared[[i]] <- splatter:::getLNormFactors(
        n.facs   = nfeature[i],
        sel.prob = de.prob[i],
        neg.prob = de.downProb[i],
        fac.loc  = de.facLoc[i],
        fac.scale = de.facScale[i]
      )
    }
    beta_shared <- log2(unlist(de.facs.shared))   # length p

    # Lists over studies for this replication
    train_list <- vector("list", nstudy)
    test_list  <- vector("list", nstudy)
    names(train_list) <- names(test_list) <- paste0("Study_", seq_len(nstudy))

    for (s in seq_len(nstudy)) {

      # --- 3. Study-specific beta_indiv and mixing with shared
      if (rho.beta < 1) {
        de.facs.indiv <- vector("list", 3)
        for (i in 1:3) {
          de.facs.indiv[[i]] <- splatter:::getLNormFactors(
            n.facs   = nfeature[i],
            sel.prob = de.prob[i],
            neg.prob = de.downProb[i],
            fac.loc  = de.facLoc[i],
            fac.scale = de.facScale[i]
          )
        }
        beta_indiv <- log2(unlist(de.facs.indiv))
      } else {
        beta_indiv <- rep(0, length(beta_shared))
      }

      # Craig-style mixture of shared and individual
      beta_s <- sqrt(rho.beta) * beta_shared +
        sqrt(1 - rho.beta) * beta_indiv

      # Study-specific intercept and SNR
      alpha_s <- if (sigma.alpha > 0) stats::rnorm(1, 0, sigma.alpha) else 0

      if (tau.snr > 0) {
        log_snr_s <- log(snr) + stats::rnorm(1, 0, tau.snr)
      } else {
        log_snr_s <- log(snr)
      }
      snr_s <- exp(log_snr_s)

      # --- 4. Generate study-specific features
      # Re-use template pcl for Study 1 to guarantee alignment; new calls for others
      if (s == 1) {
        pcl <- pcl_template
      } else {
        pcl <- trigger_InterSIM(n = nsample)
      }

      X <- pcl$feature_table %>% t() %>% as.matrix()
      n_s <- nrow(X)

      # base linear predictor
      eta_lin <- as.numeric(X %*% beta_s) + alpha_s

      # Friedman components if requested
      if (ygen.mode %in% c("Friedman", "Friedman2")) {

        nonzero_index <- which(beta_s != 0)
        if (length(nonzero_index) < 5) {
          stop("Not enough non-zero coefficients to select 5 Friedman features.")
        }
        friedman_index <- sample(nonzero_index, 5)
        X.friedman <- X[, friedman_index, drop = FALSE]

        f <- function(x) {
          10 * sin(pi * x[, 1] * x[, 2]) +
            20 * (x[, 3] - 0.5)^2 +
            10 * x[, 4] +
            5  * x[, 5]
        }

        Xbeta.friedman <- f(X.friedman)
      }

      # Final eta
      if (ygen.mode == "LM") {
        eta <- eta_lin
      } else if (ygen.mode == "Friedman") {
        eta <- Xbeta.friedman
      } else if (ygen.mode == "Friedman2") {
        eta <- eta_lin + Xbeta.friedman
      }

      pcl$sample_metadata$Xbeta <- eta
      pcl$sample_metadata$study <- paste0("Study_", s)

      ###################################
      # 5. Generate outcome from eta    #
      ###################################
      if (outcome.type == "continuous") {

        sigma2 <- as.vector(stats::var(eta) / snr_s)
        Y <- eta + stats::rnorm(n_s) * sqrt(sigma2)
        pcl$sample_metadata$Y <- as.vector(Y)

      } else if (outcome.type == "binary") {

        p_prob <- stats::plogis(eta)
        Ybin <- stats::rbinom(n_s, size = 1, prob = p_prob)
        pcl$sample_metadata$Y <- Ybin

      } else if (outcome.type == "survival") {

        h  <- as.vector(surv.hscale * exp(eta))
        X0 <- stats::rexp(n_s, rate = h)

        C <- stats::runif(n_s, cens.lower, cens.upper)

        delta <- ifelse(C >= X0, 1L, 0L)
        Tobs  <- ifelse(C >= X0, X0, C)

        pcl$sample_metadata$time   <- Tobs
        pcl$sample_metadata$status <- delta
      }

      #################################
      # 6. Train / test split         #
      #################################
      train <- test <- pcl
      tr.row <- sample.int(n_s, size = round(n_s * p.train), replace = FALSE)

      train$sample_metadata <- pcl$sample_metadata[tr.row, , drop = FALSE]
      test$sample_metadata  <- pcl$sample_metadata[-tr.row, , drop = FALSE]

      train$feature_table <- pcl$feature_table[, tr.row, drop = FALSE]
      test$feature_table  <- pcl$feature_table[, -tr.row, drop = FALSE]

      train_list[[s]] <- train
      test_list[[s]]  <- test
    }

    trainDat[[k]] <- train_list
    testDat[[k]]  <- test_list
  }

  # Return synthetic multi-study data and parameters
  list(
    trainDat    = trainDat,    # [rep][study]
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
    cens.upper  = cens.upper,
    nstudy      = nstudy,
    rho.beta    = rho.beta,
    sigma.alpha = sigma.alpha,
    tau.snr     = tau.snr
  )
}
