#' Multiview Learner Function
#'
#' This function generates synthetic multi-view data and performs various analyses, including plotting correlations and activity patterns.
#'
#' @param mainDir Character. The main directory where plots and data will be saved.
#' @param subDir Character. The subdirectory within \code{mainDir} for saving plots and data.
#' @param seed.number Numeric. Seed for random number generation.
#' @param M Integer. Number of views.
#' @param p_m Integer vector. Number of features per view.
#' @param n Integer. Train set size.
#' @param nTest Integer. Test set size.
#' @param params List. A list of parameters including:
#' \itemize{
#'   \item \code{fixed_act_pattern} Logical. Whether to use a fixed activity pattern.
#'   \item \code{K0} Integer. Number of factors.
#'   \item \code{nPatternJ} Integer. Number of joint patterns.
#'   \item \code{J0_m} Integer vector. Number of groups in the joint component per view.
#'   \item \code{K0_m} Integer vector. Number of factors in view-specific components.
#'   \item \code{S0_m} Integer vector. Number of groups in view-specific components.
#'   \item \code{nu2_0} Numeric. Variability in the joint loadings.
#'   \item \code{pi2_0} Numeric. Variability in the view-specific loadings.
#'   \item \code{snr_m} Numeric vector. Signal-to-noise ratio for each view.
#'   \item \code{snr_y} Numeric. Signal-to-noise ratio for the response.
#'   \item \code{decaying_loadings} Logical. Whether to enforce loadings decay.
#'   \item \code{decaying_type} Character. Type of decay ('sqrt' or 'log').
#'   \item \code{pJ_zero} Numeric. Probability of group-wise zero in joint component.
#'   \item \code{pJ_sign} Numeric. Probability of group-wise sign-switch in joint component.
#'   \item \code{pJ_zeroEntr} Numeric. Probability of entry-wise zero in joint component.
#'   \item \code{pS_zero} Numeric. Probability of group-wise zero in specific components.
#'   \item \code{pS_sign} Numeric. Probability of group-wise sign-switch in specific components.
#'   \item \code{pS_zeroEntr} Numeric. Probability of entry-wise zero in specific components.
#' }
#' @param Activity_pattern_plot Logical. Whether to save the activity pattern plot.
#' @param Emprical_correlation_plot Logical. Whether to save the empirical correlation plot.
#' @param Crossprd_plot Logical. Whether to save crossproduct plots.
#' @param Data_save Logical. Whether to save the generated data.
#'
#' @return A list containing the generated data and various parameters used in the simulation.
#'
#' @examples
#' mainDir <- tempdir()
#' subDir <- 'plots_simulated_data'
#' seed.number <- 1234
#' M <- 2
#' p_m <- c(500, 700)
#' n <- 100
#' nTest <- 100
#' params <- list(
#'   fixed_act_pattern = TRUE,
#'   K0 = 10,
#'   nPatternJ = 5,
#'   J0_m = c(4, 5),
#'   K0_m = c(5, 6),
#'   S0_m = c(4, 5),
#'   nu2_0 = 0.1,
#'   pi2_0 = 0.1,
#'   snr_m = NULL,
#'   snr_y = 1,
#'   decaying_loadings = FALSE,
#'   decaying_type = 'sqrt',
#'   pJ_zero = 0.6,
#'   pJ_sign = 0.4,
#'   pJ_zeroEntr = 0.5,
#'   pS_zero = 0.6,
#'   pS_sign = 0.4,
#'   pS_zeroEntr = 0.5
#' )
#' Activity_pattern_plot <- TRUE
#' Emprical_correlation_plot <- TRUE
#' Crossprd_plot <- TRUE
#' Data_save <- TRUE
#' Data <- MimeSys(mainDir, subDir, seed.number, M, p_m, n, nTest, params, Activity_pattern_plot, Emprical_correlation_plot, Crossprd_plot, Data_save)
#'
#' @export
# Function to safely close graphics devices
MimeSys <- function(
    mainDir = '/Users/siddheshkulkarni/Desktop',  # Main directory for results
    subDir = 'plots_simulated_data',  # Subdirectory under mainDir for outputs
    seed.number = 1234,  # Random seed for reproducibility
    sim_params = list(  # Parameters for manual simulation
      M = 3,  # Number of views
      p_m = c(500, 700, 1000),  # Number of features per view
      n = 100,  # Number of samples
      K0 = 10,  # Number of joint factors
      nPatternJ = 5,  # Number of joint patterns
      J0_m = c(5, 6, 7),  # Joint components per view
      K0_m = c(5, 7, 10),  # Specific factors per view
      S0_m = c(6, 7, 8),  # Specific components per view
      nu2_0 = 0.1,  # Variance for joint components
      pi2_0 = 0.1,  # Variance for specific components
      snr_m = c(1, 1, 1),  # SNR for each view
      snr_y = 1,  # SNR for response variable
      decaying_loadings = FALSE,  # Apply decaying loadings or not
      decaying_type = 'sqrt',  # Decay type ('sqrt' or 'log')
      pJ_zero = 0.6,  # Sparsity for joint components
      pJ_sign = 0.4,  # Proportion of negative joint signs
      pJ_zeroEntr = 0.5,  # Entry-level joint sparsity
      pS_zero = 0.6,  # Sparsity for specific components
      pS_sign = 0.4,  # Proportion of negative specific signs
      pS_zeroEntr = 0.5,  # Entry-level specific sparsity
      fixed_act_pattern = "YES"  # Whether to use fixed activity pattern
    ),
    user_params = list(  # Parameters for user-supplied data
      use_jive = TRUE,  # Use JIVE for analysis
      user_data = Data,  # User-provided data
      snr_data = 1,  # SNR for user data
      Heatmap_indicator = "YES"  # Whether to generate heatmaps
    ),
    Activity_pattern_plot = TRUE,  # Plot activity patterns or not
    Empirical_correlation_plot = TRUE,  # Plot empirical correlations
    Crossprd_plot = TRUE,  # Plot cross-product matrices or not
    Data_save = TRUE  # Whether to save generated data
) {
  # Helper function to safely close devices
  safe_dev_off <- function() {
    while (!is.null(dev.list())) {
      dev.off()
    }
  }

  # Helper function to cluster a matrix and return the order based on hierarchical clustering
  cluster_matrix <- function(mat) {
    # Ensure all values are finite
    mat[is.na(mat) | is.nan(mat) | is.infinite(mat)] <- 0

    # Compute the distance matrix
    dist_mat <- dist(mat)

    # Perform hierarchical clustering
    hc <- hclust(dist_mat)

    # Return the order based on clustering
    return(order.dendrogram(as.dendrogram(hc)))
  }


  # Helper function to plot clustered heatmaps
  heatmap_cor_matrix_clustered <- function(cor_matrix, title = "Heatmap") {
    # Perform hierarchical clustering on rows and columns
    row_order <- cluster_matrix(cor_matrix)
    col_order <- cluster_matrix(t(cor_matrix))

    # Reorder the matrix
    cor_matrix <- cor_matrix[row_order, col_order]

    # Plot the heatmap using ggplot2
    cor_data <- reshape2::melt(cor_matrix)
    ggplot2::ggplot(cor_data, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
      labs(title = title, x = "Variables", y = "Variables") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  # Create the directory to save the results and set random seed
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  set.seed(seed.number)

  # # Initialize SNR parameters if not provided
  # M_data <- length(user_params$user_data)
  # p_m_data <- lapply(user_params$user_data, nrow)
  # n_data <- lapply(user_params$user_data, ncol)
  #
  # if (is.null(user_params$snr_data)) {
  #   user_params$snr <- rep(1, M_data)
  # }

  # JIVE functionality - when user-provided data is used
  # Initialize SNR parameters if not provided
  if (user_params$use_jive && !is.null(user_params$user_data)) {
    M_data <- length(user_params$user_data)  # User-provided data
    p_m_data <- lapply(user_params$user_data, nrow)
    n_data <- lapply(user_params$user_data, ncol)
    if (is.null(user_params$snr_data)) {
      user_params$snr_data <- rep(1, M_data)
    }
  } else {
    M <- sim_params$M  # Simulated data
    p_m <- sim_params$p_m
    n <- sim_params$n
  }


  if (user_params$use_jive && !is.null(user_params$user_data)){
    if (!requireNamespace("r.jive", quietly = TRUE)) {
      stop("The r.jive package is required but not installed.")
    }
    library(r.jive)

    # Apply JIVE on the user-provided data
    jive_result <- jive(user_params$user_data)

    # Use jive.predict to extract joint and individual loadings
    jive_predicted <- jive.predict(user_params$user_data, jive_result)

    # Optimal factors for joint and individual components
    optimal_factors_joint <- dim(jive_predicted$joint.load)[2]
    optimal_factors_individual <- lapply(jive_predicted$indiv.load, ncol)

    # Disintegrate joint matrix into separate matrices for each view
    joint_matrix <- jive_predicted$joint.load
    N <- p_m_data
    joint_load <- list()
    start_index <- 1
    for (i in 1:length(N)) {
      num_rows <- N[[i]]
      joint_load[[i]] <- joint_matrix[start_index:(start_index + num_rows - 1), ]
      start_index <- start_index + num_rows
    }

    # Generate synthetic data from joint and individual loadings
    X_m <- lapply(1:M_data, function(m) {
      t(
        tcrossprod(joint_load[[m]], matrix(rnorm(n_data[[m]] * jive_result$rankJ), nrow = n_data[[m]])) +
          tcrossprod(jive_predicted$indiv.load[[m]], matrix(rnorm(n_data[[m]] * jive_result$rankA[m]), nrow = n_data[[m]])) +
          sqrt(user_params$snr[m]) * matrix(rnorm(p_m_data[[m]] * n_data[[m]]), ncol = n_data[[m]])  # Add noise
      )
    })

    # Generate heatmaps if indicated
    if (user_params$Heatmap_indicator == "YES") {
      png(file.path(mainDir, subDir, "Heatmaps.png"), height = 465, width = 705)
      showHeatmaps(jive_result)
      safe_dev_off()
    }

  } else {
    # Manual data generation if JIVE is not used
    M <- sim_params$M
    p_m <- sim_params$p_m
    n <- sim_params$n
    K0 <- sim_params$K0
    nPatternJ <- sim_params$nPatternJ
    J0_m <- sim_params$J0_m
    K0_m <- sim_params$K0_m
    S0_m <- sim_params$S0_m
    nu2_0 <- sim_params$nu2_0
    pi2_0 <- sim_params$pi2_0
    snr_m <- sim_params$snr_m
    snr_y <- sim_params$snr_y
    decaying_loadings <- sim_params$decaying_loadings
    decaying_type <- sim_params$decaying_type
    pJ_zero <- sim_params$pJ_zero
    pJ_sign <- sim_params$pJ_sign
    pJ_zeroEntr <- sim_params$pJ_zeroEntr
    pS_zero <- sim_params$pS_zero
    pS_sign <- sim_params$pS_sign
    pS_zeroEntr <- sim_params$pS_zeroEntr
    fixed_act_pattern <- sim_params$fixed_act_pattern

    # Generate activity pattern if fixed_act_pattern is TRUE
    if (fixed_act_pattern == "YES") {
      which_comb <- expand.grid(replicate(M + 1, c(0, 1), simplify = FALSE))
      which_comb <- apply(which_comb[rowSums(which_comb) > 1, ], 1, paste, collapse = "")
      which_comb <- sample(which_comb, nPatternJ)
      value_comb <- getBlocksDim(nPatternJ, K0, alpha = 5)
      act_J <- sapply(which_comb, function(v) as.numeric(unlist(strsplit(v, ""))))
      act_J <- t(apply(act_J, 1, function(v) rep(v, value_comb)))
    } else {
      stop("Non-fixed activity patterns are not supported.")
    }

    # Generate dimensions and sparsity patterns for joint and specific components
    dim_J <- lapply(1:M, function(m) getBlocksDim(J0_m[m], p_m[m]))
    dim_S <- lapply(1:M, function(m) getBlocksDim(S0_m[m], p_m[m]))

    # Initialize group-wise sparsity patterns for joint components
    group_sprs_J <- list()
    for (m in 1:M) {
      group_sprs_J[[m]] <- matrix(rbinom(J0_m[m] * K0, 1, 1 - pJ_zero * (sum(act_J[m, ]) / K0)), ncol = K0) *
        (2 * matrix(rbinom(J0_m[m] * K0, 1, 1 - pJ_sign), ncol = K0) - 1)

      # Correct sparsity patterns to ensure no group is entirely zero
      col_wrong <- apply(group_sprs_J[[m]], 2, function(v) min(v == 0)) * act_J[m, ]
      while (sum(col_wrong) > 0) {
        group_sprs_J[[m]][, as.logical(col_wrong)] <- matrix(rbinom(J0_m[m] * sum(col_wrong), 1, 1 - pJ_zero * (sum(act_J[m, ]) / K0)), ncol = sum(col_wrong)) *
          (2 * matrix(rbinom(J0_m[m] * sum(col_wrong), 1, 1 - pJ_sign), ncol = sum(col_wrong)) - 1)
        col_wrong <- apply(group_sprs_J[[m]], 2, function(v) min(v == 0)) * act_J[m, ]
      }
    }

    # Initialize entry-wise sparsity patterns for joint components
    elem_sprs_J <- lapply(1:M, function(m) matrix(rbinom(p_m[m] * K0, 1, 1 - pJ_zeroEntr), ncol = K0))
    for (m in 1:M) {
      sprs_tot <- elem_sprs_J[[m]] * apply(group_sprs_J[[m]], 2, function(v) rep(v, dim_J[[m]]))
      if (max(apply(sprs_tot, 1, function(v) min(v == 0))) > 0) {
        elem_sprs_J[[m]][which(apply(sprs_tot, 1, function(v) min(v == 0)) > 0), act_J[m, ]] <- 1
      }
    }

    # Initialize group-wise sparsity patterns for specific components
    group_sprs_S <- list()
    for (m in 1:M) {
      group_sprs_S[[m]] <- matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_zero), ncol = K0_m[m]) *
        (2 * matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_sign), ncol = K0_m[m]) - 1)

      # Correct sparsity patterns
      while (max(apply(group_sprs_S[[m]], 1, function(v) min(v == 0))) > 0) {
        group_sprs_S[[m]] <- matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_zero), ncol = K0_m[m]) *
          (2 * matrix(rbinom(S0_m[m] * K0_m[m], 1, 1 - pS_sign), ncol = K0_m[m]) - 1)
      }
    }

    # Generate synthetic data using manually generated loadings and specific components
    X_m <- lapply(1:M, function(m) {
      t(tcrossprod(load_J[[m]], matrix(rnorm(n * K0), nrow = n)) +
          tcrossprod(load_S[[m]], matrix(rnorm(n * K0_m[m]), nrow = n)) +
          sqrt(s2m[[m]]) * matrix(rnorm(p_m[m] * n), ncol = n))  # Add noise
    })
  }

  # Save the generated data
  if (Data_save) {
    if (user_params$use_jive && !is.null(user_params$user_data)) {

      eta <- matrix(rnorm(n_data[[1]] * optimal_factors_joint), nrow = n_data[[1]])

      phi_m <- lapply(optimal_factors_individual, function(kk) matrix(rnorm(n_data[[1]] * kk), nrow = n_data[[1]]))

      Theta <- (2 * rbinom(optimal_factors_joint, 1, 0.5) - 1) * rbeta(optimal_factors_joint, shape1 = 5, shape2 = 3) / optimal_factors_joint
      s2y <- sum(Theta^2) / user_params$snr_data
      y <- eta %*% Theta + sqrt(s2y) * rnorm(n_data[[1]])
      colnames(y) <- c("response")

      preprocess_y <- caret::preProcess(y, method = c("center", "scale"))

      y <- as.matrix(predict(preprocess_y, y))

      preprocess_X_m <- list()
      for (m in 1:M_data) {
        colnames(X_m[[m]]) <- seq(ncol(X_m[[m]]))
        preprocess_X_m[[m]] <- caret::preProcess(X_m[[m]], method = c("center", "scale"))
        X_m[[m]] <- as.matrix(predict(preprocess_X_m[[m]], X_m[[m]]))
      }
    }

    Data <- list(M=M_data,
                 n = n_data,
                 p_m = p_m_data,
                 y = y,
                 X_m = X_m,
                 preprocess_X_m = preprocess_X_m,
                 preprocess_y = preprocess_y,
                 s2y = s2y,
                 s2m = user_params$snr,
                 Theta = Theta,
                 Lambda_m = joint_load,
                 Gamma_m = jive_predicted$indiv.load,
                 eta = eta,
                 phi_m = phi_m)

    saveRDS(Data, file = file.path(mainDir, "Simulated_multiview.rds"))

  } else {

    eta <- matrix(rnorm(n * K0), nrow = n)

    phi_m <- lapply(K0_m, function(kk) matrix(rnorm(n * kk), nrow = n))

    Theta <- (2 * rbinom(K0, 1, 0.5) - 1) * rbeta(K0, shape1 = 5, shape2 = 3) / sqrt(K0)

    Theta <- unname(Theta * act_J[nrow(act_J), ])

    s2y <- sum(Theta^2) / snr_y
    y <- eta %*% Theta + sqrt(s2y) * rnorm(n)
    colnames(y) <- c("response")

    preprocess_y <- caret::preProcess(y, method = c("center", "scale"))

    y <- as.matrix(predict(preprocess_y, y))

    preprocess_X_m <- list()
    for (m in 1:M) {
      colnames(X_m[[m]]) <- seq(ncol(X_m[[m]]))
      preprocess_X_m[[m]] <- caret::preProcess(X_m[[m]], method = c("center", "scale"))
      X_m[[m]] <- as.matrix(predict(preprocess_X_m[[m]], X_m[[m]]))
    }


    Data <- list(M=M,
                 n = n,
                 p_m = p_m,
                 y = y,
                 X_m = X_m,
                 preprocess_X_m = preprocess_X_m,
                 preprocess_y = preprocess_y,
                 s2y = s2y,
                 s2m = s2m,
                 Theta = Theta,
                 Lambda_m = load_J,
                 Gamma_m = load_S,
                 eta = eta,
                 phi_m = phi_m)

    saveRDS(Data, file = file.path(mainDir, "Simulated_multiview.rds"))

  }
  #library(ggplot2)
  # Correlation comparison and visualization
  # Function to compare original and simulated correlations and save metrics as JPEG
  # Function to compare original and simulated correlations and save metrics as JPEG
  for (m in 1:M_data) {

    # Remove variable names for clean correlation matrices
    original_cor <- cor(user_params$user_data[[m]])
    rownames(original_cor) <- colnames(original_cor) <- NULL

    simulated_cor <- cor(t(X_m[[m]]))
    rownames(simulated_cor) <- colnames(simulated_cor) <- NULL

    # Compute Frobenius norm
    frobenius_norm <- sqrt(sum((original_cor - simulated_cor)^2))

    # Error handling for Kullback-Leibler Divergence and LogDet Divergence
    kl_divergence <- tryCatch({
      0.5 * (sum(diag(solve(simulated_cor) %*% original_cor)) - log(det(solve(simulated_cor) %*% original_cor)) - ncol(original_cor))
    }, error = function(e) {
      return("Couldn't calculate KL Divergence due to singularity.")
    })

    logdet_divergence <- tryCatch({
      log(det(original_cor)) - log(det(simulated_cor)) - sum(diag((solve(simulated_cor) %*% (original_cor - simulated_cor))))
    }, error = function(e) {
      return("Couldn't calculate LogDet Divergence due to singularity.")
    })

    # Create metrics text
    metrics_text <- paste0(
      "View ", m, " Comparison Metrics:\n",
      "Frobenius Norm: ", round(frobenius_norm, 4), "\n",
      "KL Divergence: ", ifelse(is.numeric(kl_divergence), round(kl_divergence, 4), kl_divergence), "\n",
      "LogDet Divergence: ", ifelse(is.numeric(logdet_divergence), round(logdet_divergence, 4), logdet_divergence), "\n"
    )

    # Plot original and simulated correlations using heatmaps
    plot_original <- heatmap_cor_matrix_clustered(original_cor, title = paste0('Original Correlation - View ', m))
    plot_simulated <- heatmap_cor_matrix_clustered(simulated_cor, title = paste0('Simulated Correlation - View ', m))

    # Open JPEG device for saving the heatmaps and metrics together
    jpeg(file.path(mainDir, subDir, paste0('comparison_metrics_view_', m, '.jpg')), width = 1200, height = 800, res = 150)

    # Arrange heatmaps side by side with the metrics text at the bottom
    gridExtra::grid.arrange(
      plot_original,
      plot_simulated,
      grid::textGrob(metrics_text, x = 0.5, hjust = 0.5, gp = grid::gpar(fontsize = 8)),  # Text at the bottom
      ncol = 2,  # Arrange in 2 columns
      layout_matrix = rbind(c(1, 2), c(3, 3))  # Put the text in the third row spanning both columns
    )

    # Close the JPEG device to save the file
    safe_dev_off()

  }

  return(Data)
}


