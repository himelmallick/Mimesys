#' This function generates simulated multi-omics data.
#'
#'Generates simulated multi-omics datasets with specified parameters,
#'including sample size, signal-to-noise ratio, DE probabilities,
#'and response variable generation mode. It also splits the data into training and testing sets.
#'
#' @param nsample Sample size
#' @param snr Signal to noise ratio
#' @param p.train Train-test split ratio
#' @param de.prob DE probability across all modalities
#' @param de.downProb Down-regulation probability
#' @param de.facLoc DE factor location
#' @param de.facScale DE factor scale
#' @param ygen.mode Y generation mode
#' @param nrep Number of repetitions
#' @param seed Random seed
#' @return List containing training and testing datasets
#' @export
#' @examples
#' simulated_data <- gen_simmba(nsample = 100)

########################################################
# Generate simulated multi-omics data for benchmarking #
########################################################

# Load necessary libraries
# library(InterSIM)  # For simulating inter-connected omics layers
# library(tidyverse)  # For data manipulation and visualization
# library(splatter)  # For generating realistic effect sizes
# library(NMF)  # For non-negative matrix factorization
# library(stringr)  # For string manipulation

# # Ensure BiocManager is installed and then install Biobase and splatter packages if not already installed
#
# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")
# if(!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("splatter")

#########################
# A wrapper for InterSIM #
#########################

# Function to trigger InterSIM to generate null multi-omics dataset

# Define the number of samples

# trigger_InterSIM <- function(n) {
#
#   # Generate null multi-omics dataset with specified parameters
#   sim.data <- InterSIM(n.sample = n,
#                        cluster.sample.prop = c(0.51, 0.49), # Sample proportion for clustering
#                        delta.methyl = 0, # Set delta to zero for null data simulation
#                        delta.expr = 0, # Set delta to zero for null data simulation
#                        delta.protein = 0) # Set delta to zero for null data simulation
#
#   # Process features into a list
#   names <- c('methyl', 'expr', 'protein')  # Define feature types
#   list_features <- vector("list", length(names))  # Create an empty list for features
#   names(list_features) <- names  # Name the list elements
#
#   # Transpose and assign simulated data to the list
#   list_features[[1]] <- t(sim.data$dat.methyl)
#   list_features[[2]] <- t(sim.data$dat.expr)
#   list_features[[3]] <- t(sim.data$dat.protein)
#
#   # Combine the features into a single data frame
#   feature_table <- as.data.frame(Reduce(rbind, list_features))
#   colnames(feature_table) <- stringr::str_to_title(colnames(feature_table))  # Format column names
#
#   # Extract and format sample metadata
#   sample_metadata <- sim.data$clustering.assignment
#   colnames(sample_metadata) <- c('subjectID', 'Y')
#   sample_metadata$subjectID <- str_to_title(sample_metadata$subjectID)
#   rownames(sample_metadata) <- sample_metadata$subjectID
#
#   # Create feature metadata
#   rowID <- rep(names, sapply(list_features, nrow))
#   feature_metadata <- cbind.data.frame(featureID = rownames(feature_table), featureType = rowID)
#   rownames(feature_metadata) <- feature_metadata$featureID
#
#   # Save datasets as a list and return it
#   pcl <- list(feature_table = feature_table,
#               sample_metadata = sample_metadata,
#               feature_metadata = feature_metadata)
#
#   return(pcl)
# }
#
# # Friedman function to generate non-linear effects
# f <- function(x) {
#   10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - .5)^2 + 10 * x[, 4] + 5 * x[, 5]
# }

# Function to generate simulated multi-omics data
gen_simmba <- function(nsample, # Sample size
                       snr = 1, # Signal to noise ratio
                       p.train = 0.7, # Train-test split ratio
                       de.prob = rep(0.1, 3), # DE probability across all modalities (vector)
                       de.downProb = rep(0.5, 3), # Down-regulation probability (vector)
                       de.facLoc = rep(1, 3), # DE factor location (vector)
                       de.facScale = rep(0.4, 3), # DE factor scale (vector)
                       ygen.mode = 'Friedman', # Y generation mode
                       nrep = 100, # Number of repetitions
                       seed = 1234) { # Random seed

  set.seed(seed)  # Set random seed

  # Initialize lists to store training and testing datasets
  trainDat <- testDat <- vector("list", nrep)
  names(trainDat) <- names(testDat) <- paste('Rep', 1: nrep, sep = '_')

  # Loop over the number of repetitions
  for (k in 1:nrep) {

    # Generate feature data using InterSIM
    pcl <- trigger_InterSIM(n = nsample)
    X <- pcl$feature_table %>% t() %>% as.matrix()

    # Initialize coefficients based on the number of features
    nfeature <- table(pcl$feature_metadata$featureType)

    # Generate realistic effect sizes using Splatter
    de.facs <- vector("list", 3)
    for (i in 1:3) {
      de.facs[[i]] <- splatter:::getLNormFactors(n.facs = nfeature[i],
                                                 sel.prob = de.prob[i],
                                                 neg.prob = de.downProb[i],
                                                 fac.loc = de.facLoc[i],
                                                 fac.scale = de.facScale[i])
    }

    # Convert to log fold changes (LFCs)
    beta0 <- log2(unlist(de.facs)) # Assuming LFC is log-normally distributed

    # Calculate residual variance sigma2 based on beta0 and snr
    mu <- as.matrix(X) %*% beta0
    sigma2 <- as.vector(var(mu) / snr)

    # Generate response variable Y using different models
    if (ygen.mode == 'LM') {
      Y = X %*% beta0 + rnorm(nsample) * sqrt(sigma2)
    } else if (ygen.mode == 'Friedman') {
      nonzero_index <- which(beta0 == 0)
      friedman_index <- sample(nonzero_index, 5)
      X.friedman <- X[, friedman_index]
      Xbeta.friedman = f(X.friedman)
      Y = Xbeta.friedman + rnorm(nsample) * sqrt(sigma2)
    } else if (ygen.mode == 'Friedman2') {
      nonzero_index <- which(beta0 == 0)
      friedman_index <- sample(nonzero_index, 5)
      X.friedman <- X[, friedman_index]
      Xbeta.friedman = f(X.friedman)
      Y = X %*% beta0 + Xbeta.friedman + rnorm(nsample) * sqrt(sigma2)
    } else {
      stop('Use `LM` for linear effects and `Friedman` or ``Friedman2`` for non-linear effects')
    }

    # Insert Y into the simulated datasets
    pcl$sample_metadata$Y <- as.vector(Y)

    # Split into training and testing sets
    train <- test <- pcl
    tr.row <- sample(1L:nsample, round(nsample * p.train), replace = FALSE)
    train$sample_metadata <- pcl$sample_metadata[tr.row, , drop = FALSE]
    test$sample_metadata <- pcl$sample_metadata[-tr.row, , drop = FALSE]
    train$feature_table <- pcl$feature_table[, tr.row, drop = FALSE]
    test$feature_table <- pcl$feature_table[, -tr.row, drop = FALSE]

    # Store train and test data
    trainDat[[k]] <- train
    testDat[[k]] <- test
  }

  # Return synthetic data and simulation parameters
  return(list(trainDat = trainDat,
              testDat = testDat,
              snr = snr,
              p.train = p.train,
              de.prob = de.prob,
              de.downProb = de.downProb,
              de.facLoc = de.facLoc,
              de.facScale = de.facScale,
              nrep = nrep,
              seed = seed))
}


# all(colnames(DD$trainDat$Rep_1$feature_table) == rownames(DD$trainDat$Rep_1$sample_metadata))
# all(colnames(DD$trainDat$Rep_2$feature_table) == rownames(DD$trainDat$Rep_2$sample_metadata))
# all(rownames(DD$trainDat$Rep_1$feature_table) == rownames(DD$trainDat$Rep_2$feature_table))



