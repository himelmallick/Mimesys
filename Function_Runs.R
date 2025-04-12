install.packages("roxygen2")

# Document your package
devtools::document()

library(MultiviewSimulator)



### Input for Functions

Benchmarking_output <- Benchmarking_multiomics(

  mainDir = '/Users/siddheshkulkarni/Desktop',

  subDir = 'plots_simulated_data',

  seed.number = 1234,

  M = 3, # Number of views (NB: the inputs are currently set up only for this specific case!!!)
  p_m = c(500, 700, 1000),  # Number of features per view
  n = 100,  # Train set size
  nTest = 100,  # Test set size
  fixed_act_patter = TRUE,  # Set a custom activity pattern in the joint component

  J0_m = c(5, 6, 7),  # Number of groups in joint component
  K0_m = c(5, 7, 10),  # Number of factors in view-specific components
  S0_m = c(6, 7, 8),  # Number of groups in view-specific components

  nu2_0 = 0.1,  # Input the joint loadings variability
  pi2_0 = 0.1,  # Input the view-specific loadings variability

  decaying_loadings = FALSE,  # Enforce loadings decay (re-weighting importance of latent axis of variation)
  decaying_type = 'sqrt',  # One of 'sqrt' and 'log'

  pJ_zero = 0.6,  # Probability of group-wise zero in joint component
  pJ_sign = 0.4,  # Probability of group-wise sign-switch in joint component
  pJ_zeroEntr = 0.5,  # Probability of entry-wise zero in joint component

  pS_zero = 0.6,  # Probability of group-wise zero in specific components
  pS_sign = 0.4,  # Probability of group-wise sign-switch in specific components
  pS_zeroEntr = 0.5, # Probability of entry-wise zero in specific components

  ### Additional Matrix Plots ##

  Activity_pattern_plot = TRUE,  #Save the Activity Pattern Plot

  Emprical_correlation_plot = TRUE, #Save the Empirical correlation plot

  Crossprd_plot=TRUE, ## crossporduct matrices

  Data_save=TRUE  ## save the data in directory

)

