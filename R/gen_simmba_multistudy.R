#' Multi-study wrapper around gen_simmba()
#'
#' @param nsample_vec Integer vector, length = n_study, samples per study
#' @inheritParams gen_simmba
#' @return list(trainDat, testDat, ...) where
#'   trainDat[[k]][[s]] is study s in replication k
gen_simmba_multistudy <- function(nsample_vec,
                                  snr = 1,
                                  p.train = 0.7,
                                  de.prob = rep(0.1, 3),
                                  de.downProb = rep(0.5, 3),
                                  de.facLoc = rep(1, 3),
                                  de.facScale = rep(0.4, 3),
                                  ygen.mode = "Friedman",
                                  nrep = 100,
                                  seed = 1234) {
  
  set.seed(seed)
  n_study <- length(nsample_vec)
  
  # outer lists: over replications
  trainDat <- testDat <- vector("list", nrep)
  names(trainDat) <- names(testDat) <- paste0("Rep_", seq_len(nrep))
  
  for (k in seq_len(nrep)) {
    # inner lists: over studies
    train_list <- vector("list", n_study)
    test_list  <- vector("list", n_study)
    names(train_list) <- names(test_list) <- paste0("Study_", seq_len(n_study))
    
    for (s in seq_len(n_study)) {
      # jitter seed so each (rep, study) combo is different
      this_seed <- seed + 1000 * k + 10 * s
      
      sim_s <- gen_simmba(
        nsample    = nsample_vec[s],
        snr        = snr,
        p.train    = p.train,
        de.prob    = de.prob,
        de.downProb = de.downProb,
        de.facLoc  = de.facLoc,
        de.facScale = de.facScale,
        ygen.mode  = ygen.mode,
        nrep       = 1,
        seed       = this_seed
      )
      
      # extract the only replication
      tr <- sim_s$trainDat[[1]]
      te <- sim_s$testDat[[1]]
      
      # add study label
      tr$sample_metadata$study <- paste0("Study_", s)
      te$sample_metadata$study <- paste0("Study_", s)
      
      train_list[[s]] <- tr
      test_list[[s]]  <- te
    }
    
    trainDat[[k]] <- train_list
    testDat[[k]]  <- test_list
  }
  
  list(
    trainDat   = trainDat,
    testDat    = testDat,
    snr        = snr,
    p.train    = p.train,
    de.prob    = de.prob,
    de.downProb = de.downProb,
    de.facLoc  = de.facLoc,
    de.facScale = de.facScale,
    nrep       = nrep,
    seed       = seed,
    nsample_vec = nsample_vec
  )
}