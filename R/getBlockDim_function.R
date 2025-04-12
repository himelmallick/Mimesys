#' Generate Block Dimensions Based on Dirichlet Distribution
#'
#' This function generates block dimensions based on the Dirichlet distribution, ensuring no block has zero size.
#'
#' @param K Integer. The number of blocks.
#' @param p Integer. The total number of elements to be distributed among the blocks.
#' @param alpha Numeric. The concentration parameter for the Dirichlet distribution. Default is 2.
#' @return Integer vector. A sorted vector of block sizes in descending order.
#' @export
#' @examples
#' # Generate block dimensions for 5 blocks and 100 elements
#' getBlocksDim(K = 5, p = 100)
getBlocksDim <- function(K, p, alpha = 2) {
  nK = rep(0, K)
  while (min(nK) == 0) {
    pK <- LaplacesDemon::rdirichlet(1, rep(alpha, K))
    nK = pmax(1, floor(p * pK))
    nK = c(nK[1:K - 1], p - sum(nK[1:K - 1]))
  }
  return(sort(nK, decreasing = TRUE))
}
