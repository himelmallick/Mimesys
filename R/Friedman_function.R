#' Generate Non-Linear Effects Using Friedman Function
#'
#' This function generates non-linear effects based on the Friedman function.
#'
#' @param x Matrix. A matrix where only the first 5 columns are used to generate non-linear effects.
#' @return Numeric vector. The generated non-linear effects.
#' @export
#' @examples
#' x <- matrix(runif(50), ncol = 5)
#' y <- f(x)
Friedman_function <- function(x) {
                      10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - .5)^2 + 10 * x[, 4] + 5 * x[, 5]
                                  }
