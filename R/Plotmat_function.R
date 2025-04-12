#' Plot a Matrix with Specified Settings
#'
#' This function plots a matrix with specified settings, creating an image file.
#'
#' @param mat_m Matrix. The matrix to be plotted.
#' @param n_colors Integer. The number of colors to use in the plot. Default is 256.
#' @param lab_col Character. Optional label for the color legend. Default is NULL.
#' @param is_cor Logical. If TRUE, the matrix is treated as a correlation matrix with values ranging from -1 to 1. Default is FALSE.
#' @param folder Character. The folder to save the plot file. Default is '~/Desktop/'.
#' @param filename Character. The filename for the plot file (without extension). Default is 'prova'.
#' @return None. The function saves the plot as a PNG file.
#' @export
#' @examples
#' # Example matrix
#' mat <- matrix(runif(100), nrow = 10)
#'
#' # Plot the matrix
#' plotmat(mat, n_colors = 256, lab_col = "Intensity", is_cor = FALSE, folder = tempdir(), filename = "matrix_plot")
plotmat <- function(mat_m, n_colors = 256, lab_col = NULL, is_cor = FALSE,
                    folder = '~/Desktop/', filename = 'prova') {

  if (is_cor) {
    zL = -1
    zR = 1
  } else {
    zM = max(abs(mat_m))
    zL = -zM
    zR = zM
    mat_m = t(mat_m)
  }

  file_path <- file.path(folder, paste(filename, '.png', sep = ""))
  png(file = file_path, height = 5, width = 5.6, res = 300, pointsize = 3.5, unit = 'cm')

  colors <- colorRampPalette(c("#AA4499", "white", "#117733"))(n_colors)
  ticks <- sapply(c(0.7 * zL, 0, 0.7 * zR), round, 2)

  if (is_cor) {
    par(pty = "s", mar = c(1, 0, 1, 4))
    image(mat_m, useRaster = TRUE, asp = 1, axes = FALSE, col = colors,
          xlim = c(0, 1), ylim = c(0, 1), zlim = c(zL, zR))
    rect(0, 0, 1, 1, border = "black")
  } else {
    par(pty = "s", mar = c(0, 0, 0, 4))
    image(mat_m, useRaster = TRUE, asp = 1, axes = FALSE, col = colors,
          xlim = c(0 - 0.1 * 15 / nrow(mat_m), 1 + 0.1 * 15 / nrow(mat_m)), ylim = c(0, 1), zlim = c(zL, zR))
    rect(0 - 0.04 * 15 / nrow(mat_m), 0, 1 + 0.04 * 15 / nrow(mat_m), 1, border = "black")
  }

  image.plot(zlim = c(zL, zR), col = colors, legend.only = TRUE, side = 4,
             axis.args = list(at = ticks, labels = TRUE, cex.axis = 0.9), legend.shrink = 0.5,
             legend.width = 0.9, legend.mar = 4.5)

  if (!is.null(lab_col)) {
    mtext(lab_col, side = 4, line = 1, cex = 1.3, las = 1, at = 0.8)
  }

  dev.off()
}
