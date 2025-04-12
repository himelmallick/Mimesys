#' Heatmap for Correlation Matrix
#'
#' This function generates a heatmap for a given correlation matrix. The heatmap
#' uses a color gradient to represent values in the correlation matrix, with blue
#' for negative values, white for neutral (zero), and red for positive values.
#'
#' @param cor_matrix A square numeric matrix representing correlations between variables.
#' @param title A character string specifying the title of the heatmap. Default is "Heatmap".
#'
#' @return A ggplot2 object representing the heatmap of the correlation matrix.
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal theme element_text
#' @export
#' @examples
#' # Example usage:
#' cor_matrix <- cor(mtcars)
#' heatmap_cor_matrix(cor_matrix, title = "Correlation Heatmap for mtcars")
heatmap_cor_matrix <- function(cor_matrix, title = "Heatmap") {
  cor_data <- reshape2::melt(cor_matrix)
  ggplot2::ggplot(cor_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    labs(title = title, x = "Variables", y = "Variables") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
