#' Create Venn Diagrams for Multiview Data
#'
#' This function generates Venn diagrams based on the number of views.
#'
#' @param M Number of views.
#' @param areas A vector of areas for each set.
#' @param n_intersections A matrix of intersections.
#' @param category A vector of category labels.
#' @param mainDir Main directory for saving plots.
#' @param subDir Subdirectory for saving plots.
#' @import VennDiagram
#' @export
#'
create_venn_diagram <- function(M, areas, n_intersections, category, mainDir, subDir) {
  if (M == 2) {
    draw.pairwise.venn(
      area1 = areas[1], area2 = areas[2], cross.area = n_intersections[2, 1],
      category = category[1:2], col = "black", fill = c('#332288', '#882255'),
      alpha = c(0.4, 0.4), scaled = TRUE, lty = 'blank', cex = 1.5, cat.cex = 1.6
    )
  } else if (M == 3) {
    draw.triple.venn(
      area1 = areas[1], area2 = areas[2], area3 = areas[3],
      n12 = n_intersections[2, 1], n13 = n_intersections[2, 2], n23 = n_intersections[2, 3],
      n123 = n_intersections[3, 1],
      category = category[1:3], col = "black", fill = c('#332288', '#882255', '#CC6677'),
      alpha = c(0.4, 0.4, 0.35), scaled = TRUE, lty = 'blank', cex = 1.5, cat.cex = 1.6
    )
  } else if (M == 4) {
    draw.quad.venn(
      area1 = areas[1], area2 = areas[2], area3 = areas[3], area4 = areas[4],
      n12 = n_intersections[2, 1], n13 = n_intersections[2, 2], n14 = n_intersections[2, 3],
      n23 = n_intersections[2, 4], n24 = n_intersections[2, 5], n34 = n_intersections[2, 6],
      n123 = n_intersections[3, 1], n124 = n_intersections[3, 2], n134 = n_intersections[3, 3],
      n234 = n_intersections[3, 4], n1234 = n_intersections[4, 1],
      category = category[1:4], col = "black", fill = c('#332288', '#882255', '#CC6677', '#117733'),
      alpha = c(0.4, 0.4, 0.35, 0.35), scaled = TRUE, lty = 'blank', cex = 1.5, cat.cex = 1.6
    )
  } else if (M == 5) {
    draw.quintuple.venn(
      area1 = areas[1], area2 = areas[2], area3 = areas[3], area4 = areas[4], area5 = areas[5],
      n12 = n_intersections[2, 1], n13 = n_intersections[2, 2], n14 = n_intersections[2, 3], n15 = n_intersections[2, 4],
      n23 = n_intersections[2, 5], n24 = n_intersections[2, 6], n25 = n_intersections[2, 7],
      n34 = n_intersections[2, 8], n35 = n_intersections[2, 9], n45 = n_intersections[2, 10],
      n123 = n_intersections[3, 1], n124 = n_intersections[3, 2], n125 = n_intersections[3, 3],
      n134 = n_intersections[3, 4], n135 = n_intersections[3, 5], n145 = n_intersections[3, 6],
      n234 = n_intersections[3, 7], n235 = n_intersections[3, 8], n245 = n_intersections[3, 9],
      n345 = n_intersections[3, 10], n1234 = n_intersections[4, 1], n1235 = n_intersections[4, 2],
      n1245 = n_intersections[4, 3], n1345 = n_intersections[4, 4], n2345 = n_intersections[4, 5],
      n12345 = n_intersections[5, 1],
      category = category[1:5], col = "black", fill = c('#332288', '#882255', '#CC6677', '#117733', '#44AA99'),
      alpha = c(0.4, 0.4, 0.35, 0.35, 0.3), scaled = TRUE, lty = 'blank', cex = 1.5, cat.cex = 1.6
    )
  } else {
    stop("Venn diagrams for more than 5 sets are not supported with this function.")
  }
}
