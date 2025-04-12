#' Install and Load Necessary Libraries
#'
#' This function checks for and installs missing CRAN and Bioconductor libraries, and then loads all specified libraries.
#'
#' @param packages Character vector. A vector of CRAN package names to install and load.
#' @return None.
#' @export
#' @examples
#' # Define the required libraries
#' libraries <- c("VennDiagram", "grid", "caret", "fields", "ggplot2", "InterSIM", "tidyverse", "splatter", "NMF", "stringr")
#'
#' # Check for and install Bioconductor packages separately
#' bioc_packages <- c("Biobase", "splatter")
#' missing_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
#' if (length(missing_bioc_packages) > 0) {
#'   BiocManager::install(missing_bioc_packages)
#' }
#'
#' # Install and load CRAN packages
#' install_and_load(libraries)
install_and_load <- function(packages= c("VennDiagram",
                                         "grid",
                                         "caret",
                                         "fields",
                                         "ggplot2",
                                         "InterSIM",
                                         "tidyverse",
                                         "splatter",
                                         "NMF",
                                         "stringr",
                                         "r.jive",
                                         "psych",
                                         "gridExtra")) {
  # Check for missing packages
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]

  # Install missing packages
  if (length(missing_packages) > 0) {
    install.packages(missing_packages)
  }

  # Load all packages
  invisible(lapply(packages, library, character.only = TRUE))
}

# Check for and install Bioconductor packages separately

bioc_packages <- c("Biobase", "splatter")

missing_bioc_packages <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(missing_bioc_packages) > 0) {
  BiocManager::install(missing_bioc_packages)
}



