#' Package initialization function
#'
#' Ensures all required packages are available when the package is loaded.
#' This helps prevent errors when the package is used in different R environments.
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  required_pkgs <- c(
    "tidyverse", "tmle", "SuperLearner", "glmnet", "gam", "foreach", 
    "nnls", "dplyr", "stringr", "tibble", "readr", "lubridate", "forcats",
    "tidyr", "purrr", "rlang", "magrittr", "stats", "utils"
  )
  
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    warning(paste0(
      "ReSliceTMLE: The following required packages are not installed: ",
      paste(missing_pkgs, collapse = ", "), 
      ". Please install them with: install.packages(c('",
      paste(missing_pkgs, collapse = "', '"), "'))"
    ))
  }
}

#' Package attachment function
#'
#' Ensures proper message and loading of namespace when the package is attached
#' with library() or require()
#'
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("ReSliceTMLE: Resampling-Based Targeted Maximum Likelihood Estimation")
}