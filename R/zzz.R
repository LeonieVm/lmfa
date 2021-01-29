#' `lmfa`: The R-package for exploring within-person changes and between person differences in measurement models in (intensive) longitudinal data.
#'
#' @details
#' The package employs a three-step estimation of continuous-time latent Markov factor analysis.
#' 
#' @section Section name:
#' You can add multiple sections.
#'
#' @import stats
#' @import mclust
#' @import foreach
#' @import doParallel
#' @import parallel
#' @import msm
#' @import expm
#' @docType package
#'
#' @name lmfa
"_PACKAGE"

.onAttach <- function(libname, pkgname) {
    # Print a welcome message when the package is attached.
    cat(rep("-", 47), "\n")
    cat("`lmfa`: The R-package for exploring within-person changes and between person differences in measurement models in (intensive) longitudinal data.", "\n")
    cat(rep("-", 47), "\n")
}