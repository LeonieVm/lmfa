#' Summary statistics for step3()
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from the main function Step3()
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#' @param ... Further arguments for the default S3 summary method
#' @examples
#' \dontrun{
#' summary(results3)
#' }
#' @export
summary.lmfa_step3 <- function(object, rounding = 4, ...){
  if(!is.numeric(rounding)) stop("rounding must be a single scalar")
  if(length(rounding)>1) stop("rounding must be a single scalar")
  if(object$convergence==1){
    cat("\n")
    cat(paste("Estimation converged after",round(object$seconds,rounding),"seconds.","\n"))
  }else{
    cat("\n")
    cat(paste("Maximum number of iterations reached without convergence.","\n"))
    cat("\n")
  }
  cat("\n")
  cat(paste("LL",round(object$LL,rounding),sep=" = "),"\n")
  cat("\n")

  cat("\n")
  print(round(object$Wald_tests,rounding))

  cat("\n")
  cat(paste("Information about the regression parameter estimates:"),"\n")
  cat("\n")
  cat(paste("For the initial state parameters, state 1 is the reference "),"\n")
  cat(paste("category. The transition intensity parameters are sorted "),"\n")
  cat(paste("by rows of the transition matrix and the staying rates serve "),"\n")
  cat(paste("as references."),"\n")
  
  cat("\n")
  print(round(object$estimates,rounding))
  cat("\n")
  
  cat(paste("Probabilities for a unit interval and covariate scores"),"\n")
  cat(paste("equal to the sample means:"),"\n")
  probabilities(object)

  cat("\n")
  cat(paste("Use the function probabilities() to calculate initial state "),"\n")
  cat(paste("and transition probabilities for any covariate score (and interval) "),"\n")
  cat(paste("of interest."),"\n")

}