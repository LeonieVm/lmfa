#' Summary statistics for step3()
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from the main function Step3()
#' @param ... Further arguments for the default S3 summary method
#' @examples
#' \dontrun{
#' summary(results3)
#' }
#' @export
summary.lmfa_step3 <- function(object, ...){
  if(object$convergence==1){
    cat("\n")
    cat(paste("Estimation converged after",round(object$seconds,2),"seconds.","\n"))
  }else{
    cat("\n")
    cat(paste("Maximum number of iterations reached without convergence.","\n"))
    cat("\n")
  }
  cat("\n")
  cat(paste("LL",round(object$LL,2),sep=" = "),"\n")
  cat("\n")

  cat("\n")
  print(round(object$WaldTests,2))

  cat("\n")
  cat(paste("Information about the regression parameter estimates:"),"\n")
  cat("\n")
  cat(paste("For the initial state parameters, state 1 is the reference "),"\n")
  cat(paste("category. The transition intensity parameters are sorted "),"\n")
  cat(paste("by rows of the transition matrix and the staying rates serve "),"\n")
  cat(paste("as references."),"\n")
  cat("\n")
  cat(paste("Use the function probabilities() to calculate initial state "),"\n")
  cat(paste("and transition probabilities for any covariate score (and interval) "),"\n")
  cat(paste("of interest."),"\n")

  cat("\n")
  print(round(object$estimates,2))
  cat("\n")
}