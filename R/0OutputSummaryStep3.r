#' Summary statistics for Step3 analysis
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
#' summary(results1)
#' }
#' @export
summary.lmfa_step3 <- function(object, ...){
  if(object$convergence==1){
    cat("\n")
    cat(paste("Estimation converged.","\n"))
  }else{
    cat("\n")
    cat(paste("Maximum number of iterations reached without convergence.","\n"))
    cat("\n")
  }
  cat("\n")
  cat(paste("LL",round(object$LL,4),sep=" = "),"\n")
  cat("\n")

  cat("\n")
  cat(paste("State 1 is the reference category for the initial state probabilities"),"\n")
  cat("\n")
  cat("\n")
  cat(paste("Transition intensities are sorted by rows of the transition matrix."),"\n")
  cat(paste("The staying rates serve as baseline"),"\n")
  cat("\n")
  
  cat("\n")
  print(object$estimates)
  print(object$WaldTests)
  cat("\n")
}