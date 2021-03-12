#' Summary statistics for Step1Step2 analysis
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from the main function Step1Step2()
#' @param ... Further arguments for the default S3 summary method
#' @examples
#' \dontrun{
#' summary(results1)
#' }
#' @export
summary.lmfa <- function(object, ...){
  if(object$convergence==1){
    cat("\n")
    cat(paste("Estimation converged after",object$n_it,"iterations.","\n"))
  }else{
    cat("\n")
    cat(paste("Maximum number of iterations reached without convergence.","\n"))
    cat("\n")
  }
  cat("\n")
  cat(paste("LL",round(object$LL,4),sep=" = "),"\n")
  cat("\n")
  cat(paste("Number of states",object$n_state,sep=": "),"\n")
  cat("\n")
  cat(paste("Number of factors: [", paste(results1a$n_fact, collapse = " "), "]", sep = ""),"\n")
  cat("\n")
}