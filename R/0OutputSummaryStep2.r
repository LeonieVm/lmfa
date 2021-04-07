#' Summary statistics for step2()
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from the main function step2()
#' @param ... Further arguments for the default S3 summary method
#' @examples
#' \dontrun{
#' summary(results2)
#' }
#' @export
summary.lmfa_step2 <- function(object, ...){

  cat("\n")
  cat(paste("R2_entropy: ",object$R2_entropy,sep=" = "),"\n")
  cat("\n")
  cat(paste("Classification errors"),"\n")
  cat("\n")
  print(object$classification_errors)
  cat("\n")
  cat(paste("Classification-error probabilities"),"\n")
  cat("\n")
  print(object$classification_errors_prob)
  cat("\n")
  cat(paste("state proportions"),"\n")
  cat("\n")
  print(object$pi_k)
  cat("\n")

  
}