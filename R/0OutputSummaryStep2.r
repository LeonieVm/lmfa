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
  cat(paste("R2_entropy: ",round(object$R2_entropy,2),sep=" = "),"\n")
  cat("\n")
  cat(paste("Classification errors:"),"\n")
  cat("\n")
  print(round(object$classification_errors,2))
  cat("\n")
  cat(paste("Classification-error probabilities:"),"\n")
  cat("\n")
  print(round(object$classification_errors_prob,2))
  cat("\n")
  cat(paste("State proportions:"),"\n")
  cat("\n")
  print(round(object$state_proportions,2),row.names = FALSE)
  cat("\n")

  
}