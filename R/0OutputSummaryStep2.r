#' Summary statistics for step2()
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from the main function step2()
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#' @param ... Further arguments for the default S3 summary method
#' @examples
#' \dontrun{
#' summary(results2)
#' }
#' @export
summary.lmfa_step2 <- function(object, rounding = 2,...){
  if(!is.numeric(rounding)) stop("rounding must be a single scalar")
  if(length(rounding)>1) stop("rounding must be a single scalar")
  n_state <- length(object$state_proportions)
  cat("\n")
  cat(paste("R2_entropy: ",round(object$R2_entropy,rounding),sep=" = "),"\n")
  cat("\n")
  cat(paste("Classification errors:"),"\n")
  cat("\n")
  classification_errors <- round(object$classification_errors,rounding)
    colnames(classification_errors) <-paste("S",1:n_state,sep="")
    rownames(classification_errors) <-paste("S",1:n_state,sep="")
  print(classification_errors)
  cat("\n")
  cat(paste("Classification-error probabilities:"),"\n")
  cat("\n")
    classification_errors_prob <- round(object$classification_errors_prob,rounding)
    colnames(classification_errors_prob) <-paste("S",1:n_state,sep="")
    rownames(classification_errors_prob) <-paste("S",1:n_state,sep="")
  print(classification_errors_prob)
  cat("\n")
  cat(paste("State proportions:"),"\n")
  cat("\n")
  state_proportions <- c(round(object$state_proportions,rounding))
  names(state_proportions) <- paste("S",1:n_state,sep="")
  print(state_proportions)
  cat("\n")

  
}