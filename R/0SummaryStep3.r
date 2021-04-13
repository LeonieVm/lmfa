#' Summary statistics for \code{step3}
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from \code{step3}.
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#' @param ... Further arguments for the default S3 summary method.
#' @examples
#' \dontrun{
#' summary(transitionmodel)
#' }
#' @export
summary.lmfa_step3 <- function(object, rounding = 4, ...){
  if(!is.numeric(rounding)) stop("rounding must be a single scalar")
  if(length(rounding)>1) stop("rounding must be a single scalar")
  n_state <- object$n_state
  cat("\n")
  cat(paste("Model estimation:"),"\n")
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
  cat("-------------------------------------------------------------")
  cat("\n")
  cat("\n")
  cat(paste("Wald tests:"),"\n")
  cat("\n")
  print(round(object$Wald_tests,rounding))

  cat("\n")
  cat("-------------------------------------------------------------")
  cat("\n")
  cat("\n")
  
  cat(paste("Parameter estimates:"),"\n")
  cat("\n")
  print(round(object$estimates,rounding))
  cat("\n")
  cat(paste("Note: For the initial state parameters, state 1 is the"),"\n")
  cat(paste("reference category. The transition intensity parameters"),"\n")
  cat(paste("are sorted by rows of the transition matrix and the"),"\n")
  cat(paste("staying rates serve as references."),"\n")
  
  cat("\n")
  cat("-------------------------------------------------------------")
  cat("\n")
  cat("\n")
  cat(paste("Probabilities:"),"\n")
  cat("\n")
 
  probabilities(object, rounding = rounding)

  cat("\n")
  cat(paste("Note: The probabilities are calculated for covariate scores"),"\n")
  cat(paste("equal to the sample means (and a unit time interval). Use"),"\n")
  cat(paste("the function probabilities() to calculate initial state and"),"\n")
  cat(paste("transition probabilities for any covariate score (and "),"\n")
  cat(paste("interval) of interest."),"\n")

  cat("\n")
  cat("-------------------------------------------------------------")
  cat("\n")
  cat("\n")
  cat(paste("State proportions:"),"\n")
  cat("\n")
  state_proportions <- c(round(object$state_proportions,rounding))
  names(state_proportions) <- paste("S",1:n_state,sep="")
  print(state_proportions)
  cat("\n")


}