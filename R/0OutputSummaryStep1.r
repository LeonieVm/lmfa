#' Summary statistics for step1()
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
summary.lmfa_step1 <- function(object, ...){
  if(object$convergence==1){
    cat("\n")
    cat(paste("Estimation converged after",round(object$seconds,2),"seconds and",object$n_it,"iterations.","\n"))
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
  cat(paste("Number of factors: [", paste(object$n_fact, collapse = " "), "]", sep = ""),"\n")
  cat("\n")
  cat("\n")
  cat(paste("Intercepts:"),"\n")
  cat("\n")
  print(object$intercept)
  cat("\n")
  cat(paste("Obliquely rotated within-state standardized loadings:"),"\n")
  cat("\n")
  print(object$loadings_w_obli)
  cat("\n")
  cat(paste("Obliquely rotated between-state standardized loadings:"),"\n")
  cat("\n")
  print(object$loadings_b_obli)
  cat("\n")
  cat(paste("Unique variances:"),"\n")
  cat("\n")
  print(object$unique_variances)
  cat("\n")
}

