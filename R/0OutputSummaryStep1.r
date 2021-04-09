#' Summary statistics for step1()
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from the main function Step1Step2()
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#' @param ... Further arguments for the default S3 summary method
#'
#' @examples
#' \dontrun{
#' summary(results1)
#' }
#' @export
summary.lmfa_step1 <- function(object, rounding = 2,...){
    if(!is.numeric(rounding)) stop("rounding must be a single scalar")
    if(length(rounding)>1) stop("rounding must be a single scalar")

  if(object$convergence==1){
    cat("\n")
    cat(paste("Estimation converged after",round(object$seconds,2),"seconds and",object$n_it,"iterations.","\n"))
  }else{
    cat("\n")
    cat(paste("Maximum number of iterations reached without convergence.","\n"))
    cat("\n")
  }
  cat("\n")
  cat(paste("LL",round(object$LL,rounding),sep=" = "),"\n")
  cat("\n")
  cat(paste("Number of states",object$n_state,sep=": "),"\n")
  cat("\n")
  cat(paste("Number of factors: [", paste(object$n_fact, collapse = " "), "]", sep = ""),"\n")
  cat("\n")
  cat("\n")
  cat(paste("Intercepts:"),"\n")
  cat("\n")
  print(round(object$intercepts,rounding))
  cat("\n")
  cat(paste("Obliquely rotated within-state standardized loadings:"),"\n")
  cat("\n")
  print(round(object$loadings_w_obli,rounding))
  cat("\n")
  cat(paste("Obliquely rotated between-state standardized loadings:"),"\n")
  cat("\n")
  print(round(object$loadings_b_obli,rounding))
  cat("\n")
  cat(paste("Factor correlations after oblique rotation:"),"\n")
  cat("\n")
  for(i in 1:object$n_state){
    cat(paste("S",i,sep=""),"\n")
    print(lapply(object$factor_correlations_obli_list,round,rounding)[[i]])
    cat("\n")
  }
  cat("\n")
  cat(paste("Unique variances:"),"\n")
  cat("\n")
  print(round(object$unique_variances,rounding))
  cat("\n")
}

