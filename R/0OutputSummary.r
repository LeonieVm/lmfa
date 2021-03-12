#' Summary
#'
#'
#'
#'
#'
#'
#' @param object The dataset (must be a dataframe and contain complete cases only).
#' @param ... The dataset (must be a dataframe and contain complete cases only).
#' @return Returns the measurement model parameters, the proportional and
#' modal state assignments, and the classification errors.
#'
#' @examples
#' \dontrun{
#' fitStep1Step2 <- Step1Step2(input_file,variable_columns,n_state,
#'                       n_fact,n_starts=25,n_initial_ite=15,n_m_step=10,
#'                       em_tolerance=1e-6,m_step_tolerance=1e-3,max_iterations=500)
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