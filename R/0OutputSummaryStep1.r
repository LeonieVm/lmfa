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
}

summary.lmfa_modelselection <- function(object, ...){
        modelcomparison <-matrix(NA,nrow=length(object),ncol=4)
        colnames(modelcomparison) <- c("LL","BIC","convergence","n_par")
        modelnames <- NULL
        for(i in 1:length(object)){
                modelcomparison[i,1] <-object[[i]]$LL
                modelcomparison[i,2] <-object[[i]]$BIC
                modelcomparison[i,3] <-object[[i]]$convergence
                modelcomparison[i,4] <-object[[i]]$n_par
                modelnames <- c(modelnames,paste("[",paste(object[[i]]$n_fact,collapse = ""),"]",sep=""))
        }
        rownames(modelcomparison) <- modelnames
        modelcomparison <-modelcomparison[order(modelcomparison[,"BIC"]),]
        modelcomparison
}