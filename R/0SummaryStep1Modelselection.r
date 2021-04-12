#' Summary statistics for \code{step1} with model selection
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from \code{step1} with model selection.
#' @param ... Further arguments for the default S3 summary method.
#' @examples
#' \dontrun{
#' summary(results1)
#' }
#' @export



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