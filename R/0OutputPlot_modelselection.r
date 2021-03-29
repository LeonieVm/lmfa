#' Plot for model selection 
#'
#'
#'
#'
#'
#'
#' @param x An x storing output from the main function Step3()
#' @param ... Further arguments for the default S3 summary method
#' @examples
#' \dontrun{
#' plot(results_modelselection)
#' }
#' @export


plot.lmfa_modelselection <- function(x, ...){
        modelcomparison <-matrix(NA,nrow=length(x),ncol=4)
        colnames(modelcomparison) <- c("LL","BIC","convergence","n_par")
        modelnames <- NULL
        for(i in 1:length(x)){
                modelcomparison[i,1] <-x[[i]]$LL
                modelcomparison[i,2] <-x[[i]]$BIC
                modelcomparison[i,3] <-x[[i]]$convergence
                modelcomparison[i,4] <-x[[i]]$n_par
                modelnames <- c(modelnames,paste("[",paste(x[[i]]$n_fact,collapse = ""),"]",sep=""))
        }
        rownames(modelcomparison) <- modelnames
        modelcomparison <-modelcomparison[order(modelcomparison[,"BIC"]),]
        
       
        par(mfrow = c(2, 1),mar = c(4.5,4.5,2,2))
        
        plot(modelcomparison[,"n_par"],modelcomparison[,"BIC"],xlab = "n_par",
             ylab = "BIC",ylim = c(min(modelcomparison[,"BIC"])-1000,max(modelcomparison[,"BIC"])+1000))
        text(modelcomparison[,"n_par"],modelcomparison[,"BIC"],  rownames(modelcomparison),
             cex=0.75,pos=3)
        plot(modelcomparison[,"n_par"],modelcomparison[,"LL"],xlab = "n_par",
             ylab = "LL", ylim = c(min(modelcomparison[,"LL"])-1000,max(modelcomparison[,"LL"])+1000))
        text(modelcomparison[,"n_par"],modelcomparison[,"LL"],  rownames(modelcomparison),
             cex=0.75,pos=3)
        
}
