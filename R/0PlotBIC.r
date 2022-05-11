#' Plot for model selection 
#'
#'
#'
#'
#'
#'
#' @param x  Output from the main function step1() (must be of class lmfa_modelselection).
#' @param ... Further arguments for the default S3 plot method.
#' @examples
#' \dontrun{
#' plot(modelselection)
#' }
#' @export


plot.lmfa_modelselection <- function(x, ...){
        if(missing(x)) stop("argument data is missing, with no default")
        
          modelcomparison <-matrix(NA,nrow=length(x),ncol=6)
          colnames(modelcomparison) <- c("LL","BIC","convergence","n_par","n_fact","n_state")
          modelnames <- NULL
          for(i in 1:length(x)){
          modelcomparison[i,1] <-x[[i]]$LL
          modelcomparison[i,2] <-x[[i]]$BIC
          modelcomparison[i,3] <-x[[i]]$convergence
          modelcomparison[i,4] <-x[[i]]$n_par
          modelcomparison[i,5] <-sum(x[[i]]$n_fact)
          modelcomparison[i,6] <-x[[i]]$n_state
          modelnames <- c(modelnames,paste("[",paste(x[[i]]$n_fact,collapse = ""),"]",sep=""))
          }
          rownames(modelcomparison) <- modelnames

    modelcomparison <- modelcomparison[modelcomparison[,"convergence"]==1,,drop=FALSE]
    
    if(nrow(modelcomparison)>1){
      modelcomparison <-modelcomparison[order(modelcomparison[,"BIC"]),]
    }
   


        #par(mfrow = c(2, 1),mar = c(4.5,4.5,2,2))
        
 

        #-------------------------------------------------------
        #                 LL
        #-------------------------------------------------------
        

        plot(modelcomparison[,"n_par"],modelcomparison[,"LL"],xlab = "n_par",
           ylab = "LL", ylim = c(min(modelcomparison[,"LL"])-1500,max(modelcomparison[,"LL"])+1500),
           col= c(rep("black",nrow(modelcomparison)) ), ...)
        text(modelcomparison[,"n_par"],modelcomparison[,"LL"],  rownames(modelcomparison),
           cex=1,pos=3)

        #-------------------------------------------------------
        #                 BIC
        #-------------------------------------------------------

        plot(modelcomparison[,"n_par"],modelcomparison[,"BIC"],xlab = "n_par",
             ylab = "BIC",ylim = c(min(modelcomparison[,"BIC"])-1500,max(modelcomparison[,"BIC"])+1500),col = c("red", rep("black",(nrow(modelcomparison)-1))),...)
        text(modelcomparison[,"n_par"],modelcomparison[,"BIC"],  rownames(modelcomparison),
             cex=1,pos=3)

       
        cat("\n")
        cat(paste('Use the navigation arrow in the "Plots" tab to switch',"\n"))
        cat(paste("between the LL and the BIC plots.","\n"))
        cat("\n")
        
#on.exit(par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0, mfrow = c(1,1)))
}
