#' Summary for \code{lmfa_chull}
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from \code{lmfa_chull}.
#' @param ... Further arguments for the default S3 summary method.
#' @examples
#' \dontrun{
#' summary(chullresults)
#' }
#' @export
summary.lmfa_chull <- function(object,...){

fitCHull <- object$fitCHull
Solution <- object$solution
Hull <- object$chull 
sumComplexModels <- object$sumComplexModels
  
        plot(fitCHull,col=c("black", "black","red"),pch=21, 
           bg="white",ylab="LL",xlab="n_par",...)

        
        cat("\n")
        cat(paste("Models on the upper boundary of the CHull:"),"\n")
        cat("\n")
        print(Hull)

        cat("\n")
        cat(paste("Selected model(s):"),"\n")
        cat("\n")
        print(Solution)

        
#calculate numerators
if(nrow(Hull)>2){

if(sumComplexModels>0){
  cat("\n")
  cat(paste("Note 1: The least and most complex models cannot be selected."),"\n")
  cat(paste("  Therefore, it is advisable to also visually inspect the CHull plot."),"\n")
        
  cat("\n")
  cat(paste("Note 2: The st value(s) of the",sumComplexModels,"best model(s) might"),"\n")
  cat(paste("  be artificially inflated. Therefore, it is advisable to also consider the"),"\n")
  cat(paste("  first less complex model."),"\n")
}else{
    cat("\n")
    cat(paste("Note: The least and most complex models cannot be selected."),"\n")
    cat(paste("  Therefore, it is advisable to also visually inspect the CHull plot."),"\n")
        
}
}else{
    cat("\n")
    cat(paste("Note: The least and most complex models cannot be selected."),"\n")
    cat(paste("  Therefore, it is advisable to also visually inspect the CHull plot."),"\n")
        
}


}


