#' A Generic Convex-Hull-Based Model Selection Method
#'
#' This function is based on the CHull function from the R package multichull
#'
#'
#' @param modelselection Output from the main function step1().
#' @param PercentageFit Required proportion of increase in fit of a more complex model.
#' @param ... Further arguments for the plot function.

#'
#'
#' @return Returns CHull output.
#'
#' @examples
#' \dontrun{
#' chull_lmfa(modelselection)
#' }


chull_lmfa <- function(modelselection,PercentageFit = 0.01,...){
    for(i in 1:length(PercentageFit)){
        CHullInput <-summary(modelselection)
        CHullInput <- CHullInput[CHullInput[,"convergence"]==1,]
        fitCHull <- CHull(CHullInput[,c("n_par","LL")],bound = "upper", PercentageFit = PercentageFit)
        plot(fitCHull,col=c("black", "black","red"),pch=21, 
           bg="white",ylab="LL",xlab="n_par",...)
    return(fitCHull)
    }
    
}
