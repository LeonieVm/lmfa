#' A Generic Convex-Hull-Based Model Selection Method
#'
#' This function is based on the CHull function from the R package multichull
#'
#'
#' @param x Output from the main function step1() (must be of class lmfa_modelselection).
#' @param PercentageFit Required proportion of increase in fit of a more complex model.
#' @param ... Further arguments for the CHull plot function.

#'
#'
#' @return Returns CHull output.
#'
#' @examples
#' \dontrun{
#' chull_lmfa(x)
#' }


chull_lmfa <- function(x,PercentageFit = 0.01,...){

    if(missing(x)) stop("argument data is missing, with no default")
    if(class(x)!="lmfa_modelselection") stop("x must be of class lmfa_modelselection")
    if(!is.numeric(PercentageFit)) stop("PercentageFit must be a single scalar")
    if(length(PercentageFit)>1) stop("PercentageFit must be a single scalar")

    for(i in 1:length(PercentageFit)){
        CHullInput <-summary(x)
        CHullInput <- CHullInput[CHullInput[,"convergence"]==1,]
        fitCHull <- CHull(CHullInput[,c("n_par","LL")],bound = "upper", PercentageFit = PercentageFit)
        plot(fitCHull,col=c("black", "black","red"),pch=21, 
           bg="white",ylab="LL",xlab="n_par",...)
    return(fitCHull)
    }
    
}
