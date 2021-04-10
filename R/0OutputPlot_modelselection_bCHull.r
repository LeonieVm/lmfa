#' A Generic Convex-Hull-Based Model Selection Method
#'
#' This function is based on the CHull function from the R package multichull
#'
#'
#' @param x Model-selection output from the function step1() (must be of class lmfa_modelselection).
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


chull_lmfa <- function(x,PercentageFit = 0.001,...){

    if(missing(x)) stop("argument data is missing, with no default")
    if(class(x)!="lmfa_modelselection") stop("x must be of class lmfa_modelselection")
    if(!is.numeric(PercentageFit)) stop("PercentageFit must be a single scalar")
    if(length(PercentageFit)>1) stop("PercentageFit must be a single scalar")

    
        CHullInput <-summary(x)
        CHullInput <- CHullInput[CHullInput[,"convergence"]==1,]
        fitCHull <- CHull(CHullInput[,c("n_par","LL")],bound = "upper", PercentageFit = PercentageFit)
        plot(fitCHull,col=c("black", "black","red"),pch=21, 
           bg="white",ylab="LL",xlab="n_par",...)

        Hull <- fitCHull$Hull
        Solution <- fitCHull$Solution
        colnames(Hull) <- c("n_par","LL","st")
        colnames(Solution) <- c("n_par","LL")
        
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
store <- matrix(NA,nrow = (nrow(Hull)-2), ncol = 2)
for(i in 2:(nrow(Hull)-1)){
  st <- Hull[i,3]
  numerator <- ((Hull[i,2]-Hull[i-1,2])/(Hull[i,1]-Hull[i-1,1]))
  store[i-1,1] <- st
  store[i-1,2] <- numerator
}
sumComplexModels <- sum(abs(store[,1]-store[,2])<10)
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

#object <- (list(solution = Solution,
#            chull = Hull))
#            object
}
