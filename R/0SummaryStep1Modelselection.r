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
#' summary(modelselection)
#' }
#' @export



summary.lmfa_modelselection <- function(object, ...){
       modelcomparison <-matrix(NA,nrow=length(object),ncol=6)
colnames(modelcomparison) <- c("LL","BIC","convergence","n_par","n_fact","n_state")
modelnames <- NULL
for(i in 1:length(object)){
  modelcomparison[i,1] <-object[[i]]$LL
  modelcomparison[i,2] <-object[[i]]$BIC
  modelcomparison[i,3] <-object[[i]]$convergence
  modelcomparison[i,4] <-object[[i]]$n_par
  modelcomparison[i,5] <-sum(object[[i]]$n_fact)
  modelcomparison[i,6] <-object[[i]]$n_state
  modelnames <- c(modelnames,paste("[",paste(object[[i]]$n_fact,collapse = ""),"]",sep=""))
}
rownames(modelcomparison) <- modelnames
modelcomparison <-modelcomparison[order(modelcomparison[,"BIC"]),]



modelcomparison2 <- c()
unistates <- unique(modelcomparison[,"n_state"])
for(i in 1:length(unistates)){
  modeli <- modelcomparison[modelcomparison[,"n_state"]==unistates[i],]
  modeli <- modeli[order(modeli[,"n_par"]),]
  local_max <- c(0)
  for(j in 2:nrow(modeli)){
    if(modeli[j,"n_par"] != modeli[j-1,"n_par"]){
      local_max <- c(local_max,sum((modeli[j,"LL"]-modeli[j-1,"LL"])<0))
    }else{
      if((j-2)<1){
        local_max <- c(local_max,0)
      }else{
        local_max <- c(local_max,sum((modeli[j,"LL"]-modeli[j-2,"LL"])<0))
      }
    }
  }
  modeli <- cbind(modeli,local_max)
  modelcomparison2 <- rbind(modelcomparison2,modeli)
}

modelcomparison2 <-modelcomparison2[order(modelcomparison2[,"BIC"]),]

objecModelselection <- (modelcomparison2[,c(1:4,7)])
print(objecModelselection[,-5])
cat("\n")
cat(paste("Note: When re-estimating models that are not convergent, the"),"\n")
cat(paste("   number of maximum iterations should be increased."),"\n")
cat("\n")

}