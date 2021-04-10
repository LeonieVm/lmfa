#' Summary statistics for step1()
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from the main function Step1Step2()
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#' @param ... Further arguments for the default S3 summary method
#'
#' @examples
#' \dontrun{
#' summary(results1)
#' }
#' @export
summary.lmfa_step1 <- function(object, rounding = 2,...){
    if(!is.numeric(rounding)) stop("rounding must be a single scalar")
    if(length(rounding)>1) stop("rounding must be a single scalar")

  if(object$convergence==1){
    cat("\n")
    cat(paste("Estimation converged after",round(object$seconds,2),"seconds and",object$n_it,"iterations.","\n"))
  }else{
    cat("\n")
    cat(paste("Maximum number of iterations reached without convergence.","\n"))
    cat("\n")
  }
  cat("\n")
  cat(paste("LL",round(object$LL,rounding),sep=" = "),"\n")
  cat("\n")
  cat(paste("Number of states",object$n_state,sep=": "),"\n")
  cat("\n")
  cat(paste("Number of factors: [", paste(object$n_fact, collapse = " "), "]", sep = ""),"\n")
  cat("\n")
  cat("\n")
  
  cat(paste("Obliquely rotated within-state standardized loadings:"),"\n")
  cat("\n")
  J <- nrow(object$intercepts)
  n_fact <- object$n_fact
  n_state <- object$n_state
  
  steps <- c()
  count <- 0
  for(i in 1:n_state){
    count <- count+n_fact[i]
    steps <- c(steps,count)
  }
  starti <- steps-(n_fact-1)

  loadings <- as.data.frame(object$loadings_w_obli)
  for(i in 1:n_state){
    loadings_new <- cbind.data.frame(loadings[,starti[i]:steps[i], drop = FALSE],"  " = c(rep("  ",J)))
    if(i==1){
      loadings_new_1 <- loadings_new
    }else{
      loadings_new_1 <- cbind.data.frame(loadings_new_1,loadings_new)
    }
  }
  
  print(round_df(loadings_new_1, digits=rounding))
  cat("\n")
  cat(paste("Obliquely rotated between-state standardized loadings:"),"\n")
  cat("\n")

  loadings <- as.data.frame(object$loadings_b_obli)
  for(i in 1:n_state){
    loadings_new <- cbind.data.frame(loadings[,starti[i]:steps[i], drop = FALSE],"  " = c(rep("  ",J)))
    if(i==1){
      loadings_new_1 <- loadings_new
    }else{
      loadings_new_1 <- cbind.data.frame(loadings_new_1,loadings_new)
    }
  }
  
  print(round_df(loadings_new_1, digits=rounding))
  cat("\n")

  cat(paste("Factor correlations after oblique rotation:"),"\n")
  cat("\n")
  
  factorcors <- object$factor_correlations_obli_list

  for(i in 1:n_state){
    if(n_fact[i]>1){
      colnames(factorcors[[i]]) <- paste("F", 1:n_fact[i], sep = "")
      rownames(factorcors[[i]]) <- paste("F", 1:n_fact[i], sep = "")
    }else{
      names(factorcors[[i]]) <- paste("F", 1:n_fact[i], sep = "")
    }
  }


  for(i in 1:n_state){
    cat(paste("S",i,sep=""),"\n")
    print(lapply(factorcors,round,rounding)[[i]])
    cat("\n")
  }

  cat("\n")
  cat(paste("Intercepts:"),"\n")
  cat("\n")
  print(round(object$intercepts,rounding))
  cat("\n")

 
  cat(paste("Unique variances:"),"\n")
  cat("\n")
  print(round(object$unique_variances,rounding))
  cat("\n")
}

