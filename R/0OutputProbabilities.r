
#' Calculates probabilities for given covariate scores and time interval.
#'
#' 
#'
#'
#'
#'
#' @param x The model estimated with step3() (must be of class lmfa_step3).
#' @param deltaT The interval for which the transition probabilities should be calculated (must be a single scalar).
#' @param transitionCovariateScores The covariate scores for which the probabilities should be calculated (must be a vector with a length equal to the number of covariates that were used for the estimation in step3()). By default scores are set to zero.
#' @param initialCovariateScores The covariate scores for which the probabilities should be calculated (must be a vector with a length equal to the number of covariates that were used for the estimation in step3()). By default scores are set to zero.
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#
#'
#' @return Returns .
#'
#' @examples
#' \dontrun{
#' Probs <- probabilities(x, 
#'          deltaT = 1, 
#'          transitionCovariateScores = NULL,
#'          initialCovariateScores = NULL,
#'          rounding = 4
#'          )
#' }
#' @export


probabilities <- function(x, deltaT = 1, transitionCovariateScores = NULL,initialCovariateScores = NULL, rounding = 4){
   
    if(missing(x)) stop("argument x is missing, with no default")
    if(class(x)!="lmfa_step3") stop("x must be of class lmfa_step3")
    if(!is.numeric(deltaT)) stop("deltaT must be a single scalar")
    if(length(deltaT)>1) stop("deltaT must be a single scalar")
    if(!is.numeric(rounding)) stop("rounding must be a single scalar")
    if(length(rounding)>1) stop("rounding must be a single scalar")


    if(!is.null(transitionCovariateScores)){
        if(!is.numeric(transitionCovariateScores)) stop("transitionCovariateScores must be numeric")
        if(length(transitionCovariateScores)!=x$n_transitionCovariates) stop("argument transitionCovariateScores must be a vector with a length equal to the number of covariates that were used for the estimation in step3()")
    }

    if(!is.null(initialCovariateScores)){
        if(!is.numeric(initialCovariateScores)) stop("initialCovariateScores must be numeric")
        if(length(initialCovariateScores)!=x$n_initialCovariates) stop("argument initialCovariateScores must be a vector with a length equal to the number of covariates that were used for the estimation in step3()")
    }

    n_state <- x$n_state
    if(is.null(transitionCovariateScores)){
        if(x$n_transitionCovariates>0){
            transitionCovariateScores <- rep(0,x$n_transitionCovariates)
        }
    }

    if(is.null(initialCovariateScores)){
        if(x$n_initialCovariates>0){
            initialCovariateScores <- rep(0,x$n_initialCovariates)
        }
    }


        ini_par <- n_state-1
        tra_par <- n_state*n_state-n_state
        interceptMatrix <- matrix(NA,nrow = ini_par,ncol = 1+length(initialCovariates))
        transitionMatrix <-  matrix(NA,nrow = tra_par,ncol = 1+length(transitionCovariates))
        
        count <- 0
        interceptMatrix[,1] <- x$estimates[1:ini_par+count,1]
        count <- count+ini_par
        if(!is.null(initialCovariates)){
                for(i in 1:length(initialCovariates)){
                        interceptMatrix[,i+1] <-x$estimates[1:ini_par+count,1]
                        count <- count+ini_par
                }      
        }
        
        transitionMatrix[,1] <- x$estimates[1:tra_par+count,1]
        count <- count+tra_par
        if(!is.null(transitionCovariates)){
                for(i in 1:length(transitionCovariates)){
                        transitionMatrix[,i+1] <-x$estimates[1:tra_par+count,1]
                        count <- count+tra_par
                }      
        }
        
        
        covVectorTr <- c(1,transitionCovariateScores)
        covVectorIni <- c(1,initialCovariateScores)
        
        #initial state probabilities
        initialProbs <- interceptMatrix %*% covVectorIni
        initialProbsExp <- initialProbs
        for(i in 1:ini_par){
                initialProbsExp[i] <-  exp(initialProbs[i])/sum(exp(initialProbs),1)
        }
        initialProbsExp <- c(1-sum(initialProbsExp),initialProbsExp)
        InitialStateProbabilities <- round(initialProbsExp,rounding)
        
        #transition probabilities
        qmatrix <- transitionMatrix %*% covVectorTr
        
        LogIntensities <- diag(n_state)
        LogIntensities[LogIntensities==0] <- qmatrix
        LogIntensities <- t(LogIntensities)
        Intensities <- exp(LogIntensities)
        
        for(i in 1: n_state){
                Intensities[i,i] <- -(sum(Intensities[i,-i]))
        }
        
        TransitionProbabilities <-expm(Intensities *deltaT)
        TransitionProbabilities <- round(TransitionProbabilities,rounding)
        
        return(list("transition probabilities" = TransitionProbabilities,
                    "initial state probabilities" = InitialStateProbabilities))
}
