
#' Caclulating initial state and transition probabilities
#'
#' \code{probabilities} calculates initial state and transition probabilities for given covariate scores and time interval for the \code{step3} estimates.
#' 
#' 
#'
#'
#'
#'
#' @param model The model estimated with \code{step3} (must be of class lmfa_step3).
#' @param deltaT The interval for which the transition probabilities should be calculated (must be a single scalar).
#' @param initialCovariateScores The covariate scores for which the probabilities should be calculated (must be a vector with a length equal to the number of covariates that were used for the estimation with \code{step3}). By default the values are set to the sample means of the covariates.
#' @param transitionCovariateScores The covariate scores for which the probabilities should be calculated (must be a vector with a length equal to the number of covariates that were used for the estimation with \code{step3}). By default the values are set to the sample means of the covariates.
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#
#'
#'
#' @examples
#' \dontrun{
#' Probs <- probabilities(model, 
#'          deltaT = 1, 
#'          initialCovariateScores = NULL,
#'          transitionCovariateScores = NULL,
#'          rounding = 4
#'          )
#' }
#' @export


probabilities <- function(model, deltaT = 1, 
initialCovariateScores = NULL, 
transitionCovariateScores = NULL,
rounding = 4){
   
    if(missing(model)) stop("argument model is missing, with no default")
    if(class(model)!="lmfa_step3") stop("model must be of class lmfa_step3")
    if(!is.numeric(deltaT)) stop("deltaT must be a scalar")
    if(deltaT < 0) stop("deltaT must be a positive scalar")
    if(length(deltaT)>1) stop("deltaT must be a single scalar")
    if(!is.numeric(rounding)) stop("rounding must be a single scalar")
    if(length(rounding)>1) stop("rounding must be a single scalar")


    if(!is.null(transitionCovariateScores)){
        if(!is.numeric(transitionCovariateScores)) stop("transitionCovariateScores must be numeric")
        if(length(transitionCovariateScores)!=model$n_transition_covariates) stop("argument transitionCovariateScores must be a vector with a length equal to the number of covariates that were used for the estimation with \code{step3}")
    }

    if(!is.null(initialCovariateScores)){
        if(!is.numeric(initialCovariateScores)) stop("initialCovariateScores must be numeric")
        if(length(initialCovariateScores)!=model$n_initial_covariates) stop("argument initialCovariateScores must be a vector with a length equal to the number of covariates that were used for the estimation with \code{step3}")
    }

    n_state <- model$n_state
    if(is.null(transitionCovariateScores)){
        if(model$n_transition_covariates>0){
            #transitionCovariateScores <- rep(0,model$n_transition_covariates)
            transitionCovariateScores <- model$transition_covariate_means
        }
    }

    if(is.null(initialCovariateScores)){
        if(model$n_initial_covariates>0){
            #initialCovariateScores <- rep(0,model$n_initial_covariates)
            initialCovariateScores <- model$initial_covariate_means
        }
    }


        ini_par <- n_state-1
        tra_par <- n_state*n_state-n_state
        interceptMatrix <- matrix(NA,nrow = ini_par,ncol = 1+length(initialCovariateScores))
        transitionMatrix <-  matrix(NA,nrow = tra_par,ncol = 1+length(transitionCovariateScores))
        
        count <- 0
        interceptMatrix[,1] <- model$estimates[1:ini_par+count,1]
        count <- count+ini_par
        if(!is.null(initialCovariateScores)){
                for(i in 1:length(initialCovariateScores)){
                        interceptMatrix[,i+1] <-model$estimates[1:ini_par+count,1]
                        count <- count+ini_par
                }      
        }
        
        transitionMatrix[,1] <- model$estimates[1:tra_par+count,1]
        count <- count+tra_par
        if(!is.null(transitionCovariateScores)){
                for(i in 1:length(transitionCovariateScores)){
                        transitionMatrix[,i+1] <-model$estimates[1:tra_par+count,1]
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
        
        #return(list("transition probabilities" = TransitionProbabilities,
        #           "initial state probabilities" = InitialStateProbabilities))



    cat(paste("1. Initial state probabilities:","\n"))
    cat("\n")
    if(!is.null(initialCovariateScores)){
        for(i in 1:length(initialCovariateScores)){
           cat(paste(names(model$initial_covariate_means)[i],"score:",round(InitialStateProbabilities[i],rounding),"\n"))
        }
    }else{
        cat(paste("(no covariates defined)","\n"))
        cat("\n")
    }
        names(InitialStateProbabilities) <- paste("S",1:n_state,sep="")
        print(InitialStateProbabilities)

    cat("\n")
    cat(paste("2. Transition probabilities:","\n"))
    cat("\n")
    cat(paste("interval length:",deltaT,"\n"))
    if(!is.null(transitionCovariateScores)){
        for(i in 1:length(transitionCovariateScores)){
           cat(paste(names(model$transition_covariate_means)[i],"score:",round(transitionCovariateScores[i],rounding),"\n"))
        }
    }else{
        cat(paste("(no covariates defined)","\n"))
        cat("\n")
    }

    colnames(TransitionProbabilities) <-paste("S",1:n_state,sep="")
    rownames(TransitionProbabilities) <-paste("S",1:n_state,sep="")
    cat("\n")
    print(round(TransitionProbabilities, rounding))
    cat("\n")

    output <-(list(initial_state_probabilities = InitialStateProbabilities,
                transition_probabilities = TransitionProbabilities,
                initialCovariateScores = initialCovariateScores,
                transitionCovariateScores = transitionCovariateScores,
                deltaT = deltaT,
                initial_covariate_means = model$initial_covariate_means,
                transition_covariate_means = model$transition_covariate_means,
                rounding = rounding
                ))

    class(output) = "lmfa_prob"            
    invisible(output)
}
