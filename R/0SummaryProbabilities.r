#' Summary for \code{probabilities}
#'
#'
#'
#'
#'
#'
#' @param object An object storing output from \code{probabilities}.
#' @param ... Further arguments for the default S3 summary method.
#' @examples
#' \dontrun{
#' summary(transitionmodel)
#' }
#' @export
summary.lmfa_prob <- function(object,...){

  
    initial_state_probabilities <- object$initial_state_probabilities
    transition_probabilities <-object$transition_probabilities
    initialCovariateScores <- object$initialCovariateScores
    transitionCovariateScores <- object$transitionCovariateScores
    deltaT <- object$deltaT
    initial_covariate_means <- object$initial_covariate_means
    transition_covariate_means <- object$transition_covariate_means
    rounding <- object$rounding

    cat(paste("1. Initial state probabilities:","\n"))
    cat("\n")
    if(!is.null(initialCovariateScores)){
        for(i in 1:length(initialCovariateScores)){
           cat(paste(names(initial_covariate_means)[i],"score:",round(initial_state_probabilities[i],rounding),"\n"))
        }
    }else{
        cat(paste("(no covariates defined)","\n"))
        cat("\n")
    }
        #names(initial_state_probabilities) <- paste("S",1:n_state,sep="")
        print(initial_state_probabilities)

    cat("\n")
    cat(paste("2. Transition probabilities:","\n"))
    cat("\n")
    cat(paste("interval length:",deltaT,"\n"))
    if(!is.null(transitionCovariateScores)){
        for(i in 1:length(transitionCovariateScores)){
           cat(paste(names(transition_covariate_means)[i],"score:",round(transitionCovariateScores[i],rounding),"\n"))
        }
    }else{
        cat(paste("(no covariates defined)","\n"))
        cat("\n")
    }

    #colnames(transition_probabilities) <-paste("S",1:n_state,sep="")
    #rownames(transition_probabilities) <-paste("S",1:n_state,sep="")
    cat("\n")
    print(round(transition_probabilities, rounding))
    cat("\n")


}