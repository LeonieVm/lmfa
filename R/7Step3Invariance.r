#' Evaluating for which individuals invariance holds across all measurement occasions
#'
#' \code{invariance} evaluates for which individuals invariance holds across all measurement occasions based on the \code{step3} output.
#' 
#' 
#'
#'
#'
#'
#' @param object  Output from the main function step3() (must be of class lmfa_step3).
#' @param identifier The name of the column containing the subject identifiers (must be a single character).
#
#'
#'
#' @examples
#' \dontrun{
#' invariance(object, identifier)
#' }
#' @export


invariance <- function(object, identifier){
        if(missing(object)) stop("argument data is missing, with no default")
        if(missing(identifier)) stop("argument identifier is missing, with no default")
        if(class(object)!="lmfa_step3") stop("model must be of class lmfa_step3")


        dataInvariance <- object$data
        invariance_time <- c()
        statemembership <- c()
        subjects <- unique(dataInvariance[,identifier])
        for(i in 1:length(subjects)){
        subject_number <- unique(dataInvariance[,identifier])[i]
        
        if(length(unique(dataInvariance[dataInvariance[,identifier]==subject_number,"Modal"]))==1){
            invariance_time <- c(invariance_time,subject_number)
            statemembership <- c(statemembership,unique(dataInvariance[dataInvariance[,identifier]==subject_number,"Modal"]))
        }
        }
        if(!is.null(statemembership)){
            statemembership <- as.matrix(statemembership,ncol=1)
        rownames(statemembership) <- invariance_time


        invariance_subjects <- list()
        for(i in 1:length(sort(unique(statemembership)))){
        invariance_subjects[[i]]<- which(statemembership==i)
        }

        names(invariance_subjects) <- paste("S",sort(unique(statemembership)),sep="")

        print(invariance_subjects)
        }else{
            cat(paste("There are no individuals for whom invariance over time holds.","\n"))
            cat("\n")
        }

}