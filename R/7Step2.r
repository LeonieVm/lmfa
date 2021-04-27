#' Obtaining state assignment and classification errors
#'
#' \code{step2} conducts step 2 of the three-step estimation of CT-LMFA. More specifically, the function extracts the classification information from the \code{step1} output.
#'
#'
#'
#'
#'
#' @param data The dataset (must be a dataframe and contain complete cases only).
#' @param model The model estimated with \code{step1} (must be of class lmfa_step1).
#
#'
#' @return Returns the posterior state probabilities, the modal state assignments, the classification errors, and the state proportions.
#'
#' @examples
#' \dontrun{
#' classification <- step2(data, model)
#' }
#' @export

step2 <- function(data, model){
    if(missing(data)) stop("argument data is missing, with no default")
    if(missing(model)) stop("argument model is missing, with no default")
    if(!is.data.frame(data)) stop("data must be a dataframe")
    
    if(class(model)!="lmfa_step1") stop("model must be of class lmfa_step1")
    if(nrow(model$classification_posteriors)!=nrow(data)) stop("data must be of same length as data used for step1()")
    
    totoal_classification_error <- 1-(sum(diag(model$classification_errors))/nrow(data))
    
      output <- list(classification_posteriors = model$classification_posteriors,
                   classification_errors = model$classification_errors,
                   classification_errors_prob = model$classification_errors_prob,
                   R2_entropy = model$R2_entropy,
                   totoal_classification_error = totoal_classification_error,
                   state_proportions = model$state_proportions,
                   data = cbind(data,model$classification_posteriors)
    )
    

    

  class(output) = "lmfa_step2"
  output
}