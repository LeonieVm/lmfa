#' Conducts steps 1 and 2 from the three-step estimation of CT-LMFA.
#'
#'
#'
#'
#'
#' @param data The dataset (must be a dataframe and contain complete cases only).
#' @param model The model estimated with step1() (must be of class lmfa_step1).
#
#'
#' @return Returns the measurement model parameters, the proportional and
#' modal state assignments, and the classification errors.
#'
#' @examples
#' \dontrun{
#' step2_results <- step2(data, model)
#' }
#' @export

step2 <- function(data, model){
    if(missing(data)) stop("argument data is missing, with no default")
    if(missing(model)) stop("argument model is missing, with no default")
    if(!is.data.frame(data)) stop("data must be a dataframe")
    
    if(class(model)!="lmfa_step1") stop("model must be of class lmfa_step1")
    if(nrow(model$classification_posterior)!=nrow(data)) stop("data must be of same length as data used in step1()")

    
      output <- list(classification_posterior = model$classification_posterior,
                   classification_errors = model$classification_errors,
                   classification_errors_prob = model$classification_errors_prob,
                   R2_entropy = model$R2_entropy,
                   data = cbind(data,model$classification_posterior)
    )
    

    

  class(output) = "lmfa_step2"
  output
}