#' Caclulating initial state and transition probabilities
#'
#' \code{probabilities} calculates initial state and transition probabilities for given covariate scores and time interval for the \code{step3} estimates.
#' 
#' 
#'
#'
#'
#' @param data The dataset (must be a dataframe and contain complete cases only).
#' @param model The model estimated with \code{step1} (must be of class lmfa_step1).
#' @param oblique Indicates whether the factor scores are obtained for the obliquely rotated loadings (TRUE) or unrotated loadings (FALSE) (must be a logical statement).
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#
#'
#'
#' @examples
#' \dontrun{
#' ESM_factorscores <- factorscores_lmfa(data, 
#'                      model, 
#'                      oblique = TRUE,
#'                      rounding = 4
#'                       )
#' }
#' @export


factorscores_lmfa <- function(data, model, oblique = TRUE, rounding = 4){
  
if(missing(data)) stop("argument data is missing, with no default")
if(missing(model)) stop("argument model is missing, with no default")
if(!is.logical(oblique)) stop("argument oblique must be a logical statement")
if(!is.numeric(rounding)) stop("rounding must be a single scalar")
if(length(rounding)>1) stop("rounding must be a single scalar")

if(class(model)!="lmfa_step1") stop("model must be of class lmfa_step1")


  #empty objects with factor scores
  factorscores <- c()
  
  #scale the raw data (i.e., the observations for the indicators)
  #raw_data_scaled <- scale(model$raw_data, center = TRUE, scale = TRUE)
  raw_data_scaled <- as.matrix(model$raw_data)
  
  #empty factor names
  factornames <- c()
  
  for(i in 1:model$n_state){
    #factor names
    factornames <- c(factornames,1:model$n_fact[i])
    
    #calculating regression weights:
    
    #obliquely rotated
    if(oblique==TRUE){
      estimated_cov_matrix <-  (tcrossprod(model$loadings_obli_list[[i]])+ diag(c(model$unique_variances_list[[i]])))
      weights_s = 
        solve(estimated_cov_matrix) %*%              #state-specific estimated covariance matrix
        model$loadings_obli_list[[i]] %*%            #state-specific obliquely rotated loadings
        model$factor_correlations_obli_list[[i]]
    }else{
      #unrotated
      estimated_cov_matrix <-  (tcrossprod(model$loadings_list[[i]])+ diag(c(model$unique_variances_list[[i]])))
      
      weights_s =
        solve(estimated_cov_matrix) %*%              #state-specific estimated covariance matrix
        model$loadings_list[[i]] %*%                 #state-specific unrotated loadings
        diag(model$n_state)                          #diagonal maxtrix (could also be deleted)
    }
    
    #calculate factor scores
    F_s <- 
      raw_data_scaled %*%                            #scaled raw data 
      weights_s                                      #the weights (depening on (un)rotated loadings)
    factorscores <- cbind(factorscores,F_s)
  }
  
  factorscores<- as.data.frame(factorscores)
  colnames(factorscores) <- c(paste("S",rep(1:model$n_state,model$n_fact),"F", factornames,sep=""))
  
  factorscores <- as.matrix(factorscores)
  factorscores <- as.data.frame(factorscores)
  return(cbind(data,factorscores))
}
