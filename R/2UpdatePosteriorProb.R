#' Update the state-membership probabilities z_ik
#'
#' Reuses the DMV calculated in the loglikelihoodfunction to save time.
#'
#'
#'
#'
#' @param Lambda_k The state-specific loading matrices.
#' @param Psi_k The state-specific unique variances.
#' @param x The dataset.
#' @param n_state Number of states. Has to be a scalar.
#' @param n_sub Number of observations. #subject or #total observations.
#' @param C_k The sample covariance matrix.
#' @param pi_k The state proportions.
#'
#' @return Returns the expected state-memberships.
#'
#' @examples
#' \dontrun{
#' z_ik <- updExpMem(Lambda_k, Psi_k, n_state, C_k, n_fact, DMV,n_sub,pi_k)
#' }
#' @export

updExpMem <- function(Lambda_k, Psi_k, n_state, C_k, n_fact, DMV,n_sub,pi_k){
  max_z_ik <-rep(0,n_sub)
  sum_pi_k <- matrix(NA,nrow=n_state,ncol=n_sub)
  z_ik <- rep(list(NA),n_state)

  for(i in 1:n_sub){
    #first create part of the denominator
    for(k in 1:n_state){
      sum_pi_k[k,i] <- log(pi_k[[k]]*DMV[k,i]) #I re-use the observation- and state-specific reponse probabilities DMV
    }
    # sum-exp-trick to prevent arithmetic underflow
    max_z_ik[i] <- max(sum_pi_k[,i])
    for(k in 1:n_state){
      sum_pi_k[k,i] <- exp(sum_pi_k[k,i]-max_z_ik[i])
    }
    #create numerator and divide by the sum of the above calculated part
    for(k in 1:n_state){
      z_ik[[k]][i] <- (sum_pi_k[k,i])/sum(sum_pi_k[,i])
    }
  }
  return(z_ik=z_ik)
}
