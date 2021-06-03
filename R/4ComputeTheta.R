#' Compute Coefficients Theta
#'
#'
#'
#' @param Beta_k Regression coefficients needed to needed to compute the factor scores from the item scores.
#' @param n_state Number of states. Has to be a scalar.
#' @param C_k The sample covariance matrix.
#' @param Lambda_k The state-specific loading matrices.
#' @param n_fact Number of factors (vector of length n_state).
#'
#' @return Returns the Theta coefficients.
#'
#' @noRd


comTheta <- function(Beta_k, n_state, C_k, Lambda_k, n_fact){
  Theta_k <- rep(list(NA),n_state)
  for(k in 1:n_state){
    Theta_k[[k]] <- diag(n_fact[k])-
      crossprod(t(Beta_k[[k]]),Lambda_k[[k]])+
      tcrossprod(crossprod(t(Beta_k[[k]]),C_k[[k]]),Beta_k[[k]])
  }
  return(Theta_k)
}
