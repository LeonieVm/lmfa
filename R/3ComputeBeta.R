#' Compute Regression coefficients Beta
#'
#' Regression coefficients needed to needed to compute the factor scores from the item scores.
#'
#'
#'
#' @param Lambda_k The state-specific loading matrices.
#' @param Psi_k The state-specific unique variances.
#' @param n_state Number of states. Has to be a scalar.
#' @param n_fact Number of factors for each state. Has to be a vector of length n_state.
#'
#' @return Returns the Beta coefficients.
#'
#' @examples
#' \dontrun{
#' Beta_k <- comBetas(Lambda_k, Psi_k, n_state,n_fact)
#' }
#' @noRd

comBetas <- function(Lambda_k, Psi_k, n_state,n_fact){
  Beta_k <-rep(list(NA),n_state)
  for(k in 1:n_state){
    # calculate inverse once and reuse it because this saves time
    inversePsi_k <- diag(1/diag(Psi_k[[k]])) # same as chol2inv(chol(Psi_k[[k]])) but faster for large matrices
    # Woodbury Identity for a more efficient computation of the inverse
    # note that transposes may look different than in the formula.
    # this is because I use a crossprod function that is slightly faster
    # than the %*% function but requires some extra transposes sometimes
    Beta_k[[k]] <- crossprod(Lambda_k[[k]],(
      inversePsi_k-
        tcrossprod(tcrossprod(crossprod(t(crossprod(t(inversePsi_k),Lambda_k[[k]])),
                                        chol2inv(chol(diag(n_fact[k])+
                                                        tcrossprod(crossprod(Lambda_k[[k]],inversePsi_k),t(Lambda_k[[k]]))))),
                              Lambda_k[[k]]),
                   t(inversePsi_k))
    ))
  }
  return(Beta_k)
}


#same code but much slower
# t(Lambda_k[[k]])%*%(
#   inversePsi_k-
#     inversePsi_k%*%Lambda_k[[k]]%*%
#     chol2inv(chol(diag(n_fact[k])+t(Lambda_k[[k]])%*%inversePsi_k%*%Lambda_k[[k]]))%*%
#     t(Lambda_k[[k]])%*%
#     inversePsi_k
# )
