#' Compute Multivariate Normal Densities
#'
#' Needed to calculate the posterio-state membership probabilities and the loglikelihood.
#' The function is based on code in dMvn() that is used internally by the mixture of normals model in bayesm.
#'
#' @param x The dataset.
#' @param n_sub Number of observations. #subject or #total observations.
#' @param Lambda_k The state-specific loading matrices.
#' @param Psi_k The state-specific unique variances.
#' @param n_state Number of states. Has to be a scalar.
#' @param J Number of items.
#' @param nu_k The state-specific intercepts.
#'
#' @return Returns the DMV coefficients.
#'
#' @noRd

DMV <- function(x, Lambda_k, Psi_k, n_state, J, n_sub, nu_k) {
  Sigma_k <- rep(list(NA), n_state)
  rooti_k <- rep(list(NA), n_state)
  for (k in 1:n_state) {
    Sigma_k[[k]] <- (tcrossprod(Lambda_k[[k]]) + Psi_k[[k]])
    # rooti_k[[k]] <- backsolve(chol(Sigma_k[[k]]),diag(J))
  }
  saveDMV <- matrix(NA, nrow = n_state, ncol = n_sub) # empty matrix
  for (i in 1:n_sub) { # fill matrix for all observations
    for (k in 1:n_state) { # and for all states
      # quads <- colSums((crossprod(rooti_k[[k]],(t(x[i,])-nu_k[[k]])))^2)
      # saveDMV[k,i] <- exp(-(J/2)*log(2*pi) + sum(log(diag(rooti_k[[k]]))) - .5*quads) #the middle part is the same as -log(det(Sigma_k[[k]]))/2 but it is apparently faster
      saveDMV[k, i] <- mvnpdfC(as.matrix(unlist(x[i, ])), nu_k[[k]], Sigma_k[[k]], Log = FALSE)
    }
  }

  return(saveDMV)
}
