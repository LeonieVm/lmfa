#' Update Loadings and Unique Variances
#'
#'
#'
#'
#'
#' @param n_state Number of states. Has to be a scalar.
#' @param C_k The sample covariance matrix.
#' @param n_fact Number of factors (vector of length n_state).
#' @param Beta_k Regression coefficients needed to needed to compute the factor scores from the item scores.
#' @param Theta_k Expectation of squared factor scores given the data.
#' @param J Number of items.
#'
#' @return Returns the loading and unique variance parameter.
#'
#' @examples
#' \dontrun{
#' LambPsi<-updLambPsi(n_state, C_k, n_fact, Beta_k, Theta_k)
#' Lambda_k <- LambPsi$Lambda_k
#' Psi_k <- LambPsi$Psi_k
#' }
#' @noRd


updLambPsi <- function(n_state, C_k, n_fact, Beta_k, Theta_k,residualVariance,J){
  Lambda_kN <- rep(list(NA),n_state)
  Psi_kN <- rep(list(NA),n_state)
  act_constraints <-0
  #Update Lambda and Psi in one loop
  for(k in 1:n_state){
    Lambda_kN[[k]]<- tcrossprod(tcrossprod(C_k[[k]],(Beta_k[[k]])),chol2inv(chol(Theta_k[[k]])))
    diagonalElementsPsi <- diag(C_k[[k]]-tcrossprod(crossprod(t(Lambda_kN[[k]]),Beta_k[[k]]),C_k[[k]]))
    #make sure there are no negative variances
    
    for(rv in 1:J){
      if(diagonalElementsPsi[rv]<residualVariance[rv]){
        act_constraints <- act_constraints+1
        diagonalElementsPsi[rv] <- residualVariance[rv]
      }
    }
    
    
    Psi_kN[[k]]<- diag(diagonalElementsPsi)
  }
  #replace the old lambda and psi scores
  Lambda_k <- Lambda_kN
  Psi_k <- Psi_kN

  return(list(Lambda_k=Lambda_k,Psi_k=Psi_k,act_constraints=act_constraints))
}


# updLambPsi <- function(n_state, C_k, n_fact, Beta_k, Theta_k,ActivatedVarianceConstraints){
#   Lambda_kN <- rep(list(NA),n_state)
#   Psi_kN <- rep(list(NA),n_state)
#   act_constraints <-0
#   #Update Lambda and Psi in one loop
#   for(k in 1:n_state){
#     Lambda_kN[[k]]<- tcrossprod(tcrossprod(C_k[[k]],(Beta_k[[k]])),chol2inv(chol(Theta_k[[k]])))
#     diagonalElementsPsi <- diag(C_k[[k]]-tcrossprod(crossprod(t(Lambda_kN[[k]]),Beta_k[[k]]),C_k[[k]]))
#     #make sure there are no negative variances
#     act_constraints <- act_constraints+sum(diagonalElementsPsi<0.0003)
#     diagonalElementsPsi <- ifelse(diagonalElementsPsi<0.0003,0.0003,diagonalElementsPsi)
#     Psi_kN[[k]]<- diag(diagonalElementsPsi)
#   }
#   #replace the old lambda and psi scores
#   Lambda_k <- Lambda_kN
#   Psi_k <- Psi_kN
#   
#   return(list(Lambda_k=Lambda_k,Psi_k=Psi_k,act_constraints=act_constraints,Constraints=ActivatedVarianceConstraints))
# }

# updLambPsi <- function(n_state, C_k, n_fact, Beta_k, Theta_k,ActivatedVarianceConstraints){
#   Lambda_kN <- rep(list(NA),n_state)
#   Psi_kN <- rep(list(NA),n_state)
#   act_constraints <-0
#   #Update Lambda and Psi in one loop
#   for(k in 1:n_state){
#     Lambda_kN[[k]]<- tcrossprod(tcrossprod(C_k[[k]],(Beta_k[[k]])),chol2inv(chol(Theta_k[[k]])))
#     diagonalElementsPsi <- diag(C_k[[k]]-tcrossprod(crossprod(t(Lambda_kN[[k]]),Beta_k[[k]]),C_k[[k]]))
#     
#     if(iteration==1){
#       #make sure there are no negative variances
#       act_constraints <- act_constraints+sum(diagonalElementsPsi<0.0003)
#       diagonalElementsPsi <- ifelse(diagonalElementsPsi<0.0003,0.0003,diagonalElementsPsi)
#       Psi_kN[[k]]<- diag(diagonalElementsPsi)
#       ActivatedVarianceConstraints[[k]] <- which(diagonalElementsPsi<0.0003)
#     }else{
#       if(estimation[iteration-1,"act_constraints"]==0){
#         #make sure there are no negative variances
#         act_constraints <- act_constraints+sum(diagonalElementsPsi<0.0003)
#         ActivatedVarianceConstraints[[k]] <- which(diagonalElementsPsi<0.0003)
#         diagonalElementsPsi <- ifelse(diagonalElementsPsi<0.0003,0.0003,diagonalElementsPsi)
#         Psi_kN[[k]]<- diag(diagonalElementsPsi)
#         
#       }else{
#         diagonalElementsPsi[ActivatedVarianceConstraints[[k]]] <- 0.0003
#         Psi_kN[[k]]<- diag(diagonalElementsPsi)
#         act_constraints <-999
#       }
#     }
#   }
#   #replace the old lambda and psi scores
#   Lambda_k <- Lambda_kN
#   Psi_k <- Psi_kN
#   
#   return(list(Lambda_k=Lambda_k,Psi_k=Psi_k,act_constraints=act_constraints,Constraints=ActivatedVarianceConstraints))
# }




