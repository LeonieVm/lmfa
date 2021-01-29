#' Initialize all Parameters
#'
#' Initial state-membership probabilities are based on mclust and random assignment.
#'
#'
#' @param x The dataset.
#' @param n_sub Number of observations.
#' @param n_state Number of states. Has to be a scalar.
#' @param n_fact Number of factors for each state. Has to be a vector of length n_state.
#' @param J Number of items.
#' @param startval Indicating whether startvalues are based on random assignment or mclust.
#'
#'
#' @return Returns list with all the state-specific initialized parameters.
#'
#' @examples
#' \dontrun{
#' InitialValues <-initializeStep1(
#' x = x, 
#' n_sub = 500, 
#' n_state = 2, 
#' n_fact = c(2,2), 
#' J = 20,
#' startval = "random")
#' z_ik<- InitialValues$z_ik
#' N_k<- InitialValues$N_k
#' pi_k<- InitialValues$pi_k
#' mu_k<- InitialValues$mu_k
#' C_k<- InitialValues$C_k
#' Lambda_k<- InitialValues$Lambda_k
#' Psi_k<- InitialValues$Psi_k
#' }
#' @noRd


initializeStep1 <- function(x,n_sub,n_state,n_fact,J,startval="random"){
  
  
 if(startval=="random"){
    #state-membership probabilities random
    probVector<-runif(n_state,min = 0.1,max=1)
    probVector<-probVector/sum(probVector)
    ini_random <- sample(1:n_state,replace = T,
                         size = nrow(x),
                         prob = probVector)

    z_ik <- rep(list(ini_random),
                n_state) #empty list for posterior state probabilities (as start values)
    for(i in 1:n_state){ #fill list
      z_ik[[i]] <- ifelse(z_ik[[i]]==i,1,0) #give a probability of 1 to the assigned state
    }
  }else{
    #state-membership probabilities with mclust; kxn_sub matrices; 
    #usually mclust is better! (see comment paper Jeroen 2011)
    ini_mclust <- Mclust(x, G =n_state,verbose=FALSE)
    
    ini_mclust <-ini_mclust$classification
    z_ik <- rep(list(ini_mclust),n_state) #empty list for posterior state probabilities (as start values)
    for(i in 1:n_state){ #fill list
      z_ik[[i]] <- ifelse(z_ik[[i]]==i,1,0) #give a probability of 1 to the assigned state
    }
    
  }
  

  #number of observations per state; kx1 vector
  N_k <- lapply(z_ik,sum) #does the same as summing over all observations that are in one state k
  #Reduce("+", N_k) #just a test to see if the number is correct

  #state proportions; kx1 vector
  pi_k <- lapply(z_ik,mean)
  #Reduce("+", pi_k) #just a test to see if the number is correct

  #state-specific intercepts KxJ matrix
  mu_k <-lapply(z_ik, FUN = function(yy) yy * x) #product of state-membership probabilities and data
  mu_k <-lapply(mu_k, colSums) #sum over all observations but not items
  mu_k <-mapply("/",mu_k,N_k,SIMPLIFY = FALSE) #divide by number of observations per state

  #state-specific sample covariance matrix; JxJ matrix
  C_k <- rep(list(NA),n_state) #empty list

  for(sc in 1:n_state){ #fill list
    C_k[[sc]] <- cov.wt(x,z_ik[[sc]],method = 'ML',center = T )$cov #existing function to obtain the weighted cov matrix
  }

  #eigendecomposition; P = matrix of eigenvectors, D = diagonal matrix with eigenvalues
  P_k <- lapply(C_k, function(x) eigen(x)$vectors)  #eigenvectors
  D_k <- lapply(C_k, function(x) eigen(x)$values) #eigenvalues
  
  #save all eigenvalues to use the disregarded onces for calculating the average of the J ??? ff smallest eigenvalues
  D_k_all <-D_k
  
  #get the first ff eigenvalues and vectors per state (i.e., depending on the state-specific number of factors)
  for(k in 1:n_state){ 
    ff <- n_fact[k] #important to define the state-specific number of factors here to allow for different numbers
    P_k[[k]] <- P_k[[k]][,1:ff]#eigenvectors
    D_k[[k]] <- diag(D_k[[k]][1:ff],nrow = ff)#eigenvalues
  }

  #use the eigenvalues and vectors to obtain unique variances and loadings
  aver <-  rep(list(NA),n_state)      #empty list average
  d_k <- rep(list(NA),n_state)        #empty list sqrt of eigenvalues
  Psi_k <- rep(list(NA),n_state)      #empty list unique variances
  for(k in 1:n_state){
    ff <- n_fact[k] #important to define the state-specific number of factors here to allow for different numbers
    aver[[k]] <-  mean(D_k_all[[k]][(ff+1):length(D_k_all[[k]])]) #average of the disregarded eigenvalues
    Psi_k[[k]] <- aver[[k]]*diag(J) #unique variances based on this average
    d_k[[k]] <- sqrt(D_k[[k]]-aver[[k]]*diag(ff)) #sqrt of eigenvalues but not as in normal PCA but as in proportinal PCA
  }
  
  #obtain the loadings
  Lambda_k <- rep(list(NA),n_state) #empty list
  for(k in 1:n_state){ #fill list by calculating the cross-product of the first ff eigenvectors and the sqrt of eigenvalues
    Lambda_k[[k]] <- crossprod(t(P_k[[k]]),d_k[[k]]) #faster version of P_k[[k]]%*%d_k[[k]]
  }
  
  

  return(list(z_ik=z_ik,
              N_k=N_k,
              pi_k=pi_k,
              mu_k=mu_k,
              C_k=C_k,
              Lambda_k=Lambda_k,
              Psi_k=Psi_k))
}
