#' Conducts step 3 from the three-step estimation of CT-LMFA.
#'
#' 
#'
#'
#'
#'
#' @param data The dataset (must be a dataframe).
#' @param timeintervals The name of the column containing the intervals (must be a single character).
#' @param identifier The name of the column containing the subject identifiers (must be a single character).
#' @param n_state The number of states that should be estimated (must be a single scalar).
#' @param postprobs The posterior state-membership probabilities (must be a dataframe with n_state columns and of same length as the data).
#' @param transitionCovariates The covariate(s) for the transition intensities (must be a (vector of) character(s)).
#' @param initialCovariates The covariate(s) for the initial state probabilities (must be a (vector of) character(s)).
#' @param i.method The type of optimization method that should be used (must be "BFGS" or "CG")
#' @param i.maxit The maximum number of iterations that should be used (must be a single scalar and larger than n_initial_ite).
#' @param i.reltol The tolerance to evaluate convergend that should be used (must be a single scalar).
#' @param i.fnscale An overall scaling to be applied to the value of fn (a function to be minimized) and gr (a function to return the gradient for the "BFGS" and "CG" methods) during optimization (see optim() docomentation for details). In this package it has to be a positive integer.
#' @param n_q The number of start values for the transition intensity parameters that should be used (must be a single scalar).
#' @param n_initial_ite The number of initial iterations for the different start sets that should be used (must be a single scalar).
#' @param rounding The number of decimals to which the results should be rounded (must be a single scalar).
#
#'
#' @return Returns .
#'
#' @examples
#' \dontrun{
#' step3_results <- step3(data,
#'                  identifier,
#'                  n_state,
#'                  postprobs,
#'                  timeintervals = NULL,
#'                  transitionCovariates = NULL,
#'                  initialCovariates = NULL,
#'                  i.method = "BFGS",
#'                  i.maxit = 10000,
#'                  i.reltol = 1e-10,
#'                  i.fnscale = 1,
#'                  n_q = 25,
#'                  n_initial_ite = 15,
#'                  rounding = 4
#'                  )
#' }
#' @export



step3 <- function(data,
                  identifier,
                  n_state,
                  postprobs,
                  timeintervals = NULL,
                  transitionCovariates = NULL,
                  initialCovariates = NULL,
                  i.method = "BFGS",
                  i.maxit = 10000,
                  i.reltol = 1e-10,
                  i.fnscale = "proxi",
                  n_q = 25,
                  n_initial_ite = 10,
                  rounding = 4
                 ){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                                   Step 3
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   #Indicates whether the covariate at t (FALSE) or t-1 (TRUE) should be used (previousCov must be a single logical statement).
   previousCov = FALSE #option that might be added in the future
  
  
  if(missing(data)) stop("argument data is missing, with no default")
  if(missing(identifier)) stop("argument identifier is missing, with no default")
  if(missing(n_state)) stop("argument n_state is missing, with no default")
  if(missing(postprobs)) stop("argument postprobs is missing, with no default")
  #if(nrow(data)!=fitStep1Step2$number_of_timepoints) stop("data must have the same length as the data used in fitStep1Step2")
  if(nrow(postprobs)!=nrow(data)) stop("postprobs must have the same length as data")
  if(ncol(postprobs)!=n_state) stop("the number of columns of postprobs must be of length n_state")
  if(!is.data.frame(data)) stop("data must be a dataframe")
  if(!is.null(transitionCovariates)) if(!is.character(transitionCovariates)) stop("transitionCovariates must be a (vector of) character(s)")
  if(!is.null(initialCovariates)) if(!is.character(initialCovariates)) stop("initialCovariates must be a (vector of) character(s)")
  if(sum(transitionCovariates %in% names(data)) != length(transitionCovariates)) stop("all covariates in transitionCovariates must exist in data")
  if(sum(initialCovariates %in% names(data)) != length(initialCovariates)) stop("all covariates in initialCovariates must exist in data")
  if(i.method != "BFGS") if(i.method != "CG") stop('i.method must be "BFGS" or "CG"')
  if(!is.numeric(i.maxit)) stop("i.maxit must be a single scalar")
  if(length(i.maxit)>1) stop("i.maxit must be a single scalar")
  if(!is.numeric(i.reltol)) stop("i.reltol must be a single scalar")
  if(length(i.reltol)>1) stop("i.reltol must be a single scalar")
  if(length(i.fnscale)>1) stop("i.fnscale must be a single statement")
  if(!is.numeric(i.fnscale) & i.fnscale!="proxi") stop("i.fnscale must be a single scalar")
  if(i.fnscale<1) stop("i.fnscale must be a positive scalar equal to or larger than 1")
  #if(!is.logical(i.center)) stop("i.center must be a single logical statement")
  #if(length(i.center)>1) stop("i.center must be a single logical statement")
  if(!is.numeric(n_q)) stop("n_q must be a single scalar")
  if(length(n_q)>1) stop("n_q must be a single scalar")
  if(!is.numeric(n_initial_ite)) stop("n_initial_ite must be a single scalar")
  if(length(n_initial_ite)>1) stop("n_initial_ite must be a single scalar")
  if(!is.logical(previousCov)) stop("previousCov must be a single logical statement")
  if(length(previousCov)>1) stop("previousCov must be a single logical statement")
  if(!is.numeric(rounding)) stop("rounding must be a single scalar")
  if(length(rounding)>1) stop("rounding must be a single scalar")

  
  if(i.maxit <= n_initial_ite) stop("i.maxit must be larger than n_initial_ite")
  # just a warning for non-specified interval column
  if(is.null(timeintervals)) warning("intervals are assumed to be equidistant because no timeintervals has been specified")
  ptm <- proc.time()
  # i.center Indicates whether covariates are centered at their means during the maximum likelihood estimation (TRUE) or not (FALSE). Centering usually improves stability of the numerical optimisation.
   i.center <- TRUE #option that does not work if not centered: therefore, we always center but report non-centered results. This is not a problem as long as no restictons are made on the intercept (which is not possible in kmfa)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Obtain all necessary elements (from user input or from step 1 and 2).
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Obtain the time_column from timeintervals
  #identifier <- fitStep1Step2$identifier
  n_cases <- length(unlist(unique(data[,identifier])))
  #n_state <- ncol(postprobs)
 
  
  newData <- c()
  if(!is.null(timeintervals)){
    for(i in 1:n_cases){
      datai <- subset(data,get(noquote(identifier))==unlist(unique(data[,identifier]))[i])
      int <- c(0)
      for(ti in 2:nrow(datai)){
        int[ti] <- int[ti-1]+unlist(datai[,timeintervals])[ti]
      }
      datai$time <- int
      newData<-rbind(newData, datai)
    }
  }else{
    for(i in 1:n_cases){
      datai <- subset(data,get(noquote(identifier))==unique(data[,identifier])[i])
      datai$time <- seq(1:nrow(datai))
      newData<-rbind(newData, datai)
    }
  }
  
  # Define the time_column.
  time_column <- "time"


# Obtain the modal state assignments.
  modal_data <- max.col(postprobs)

  Posteriors <-cbind.data.frame(modal_data,postprobs)
  colnames(Posteriors) <- c("Modal", paste("State",1:n_state,sep=""))

  ModalClassificationTable <- matrix(NA,ncol=n_state,nrow=n_state)
  for(i in 1:n_state){
    for(j in 1:n_state){
      ModalClassificationTable[j,i] <- sum((Posteriors[Posteriors$Modal==i,j+1]))
    }
  }

  # Calculate probabilities based on counts.
  W_mod <-ModalClassificationTable/rowSums(ModalClassificationTable)

  # Add the classification from step 1 and 2 to the internally-used dataset.
  newData$State <- c(Posteriors[,"Modal"])

  
  if(sum(W_mod<1e-7)>0){
    W_mod[which(W_mod<1e-7)] <- 1e-7
  }
  
  W_mod <- W_mod/rowSums(W_mod)
  
  #add a proxi for the LL in step 3 (for better fnscale start)
  if(i.fnscale == "proxi"){
  ll_proxi <- 0
  for(i in 1:n_state){
    ll_proxi <- ll_proxi+rowSums(ModalClassificationTable)[i]*
      log(rowSums(ModalClassificationTable)[i]/sum(ModalClassificationTable))
  }
  ll_proxi <-ll_proxi*-2
  i.fnscale <- ll_proxi
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Fixing the right parameters
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  # Define the positions that have to be fixed to equal the entries of the 
  # responseprobabilities in W_mod. Because we work with probabilities that
  # sum to 1, not all n_state*n_state probabilities are fixed but only
  # n_state*n_state-n_state probabilities. Note that the number of covariates 
  # still has to be added.This is based on the additionalCounts defined below.
  n_state_square <- n_state*n_state
  first_element <- n_state_square-n_state+1
  last_element <- n_state_square-n_state+n_state_square-n_state
  fixed_responseprobabilities <- c(first_element:last_element)
  
  #include the transition intensity covariates
  if(length(transitionCovariates)>0){
    #a way to make an equation for the covariate values
    defineCovariates <- as.formula(paste("~", 
                                         paste(transitionCovariates, sep="", 
                                               collapse=" + ")))
  }else{
    #if no covariates were defined
    defineCovariates <-NULL
  }
  additionalCounts <- length(c(fixed_responseprobabilities))*
    length(transitionCovariates)
  
  # The initial state probability covariates are included after covariates are
  # possibly shifted.
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Decide if covariates should be shifted (i.e., consider time-point t) or not
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # If data is shifted: make sure that the initial state covariate is not affected
  # in case the same variable as for predicting the intensities is used.
  
  if(previousCov == FALSE & sum(transitionCovariates%in%initialCovariates)>0){
    for(i in 1:length(initialCovariates)){
      newData[,paste(initialCovariates[i],"x",sep="")] <- 
        newData[,initialCovariates[i]]
      initialCovariatesIntern <- paste(initialCovariates,"x",sep="")
    }
  }else{
    initialCovariatesIntern <-initialCovariates
  }
  
  
  #include the initial state probability covariates
  if(length(initialCovariates)>0){
    #a way to make an equation for the covariate values
    defineInitialCovariates <- as.formula(paste("~", 
                                                paste(initialCovariatesIntern, 
                                                      sep="", collapse=" + ")))
  }else{
    #if no covariates were defined
    defineInitialCovariates <-NULL
  }
  
  # Then shift the data.
  if(length(transitionCovariates)>0){
    #if one want to use the previous covariate, then the data is already correct.
    if(previousCov==FALSE){
      if(length(transitionCovariates)==1){
        newData[,transitionCovariates] <- c(unlist(
          newData[-1,transitionCovariates]),NA)
      }else{
        newData[,transitionCovariates] <- 
          rbind(newData[-1,transitionCovariates],
            rep(NA,length(transitionCovariates)))
      }
    }
  }

  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Multistartprocedure
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  cat("\n")
  cat("1.Initializing...")
  # Obtain a list of initial values (considering the average time-interval)
  # Caclulate average time-interval.
  myid.uni <- unlist(unique(newData[,identifier]))
  myid.uni_length <-length(unlist(myid.uni))
  newData2 <- c()
  for (i in 1:myid.uni_length) {
    temp<-subset(newData, get(noquote(identifier))==myid.uni[i])
    timeDiff <- unlist(temp[-1,time_column]-temp[-nrow(temp),time_column])
    temp$deltaT <- c(NA,as.numeric(timeDiff))
    newData2<-rbind(newData2, temp)
  }
  # Obtain good startvalues based on average
  averageInt <- mean(na.omit(newData2$deltaT))
  unitInt <- 1
  ScalingFactor <- unitInt/averageInt
  
  initialQm <- list()
  initialPm <- list()
  probabilityVector <- runif(n_q,min = 0.5,max=1)
  for(i in 1:n_q){
    Pm <-diag(n_state)
    Pm[Pm==1] <-probabilityVector[i] 
    Pm[Pm==0] <-(1-probabilityVector[i])/(n_state-1 )
    initialPm[[i]] <-Pm
    Qm <- logm(Pm)
    Qm <- Qm*ScalingFactor
    initialQm[[i]] <- Qm
  }
  
  
  # Do n_initial_ite iterations and store the transition intensities
  # We briefly put-off the warnings because there will be one if the
  # model does not converge (and it won't with so few iterations).
  
bestloglik <- list()
q_bestloglik <- list()
identifier<<-identifier #is it problematic to add something to the global environment if I also remove it from inside the function again?
for(i in 1:n_q){
  step3Results <- NULL
  try(
    step3Results <-  suppressWarnings(msm(
    as.formula(paste("State", "~",time_column, sep="")), 
    subject = get(noquote(identifier)),
    data = as.data.frame(newData),
    qmatrix = initialQm[[i]],
    ematrix = W_mod,
    est.initprobs=TRUE,
    hessian = TRUE,
    fixedpars = c(fixed_responseprobabilities)+additionalCounts,
    method  = i.method,
    control = list(maxit = n_initial_ite,reltol = i.reltol,fnscale = i.fnscale),
    center = i.center,
    covariates = defineCovariates,
    initcovariates = defineInitialCovariates))
    , silent = TRUE)
  
  #if(!is.null(step3Results)){
    q_bestloglik[[i]] <-step3Results$Qmatrices$baseline
    bestloglik[[i]] <-step3Results$minus2loglik*-2
  #}
}

if(length(bestloglik)!=0){
  if(length(which(sapply(q_bestloglik, is.null)))>0){
    bestloglik <- bestloglik[-which(sapply(q_bestloglik, is.null))]
    q_bestloglik <- q_bestloglik[-which(sapply(q_bestloglik, is.null))]
  }
}else{
  q_bestloglik <- NULL
}

# Consider the transition intensities from the best start set.
if(is.null(q_bestloglik)){
  stop("numerical overflow; consider changing scale of argument timeintervals and/or using a different value for the argument i.fnscale ")
}else{
  Qm <- q_bestloglik[[which.max(sapply(bestloglik, max))]]
}
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Final Analysis with best startset.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  cat("\n")
  cat(paste("2.Analyzing data with the best start set..."))
  step3Results <-suppressWarnings(msm(
    as.formula(paste("State", "~",time_column, sep="")), 
    subject = get(noquote(identifier)), 
    data = as.data.frame(newData),
    qmatrix = Qm,
    ematrix = W_mod,
    est.initprobs=TRUE,
    hessian = TRUE,
    fixedpars = c(fixed_responseprobabilities)+additionalCounts,
    method  = i.method,
    control=list(maxit = i.maxit,reltol = i.reltol,fnscale= i.fnscale),
    center = i.center,
    covariates = defineCovariates,
    initcovariates = defineInitialCovariates))
                    
#Is the following problematic ? 
CleanEnvir <- function(x) {rm(list=deparse(substitute(x)),envir=.GlobalEnv)}
on.exit(CleanEnvir(identifier))



  requiredTime <- as.numeric((proc.time() - ptm)[3])
  eigenvalues <- eigen(step3Results$paramdata$opt$hessian, only.values = TRUE)$values
  n_eig <- nrow(step3Results$paramdata$opt$hessian)
  for(i in 1: n_eig){
    if(abs(eigenvalues[i])<1e-8) {
      eigenvalues[i]<-0
    }
  }    
  if(any(eigenvalues<=0)){
    warning("Optimisation has probably not converged to the maximum likelihood - Hessian is not positive definite.")
  }
 
  #Has the model reached convergence? (we use 1 as convergence and 0 as non-convergence; thus, we turn the number around)
  if(step3Results$opt$convergence==1){
    convergence <- 0
  }else{
    convergence <- 1
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Extracting the results
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------------------------------------#
  #            Preparation
  #--------------------------------------#

  namesTransitions <- NULL
  for(i in 1:n_state){
        for(j in 1:n_state){
                if(i!=j){
                        namesTransitions <- c(namesTransitions,(paste(i,j,sep = "|")))
                }
        }
  }  

  namesInitial <- NULL
  for(i in 2:n_state){
   namesInitial <- c(namesInitial,i)
  }

  # What is the number of estimated parameters?
  # -For the intensities:
  numberParInt <- additionalCounts+length(c(fixed_responseprobabilities))
  # -For the initial state probabilities
  numberParIni <-n_state-1 + length(initialCovariates)*(n_state-1)
  # -Total number:
  n_par_step3 <- numberParInt+numberParIni
  
  # Create an empty parameter matrix
  parameterEstimates <- matrix(NA,nrow=n_par_step3,ncol=4)
  colnames(parameterEstimates) <- c("coef",
                                    "s.e.",
                                    "z-value",
                                    "p-value")
  # covariate names initial (in case covariate is used twice)
  if(length(initialCovariates)>0){
    iniName <- paste(initialCovariates,"(initial)")
  }else{
    iniName <- NULL
  }
  
  rownames(parameterEstimates) <- 
    #initial state probabilities
  c(paste("initial state parameters",
            namesInitial),
  #cov on initial state probabilities
  paste(rep(iniName,
            each=n_state-1),
        rep(namesInitial,
            length(iniName))),
  
  #intensities
  paste("transition parameters",
        namesTransitions),
  #cov on intensities
  paste(rep(transitionCovariates,
            each=length(fixed_responseprobabilities)),
        rep(namesTransitions,
            length(transitionCovariates))))
  #--------------------------------------#
  #         Parameter estimates
  #--------------------------------------#
  # What are the estimates for the initial state parameters?
  parameterEstimates[1:numberParIni,1] <- 
    step3Results$paramdata$opt$par[(1+numberParInt):n_par_step3]  
  
  # What are the estimates for the transition intensity parameters?
  parameterEstimates[(numberParIni+1):(numberParIni+numberParInt),1] <-
    step3Results$paramdata$opt$par[1:numberParInt]  
  
  #--------------------------------------#
  #            Standard errors
  #--------------------------------------#
  
  # We first calculate the estimated variance-covariance matrix
  # based on the hessian. As the optim() function minimizes (-2)*log(likelihood),
  # we do NOT take the negative of the hessian but multiply it by 0.5
  # (or times 2 after taking the inverse).
  estimatedCovmatrix <- solve(0.5 * 
                                step3Results$paramdata$opt$hessian,
                              silent = TRUE)
  
  # What are the standard errors for the initial state parameters?
  parameterEstimates[1:numberParIni,2] <- 
    suppressWarnings(sqrt(diag(estimatedCovmatrix)))[(1+numberParInt):n_par_step3]
  
  # What are the standard errors for the transition intensity parameters?
  parameterEstimates[(numberParIni+1):(numberParIni+numberParInt),2] <- 
    suppressWarnings(sqrt(diag(estimatedCovmatrix)))[1:numberParInt]
  

#---------------------------------------------------------------#
#Recalculating parameter estimates & SEs (based on delta method)
#---------------------------------------------------------------#



#.........................
#initial state probability
#.........................
if(!is.null(initialCovariates)){
  newSEsIni <- NULL
  if(length(initialCovariates)==1){
    meanvectorInitial <- mean(data[,initialCovariates])
  }else{
    meanvectorInitial <- colMeans(data[,initialCovariates])
  }
  #``````````
  #Parameters
  #``````````
  ParInitialState <- parameterEstimates[1:((n_state-1+
                                              length(initialCovariates)
                                            +(n_state-1)-1)),1]

  ParInitialState <- cbind(ParInitialState,
                           c(rep(0,n_state-1),
                             rep(1:length(initialCovariates),
                                 each=n_state-1)))

  newInterceptInitial <- ParInitialState[ParInitialState[,2]==0,1]

  for(i in 1:length(initialCovariates)){
    newInterceptInitial <- newInterceptInitial-
      meanvectorInitial[i]*ParInitialState[ParInitialState[,2]==i,1]
  }

  #``````````
  #SEs
  #``````````
  #make a matrix with parameters for SE recalculation
  B_matrix <- NULL
  for(i in 1:(n_state-1)){
    B <- ParInitialState[(1:n_state-1),1][i] #take the old intercepts
    for(j in 1:length(initialCovariates)){
      B <-c(B,ParInitialState[ParInitialState[,2]==j,1][i])
    }
    B_matrix <-rbind(B_matrix,B)
  }
  #get the gradients
  grad <- c(1,meanvectorInitial)
  #get correct values of estimated cov matrix

  #note that SEs are at the end; thus, order is different than in parameters
  beginInitial <- (n_state*n_state-n_state+
                     (n_state*n_state-n_state)*
                     length(transitionCovariates))+1
  ParCov <- estimatedCovmatrix[beginInitial:(beginInitial+n_state-1+
                                               (n_state-1)*
                                               length(initialCovariates)-1),
                               beginInitial:(beginInitial+n_state-1+
                                               (n_state-1)*
                                               length(initialCovariates)-1)]
  whichPar <- seq(1, ncol(ParCov), (n_state-1))
  for(i in 1:(n_state-1)){

    vb <- ParCov[whichPar,whichPar]
    B <- B_matrix[i,]
    vG <- grad %*% vb %*% grad
    newSEsIni <- c(newSEsIni,sqrt(vG))
    whichPar <-whichPar+1
  }

parameterEstimates[1:(n_state-1),1:2] <- cbind(newInterceptInitial,newSEsIni)

}
#.........................
#transition intensities
#.........................

if(!is.null(transitionCovariates)){
  newSEs <- NULL
  if(length(transitionCovariates)==1){
    meanvectorTransition <- mean(data[,transitionCovariates])
  }else{
    meanvectorTransition <- colMeans(data[,transitionCovariates])
  }
  #``````````
  #Parameters
  #``````````
  startTran <- ((n_state-1+
                   length(initialCovariates)
                 *(n_state-1)))
  endTran <- startTran+(n_state*n_state-n_state+(
    (n_state*n_state-n_state)*
      length(transitionCovariates)))

  ParTransition <- parameterEstimates[(startTran+1):(endTran),1]

  ParTransition <- cbind(ParTransition,
                           c(rep(0,n_state*n_state-n_state),
                             rep(1:length(transitionCovariates),
                                 each=n_state*n_state-n_state)))

  newTransition <- ParTransition[ParTransition[,2]==0,1]

  for(i in 1:length(transitionCovariates)){
    newTransition <- newTransition-
      meanvectorTransition[i]*ParTransition[ParTransition[,2]==i,1]
  }

  #``````````
  #SEs
  #``````````
  #make a matrix with parameters for SE recalculation
  B_matrix <- NULL
  for(i in 1:(n_state*n_state-n_state)){
    B <- ParTransition[1:(n_state*n_state-n_state),1][i] #take the old intercepts
    for(j in 1:length(transitionCovariates)){
      B <-c(B,ParTransition[ParTransition[,2]==j,1][i])
    }
    B_matrix <-rbind(B_matrix,B)
  }
  #get the gradients
  grad <- c(1,meanvectorTransition)
  #get correct values of estimated cov matrix


  ParCov <- estimatedCovmatrix[1:(n_state*n_state-n_state+
                                    (n_state*n_state-n_state)*
                                    length(transitionCovariates)),
                               1:(n_state*n_state-n_state+
                                    (n_state*n_state-n_state)*
                                    length(transitionCovariates))]
  whichPar <- seq(1, ncol(ParCov), (n_state*n_state-n_state))
  for(i in 1:(n_state*n_state-n_state)){

    vb <- ParCov[whichPar,whichPar]
    B <- B_matrix[i,]
    vG <- grad %*% vb %*% grad
    newSEs <- c(newSEs,sqrt(vG))
    whichPar <-whichPar+1
  }
parameterEstimates[(startTran+1):(startTran+ 
                                (n_state*n_state-n_state)),1:2] <-cbind(newTransition,newSEs)

}






  parameterEstimates[which(is.nan(parameterEstimates))] <- 1000
  
  #--------------------------------------#
  #            z-values
  #--------------------------------------#
  
  # To obtain the z-values, we divide the coefficients by
  # their standard errors
  
  parameterEstimates[,3] <- parameterEstimates[,1]/parameterEstimates[,2]
  
  #--------------------------------------#
  #            p-values
  #--------------------------------------#
  
  parameterEstimates[,4] <- 2*pnorm(-abs(parameterEstimates[,3]))
  
  
  #--------------------------------------#
  #            Wald Test Statistics
  #--------------------------------------#
  
  # We again use the estimated variance-covariance matrix.
  
  parWald <- 1:length(fixed_responseprobabilities)
  whichParameters <- list()
  count <- 0
  
  
  # transition intensities
  whichParameters[[1]] <- 1:length(fixed_responseprobabilities)
  
  # transition intensity covariates
  if(length(transitionCovariates)>0){
    for(i in 1:length(transitionCovariates)){
      whichParameters[[1+i]] <- whichParameters[[1+i-1]]+
        length(fixed_responseprobabilities)
      count <- count+1
    }
  }
  
  # initial state probabilities
  whichParameters[[1+1+ count]] <- (max(
    whichParameters[[1+1+ count-1]])+1):(max(
      whichParameters[[1+1+ count-1]])+n_state-1)
  
  # initial state probability covariates
  if(length(initialCovariates)>0){
    for(i in 1:length(initialCovariates)){
      whichParameters[[1+1+count+i]] <- 
        whichParameters[[1+1+count+i-1]]+n_state-1
    }
  }
  
  allWaldStat <- list()
  allWaldDf <- list()
  # note that there is a 2 because there are always initial state and
  # transition parameters. Then we add possible covariate parameters.
  for(i in 1:(2+length(transitionCovariates)+
              length(initialCovariates))){
    parWald <- whichParameters[[i]]
    ParameterCovarianceMatrix_cov <- 
      estimatedCovmatrix[parWald,parWald]
    parametersWald <- c(step3Results$paramdata$opt$par[parWald])
    #based on the parameters and the values in the variance-covariance matrix
    #we obtain the Wald test statistics.
    WStat <- t(parametersWald)%*%
      solve(ParameterCovarianceMatrix_cov)%*%(parametersWald)
    allWaldStat[[i]] <- WStat
    allWaldDf[[i]] <- length(parWald)
  }
  
  waldMatrix <- matrix(NA, nrow=(2+length(transitionCovariates)+
                                   length(initialCovariates)),ncol=3)
  for(i in 1:nrow(waldMatrix)){
    waldMatrix[i,1] <- allWaldStat[[i]]
    waldMatrix[i,2] <- allWaldDf[[i]]
    #calculate p-values
    waldMatrix[i,3] <- pchisq(
      waldMatrix[i,1],
      df=waldMatrix[i,2],
      lower.tail = F)
  }

  rownames(waldMatrix) <- c("transition intercepts",
                            transitionCovariates,
                            "initial state",
                            iniName)
  colnames(waldMatrix) <- c("Wald(0)",
                            "df",	
                            "p-value")
  
  # reorder just to have the same order as LG
  waldMatrix <- waldMatrix[c("initial state",
                             iniName,
                             "transition intercepts",
                             transitionCovariates),]
  #we do not show the Wald tests for the intercepts because these are based on the wrong SE values/worng covariance matrix entries
    WaldMatrixNoIntercepts <-waldMatrix[-c(which(rownames(waldMatrix)=="initial state"),
                      which(rownames(waldMatrix)=="transition intercepts")),]

  #--------------------------------------#
  #          hessian/covmatrix names
  #--------------------------------------#
  hessianAndCovNames<-
    #initial state probabilities
    c(
      #intensities
      paste("transition intercepts",
            1:length(fixed_responseprobabilities)),
      #cov on intensities
      paste(rep(transitionCovariates,
                each=length(fixed_responseprobabilities)),
            rep(1:length(fixed_responseprobabilities),
                length(transitionCovariates))),
      #initial state          
      paste("initial state",
            1:(n_state-1)),
      #cov on initial state probabilities
      paste(rep(iniName,
                each=n_state-1),
            rep(1:(n_state-1),
                length(iniName))))

    printHessian <- step3Results$paramdata$opt$hessian

    rownames(estimatedCovmatrix) <- hessianAndCovNames
    rownames(printHessian) <- hessianAndCovNames
    colnames(estimatedCovmatrix) <- hessianAndCovNames
    colnames(printHessian) <- hessianAndCovNames
    
  #--------------------------------------#
  #           Vitberi
  #--------------------------------------#
    classification_posterior <- as.matrix(viterbi.msm(step3Results)[,-c(1:2)])
    colnames(classification_posterior) <- c("ModalStep2","Modal",colnames(postprobs))
    classification_posterior <- classification_posterior[,-1] #we do not need the previous assignment anymore
   

    #add updated state proportions
    pi_k <-list()
    for(i in 1:n_state){
      pi_k[[i]] <- as.numeric((table(classification_posterior[,1])/(nrow(classification_posterior)))[i])
    }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                         Return Step 3 Results
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  

  output <-list(
              seconds=round(requiredTime, rounding),
              convergence = convergence,
              LL=round(step3Results$minus2loglik/-2, rounding),
              WaldTests=round(WaldMatrixNoIntercepts,rounding),
              estimates=round(parameterEstimates,rounding), 
              classification_posterior=as.data.frame(classification_posterior),
              pi_k = lapply(pi_k,round,rounding),
              n_transitionCovariates = length(transitionCovariates),
              n_initialCovariates = length(initialCovariates),
              n_state = n_state,
              data = cbind(data, classification_posterior)
              #hessian=printHessian,
              #cov.matrix = estimatedCovmatrix
              )

  class(output) = "lmfa_step3"


  if(output$convergence==1){
    cat("\n")
    cat(paste("Estimation converged.","\n"))
  }else{
    cat("\n")
    cat(paste("Maximum number of iterations reached without convergence.","\n"))
    cat("\n")
  }
  cat("\n")
  cat(paste("LL",round(output$LL,rounding),sep=" = "),"\n")
  cat("\n")

  output
}
