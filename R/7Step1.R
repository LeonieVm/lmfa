#' Estimating state-specific measurement models
#'
#' \code{step1} conducts step 1 of the three-step estimation of LMFA and thus the estimation of the measurement models. It is possible to estimate the parameters for a single model or for a range of models and to conduct model selection.
#'
#'
#'
#'
#'
#'
#' @param data The dataset (must be a dataframe and contain complete cases only).
#' @param indicators The variable names of the indicators (must be a vector of characters).
#' @param n_state The number of states that should be estimated (must be a single scalar).
#' @param n_fact The number of factors per state that should be estimated (must be a numeric vector of length n_state).
#' @param modelselection Indicates whether model selection should be performed or not. If TRUE, the arguments n_state_range and n_fact_range are required (must be a logical statement).
#' @param n_state_range A vector indicating the number of states that should be considered in the model selection.
#' @param n_fact_range A vector indicating the number of factors that should be considered in the model selection.
#' @param n_starts The number of random starts that should be used (must be a single scalar).
#' @param n_initial_ite The number of initial iterations for the best starts (must be a single scalar).
#' @param n_m_step The number of M-step iterations that should be used when parameters still change more than defined by the m_step_tolerance (must be a single scalar).
#' @param em_tolerance The convergence criterion for parameters and loglikelihood (must be a single scalar and smaller than m_step_tolerance).
#' @param m_step_tolerance The criterion for stopping the n_m_step M-step iterations (must be a single scalar).
#' @param max_iterations The maximum number of iterations (must be a single scalar and larger than n_initial_ite).
#' @param n_mclust The number of mclust starts (must be a single scalar and at least equal to 2).
#'
#' @return Returns the state-specific measurement model parameters and model fit information (for one or multiple estimated model).
#'
#' @examples
#' \dontrun{
#' step1_results <- step1(data,
#'                  indicators,
#'                  n_state = NULL,
#'                  n_fact = NULL, 
#'                  modelselection = FALSE, 
#'                  n_state_range = NULL, 
#'                  n_fact_range = NULL,
#'                  n_starts = 25,
#'                  n_initial_ite = 10,
#'                  n_m_step = 10,
#'                  em_tolerance = 1e-8, 
#'                  m_step_tolerance = 1e-3, 
#'                  max_iterations = 1000,
#'                  n_mclust = 5
#'                  )
#' }
#' @export

step1 <- function(data,
                  indicators,
                  n_state = NULL,
                  n_fact = NULL, 
                  modelselection = FALSE, 
                  n_state_range = NULL, 
                  n_fact_range = NULL,
                  n_starts = 25,
                  n_initial_ite = 15,
                  n_m_step = 10,
                  em_tolerance = 1e-8, 
                  m_step_tolerance = 1e-3, 
                  max_iterations = 1000,
                  n_mclust = 5 
                  ){
  rounding = 12
  if(missing(data)) stop("argument data is missing, with no default")
  if(missing(indicators)) stop("argument indicators is missing, with no default")
  if(length(modelselection)>1) stop("modelselection must be a single logical statement")
  if(!is.logical(modelselection)) stop("argument modelselection must be a logical statement")
  #no modelselection
  if(modelselection == FALSE){
  if(is.null(n_state)) stop("argument n_state is missing, with no default")
  if(is.null(n_fact)) stop("argument n_fact is missing, with no default")
  if(!is.numeric(n_state)) stop("n_state must be a single scalar")
  if(length(n_state)>1) stop("n_state must be a single scalar")
  if(length(n_starts)>1) stop("n_starts must be a single scalar")
  if(!is.numeric(n_fact)) stop("n_fact must be a numeric vector")
  if(length(n_fact)!=n_state) stop("n_fact must be of length n_state")
  if(!is.null(n_state_range)){
    cat("\n")
    cat("argument n_state_range will be overwritten by n_state because modelselection == FALSE") 
    cat("\n")
  } 
  if(!is.null(n_fact_range)){
    cat("\n")
    cat("argument n_fact_range will be overwritten by n_fact because modelselection == FALSE") 
    cat("\n")
  }
  #modelselection
  }else{
  if(is.null(n_state_range)) stop("argument n_state_range is missing, with no default")
  if(is.null(n_fact_range)) stop("argument n_fact_range is missing, with no default")
  if(!is.null(n_state)){
    cat("\n")
  cat("argument n_state will be overwritten by n_state_range because modelselection == TRUE") 
  cat("\n")
  } 
  if(!is.null(n_fact)){
    cat("\n")
    cat("argument n_fact will be overwritten by n_fact_range because modelselection == TRUE") 
    cat("\n")
  } 

  }

  if(!is.data.frame(data)) stop("data must be a dataframe")
  if(!is.character(indicators)) stop("indicators must be a vector of characters")
  if(!is.numeric(n_initial_ite)) stop("n_initial_ite must be a single scalar")
  if(length(n_initial_ite)>1) stop("n_initial_ite must be a single scalar")
  if(!is.numeric(n_m_step)) stop("n_m_step must be a single scalar")
  if(length(n_m_step)>1) stop("n_m_step must be a single scalar")
  if(!is.numeric(em_tolerance)) stop("em_tolerance must be a single scalar")
  if(length(em_tolerance)>1) stop("em_tolerance must be a single scalar")
  if(em_tolerance<0) stop("em_tolerance must be a positive scalar")
  if(!is.numeric(m_step_tolerance)) stop("m_step_tolerance must be a single scalar")
  if(length(m_step_tolerance)>1) stop("m_step_tolerance must be a single scalar")
  if(m_step_tolerance<0) stop("m_step_tolerance must be a positive scalar")
  if(length(max_iterations)>1) stop("max_iterations must be a single scalar")
  if(!is.numeric(max_iterations)) stop("max_iterations must be a single scalar")
  if(max_iterations<0) stop("max_iterations must be a positive scalar")
  if(sum(complete.cases(data[,indicators])==FALSE)>0) stop("data must contain complete cases with regard to the indicators only")
  if(max_iterations <= n_initial_ite) stop("max_iterations must be larger than n_initial_ite")
  if(n_initial_ite<0) stop("n_initial_ite must be a positive scalar")
  if(em_tolerance >= m_step_tolerance) stop("em_tolerance must be smaller than m_step_tolerance")
  if(!is.numeric(rounding)) stop("rounding must be a single scalar")
  if(length(rounding)>1) stop("rounding must be a single scalar")
  if(!is.numeric(n_mclust)) stop("n_mclust must be a single scalar")
  if(length(n_mclust)>1) stop("n_mclust must be a single scalar")
  if(n_mclust<2) stop("n_mclust must be larger than 1")
  
  if(!is.null(n_state_range)){
    if(sum(round(n_state_range)!=n_state_range)>0) stop("n_state_range must be a sequence of integers")
  }

  if(!is.null(n_fact_range)){
    if(sum(round(n_fact_range)!=n_fact_range)>0) stop("n_fact_range must be a sequence of integers")
  }

  if(!is.null(n_state)){
    if(round(n_state)!=n_state) stop("n_state must be an integer")
  }

  if(!is.null(n_fact)){
    if(sum(round(n_fact)!=n_fact)>0) stop("n_fact must be a vector of integers")
  }

  if(round(n_starts)!=n_starts) stop("n_starts must be an integer")
  if(round(n_initial_ite)!=n_initial_ite) stop("n_initial_ite must be an integer")
  if(round(n_m_step)!=n_m_step) stop("n_m_step must be an integer")
  if(round(max_iterations)!=max_iterations) stop("max_iterations must be an integer")
  if(round(n_mclust)!=n_mclust) stop("n_mclust must be an integer")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Input: b) defined in the package.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  # Perhaps set a seed in order to replicate results (e.g., same starting values).
  

  # Obtain the columns with the variables
  x <- data[,indicators]
  x <- as.data.frame(x)
  if(sum(is.na(x)>0)) stop("data contains missing values on indicator variables that must be removed")
  
 
  #*******************************************************************************#
  # NOTE: Usualy, in mixture factor analysis, this would be the number
  # of subjects but here this number represents the independently treated total
  # number of observations.
  #*******************************************************************************#
  # Number of observations.
  n_sub <- nrow(x)

  # Number of items.
  J <- ncol(x)
  if(sum(apply(x, 2, is.numeric))!=J) stop("indicator variables must be numeric")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                        Setup Modelselection
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if(modelselection == TRUE){
        n_fact_range <- unique(n_fact_range)
        allModels <- list()
        all_estimated_models <-list()
        n_models <- 0
        for(i in 1:length(n_state_range)){
                MMcombi <- combinations_K_F_k(length(n_fact_range), n_state_range[i], n_fact_range)
                allModels[[i]] <- MMcombi
                n_models <- n_models+nrow(MMcombi)
        }
        ModelMatrix <- matrix(NA,nrow = n_models,ncol = max(n_state_range))
        skip <- 0
        for(i in 1:length(n_state_range)){
                n_state_specific_model <- nrow(allModels[[i]])
                ModelMatrix[(1:n_state_specific_model+skip),1:ncol(allModels[[i]])] <- allModels[[i]]
                skip <- skip+n_state_specific_model
        }  
}else{
        n_models <- 1
}
  
  cat("\n")
  cat("-------------------------------------------------------------")
  cat("\n")
  # for all models in the model selection (if modelselection == FALSE, only one model is estimated)
  for(comparingmodels in 1:n_models){
  ptm <- proc.time()
  cat("\n")
  cat(paste("Model",comparingmodels,"out of", n_models,sep=" "),"\n")
  if(modelselection==TRUE){
    currentmodel <- ModelMatrix[comparingmodels,]
    currentmodel <- currentmodel[!is.na(currentmodel)]
    n_state <- length(currentmodel)
    n_fact <- currentmodel
  }

  cat("\n")
  cat(paste("Number of states",n_state,sep=": "),"\n")
  cat("\n")
  cat(paste("Number of factors: [", paste(n_fact, collapse = " "), "]", sep = ""),"\n")
  cat("\n")

  # List of multistart procedure results.
  MultistartResults1 <- rep(list(list(NA)),n_starts*10)
  MultistartResults2 <- rep(list(list(NA)),n_starts)
  MclustResults <- rep(list(list(NA)),3)

  # Prepare storing iterations.
  estimation <- matrix(c(NA),ncol=3,nrow=max_iterations)
  colnames(estimation) <- c("interation","loglik","activ. constraints")

  # Define minimum amount of residual variances.
  residualVariance <- rep(NA,J)
  for(j in 1:J){
    residualVariance[j] <- (var(x[,j]) * (n_sub - 1) / n_sub)*1.0e-6
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                        Step 1: Measurement Models
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  
  cat("\n")
  cat("1.Initializing...")
  # Start parallelization.
  getDoParWorkers()
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)

  iteration <- 1
  SumParameterChange <- 100

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Initialize based on mclust 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ini_mclust_list <- list()

  for(mcluststarts in 1:n_mclust){
    ini_mclust <- Mclust(x, G = n_state, verbose=FALSE)
    ini_mclust <- ini_mclust$classification
    ini_mclust_list[[mcluststarts]] <- ini_mclust


    # Self-created function (see '1InitializeEM.R').
    InitialValues <- initializeStep1(x,n_sub,n_state,n_fact,
                                    J,startval="mclust",
                                    RandVec=RandVec,
                                    ini_mclust = ini_mclust,
                                    ini_mclust_specific = ini_mclust_specific)#mclust;


    # Extract all parameters that are going to be updated in the EM algorithm.
    z_ik<- InitialValues$z_ik         #expected state-membership-probabilities
    N_k<- InitialValues$N_k           #sample size per state
    pi_k<- InitialValues$pi_k         #state proportions
    nu_k<- InitialValues$nu_k         #state-specific intercepts
    C_k<- InitialValues$C_k           #sample covariance matrix
    Lambda_k<- InitialValues$Lambda_k #state-specific loading matrices
    Psi_k<- InitialValues$Psi_k       #state-specific unique variances

    #-------------------------------------------------------------------------------#
    # Compute: Regression Weights (Beta)
    #-------------------------------------------------------------------------------#

    # Self-created function (see '3ComputeBeta.R').
    Beta_k <- comBetas(Lambda_k, Psi_k, n_state,n_fact)

    #-------------------------------------------------------------------------------#
    # Compute: Exp. of crossproduct of factor scores given data (Theta)
    #-------------------------------------------------------------------------------#

    # Self-created function (see '4ComputeTheta.R').
    Theta_k <- comTheta(Beta_k, n_state, C_k, Lambda_k, n_fact)

    #===============================================================================#
    # Update: Loadings Lambda_k and unique variances Psi_k
    #===============================================================================#

    # Self-created function (see '5UpdateLoadingUnique.R')
    LambPsi<-updLambPsi(n_state, C_k, n_fact, Beta_k, Theta_k,
                        residualVariance,J)
    Lambda_k <- LambPsi$Lambda_k
    Psi_k <- LambPsi$Psi_k

    #-------------------------------------------------------------------------------#
    # Compute: Observed-data loglikelihood
    #-------------------------------------------------------------------------------#

    #*******************************************************************************#
    # NOTE: FIRST calculate the DMV (response probabilities) to be able to re-use
    # them because this takes quite some time.
    # THEN calculate the likelihood while reusing the response probabilities DMV.
    #*******************************************************************************#

    # Self-created function (see '6ComputeResponseprobSaveDMV.R').
    saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,nu_k)

    # Obtain the observed-data loglikelihood.
    logli <- c(rep(0,n_sub))
    for(i in 1:n_sub){
      for(k in 1:n_state){
        logli[i] <-logli[i]+pi_k[[k]]*saveDMV[k,i]
      }
    }
    total_logl <-sum(log(logli))


    #*******************************************************************************#
    # NOTE: Get all relevant parameters in one list to continue with the ones that
    # belong to the best start sets according to the loglikelihood.
    #*******************************************************************************#

    # Make DMV a list for convenient storage.
    DMV_list <- rep(list(NA),n_state)
    for(k in 1:n_state){ DMV_list[[k]] <- saveDMV[k,]}

    AllParameters <-list(pi_k,              #state proportions
                        nu_k,              #state-specific intercepts
                        Lambda_k,          #state-specific loading matrices
                        lapply(Psi_k,diag),#state-specific unique variances
                        total_logl,        #loglikelihood value
                        DMV_list,          #state-specific resp. probabilities
                        C_k)               #sample covariance matrix

    # Store all mclust parameters
    MclustResults[[mcluststarts]] <- AllParameters

  }

  # Evaluate which mclust works best.
  # Extract the likelihood values.
  loglikMulti <- as.data.frame(matrix(unlist(lapply(MclustResults,
                                                    function(x) {x[[5]]})),ncol=1))
  row.names(loglikMulti) <- 1:n_mclust

  # Obtain the number of the best mclust starts and use them for creating random deviations.
  best_mclust <- order(loglikMulti,decreasing = T)[1]
  ini_mclust <- ini_mclust_list[[best_mclust]]
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Initialize n_starts*10 random deviations from best mclust assignments 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  
  ini_mclust_random <- matrix(NA,ncol = (n_starts*10), nrow =n_sub )
  for(multistart in 1:(n_starts*10)){
    if(multistart == 1){
        change_ini_mclust <- ini_mclust
    }else{
      change_ini_mclust <-ini_mclust_random[,multistart-1]
    }
    change_ini_mclust[sample(1:n_sub,size = 0.30*n_sub)] <- sample(1:n_state,(0.30*n_sub),replace=TRUE)
    ini_mclust_random[,multistart] <- change_ini_mclust
  }
  

  if(n_starts>0){ #otherwise only mclust is used
  multistart <- NA #just because the CRAN check would otherwise 
  # say that there is no visible binding for global variable 'multistart'.
    MultistartResults1 <- foreach(multistart = 1: (n_starts*10),
                                  .packages = c("doParallel",
                                              "mclust",
                                              "NPflow"),
                                  .export = c("RandVec",
                                            #"ini_mclust",
                                            "initializeStep1",
                                            "updExpMem",
                                            "comBetas",
                                            "comTheta",
                                            "updLambPsi",
                                            "DMV"),
                                  .verbose = FALSE)%dopar%{
    
    ini_mclust_specific <- c(ini_mclust_random[,multistart])
    # Self-created function (see '1InitializeEM.R').
    InitialValues <- initializeStep1(x,n_sub,n_state,
                                     n_fact,J,startval = "MCrandom",
                                     RandVec = RandVec,
                                     ini_mclust = ini_mclust,
                                     ini_mclust_specific = ini_mclust_specific)#random;

    # Extract all parameters that are going to be updated in the EM algorithm.
    z_ik<- InitialValues$z_ik         #expected state-membership-probabilities
    N_k<- InitialValues$N_k           #sample size per state
    pi_k<- InitialValues$pi_k         #state proportions
    nu_k<- InitialValues$nu_k         #state-specific intercepts
    C_k<- InitialValues$C_k           #sample covariance matrix
    Lambda_k<- InitialValues$Lambda_k #state-specific loading matrices
    Psi_k<- InitialValues$Psi_k       #state-specific unique variances

    #-------------------------------------------------------------------------------#
    # Compute: Regression Weights (Beta)
    #-------------------------------------------------------------------------------#

    # Self-created function (see '3ComputeBeta.R').
    Beta_k <- comBetas(Lambda_k, Psi_k, n_state,n_fact)

    #-------------------------------------------------------------------------------#
    # Compute: Exp. of crossproduct of factor scores given data (Theta)
    #-------------------------------------------------------------------------------#

    # Self-created function (see '4ComputeTheta.R')
    Theta_k <- comTheta(Beta_k, n_state, C_k, Lambda_k, n_fact)

    #===============================================================================#
    # Update: Loadings Lambda_k and unique variances Psi_k
    #===============================================================================#

    # Self-created function (see '5UpdateLoadingUnique.R')
    LambPsi<-updLambPsi(n_state, C_k, n_fact, Beta_k, Theta_k,
                        residualVariance,J)
    Lambda_k <- LambPsi$Lambda_k
    Psi_k <- LambPsi$Psi_k

    #-------------------------------------------------------------------------------#
    # Compute: Observed-data loglikelihood
    #-------------------------------------------------------------------------------#

    #*******************************************************************************#
    # NOTE: FIRST calculate the DMV (response probabilities) to be able to re-use
    # them because this takes quite some time.
    # THEN calculate the likelihood while reusing the response probabilities DMV.
    #*******************************************************************************#

    # Self-created function (see '6ComputeResponseprobSaveDMV.R').
    saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,nu_k)

    # Obtain the observed-data loglikelihood.
    logli <- c(rep(0,n_sub))
    for(i in 1:n_sub){
      for(k in 1:n_state){
        logli[i] <-logli[i]+pi_k[[k]]*saveDMV[k,i]
      }
    }
    total_logl <-sum(log(logli))

    #*******************************************************************************#
    # NOTE: Get all relevant parameters in one list to continue with the ones that
    # belong to the best start sets according to the loglikelihood.
    #*******************************************************************************#

    # Make DMV a list for convenient storage.
    DMV_list <- rep(list(NA),n_state)
    for(k in 1:n_state){ DMV_list[[k]] <- saveDMV[k,]}

    AllParameters <-list(pi_k,              #state proportions
                         nu_k,              #state-specific intercepts
                         Lambda_k,          #state-specific loading matrices
                         lapply(Psi_k,diag),#state-specific unique variances
                         total_logl,        #loglikelihood value
                         DMV_list,          #state-specific resp. probabilities
                         C_k)               #sample covariance matrix

    return(AllParameters)
   }
  }
  
  # Add the mclust solutions to the (existing multistart) results.
  for(i in 1:n_mclust){
    MultistartResults1[[(length(lengths(MultistartResults1))+1)]] <- MclustResults[[i]]
  }
  # Extract the likelihood values.
  loglikMulti <- as.data.frame(matrix(unlist(lapply(MultistartResults1,
                                                    function(x) {x[[5]]})),ncol = 1))
  row.names(loglikMulti) <- 1:((n_starts*10)+n_mclust) #plus n_mclust for the mclust starts

  # Obtain the number of the best starts.
  best <- order(loglikMulti,decreasing = T)[1:(n_starts)]
  stopCluster(cl)
  stopImplicitCluster()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Do n_initial_ite iterations for the best 10 percent of the start sets
  # (in parallel).
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Start parallelization.
  getDoParWorkers()
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  multistart2 <- NA #just because the CRAN check would otherwise 
  # say that there is no visible binding for global variable 'multistart2'.
  cat("\n")
  cat(paste("2.Iterating through the best 10 %",
            "of the start sets..."))
  MultistartResults2 <- foreach(multistart2=1:(n_starts),
                     .packages=c("doParallel",
                                 "NPflow"),
                     .export=c("updExpMem",
                               "comBetas",
                               "comTheta",
                               "updLambPsi",
                               "DMV"),
                     .verbose = FALSE)%dopar%{

   # Extract parameter values belonging to the best loglikelihood values.
   TakeOverResults <- MultistartResults1[[c(best)[multistart2]]]
   pi_k <- TakeOverResults[[1]]
   nu_k<- TakeOverResults[[2]]
   Lambda_k <- TakeOverResults[[3]]
   Psi_k<- TakeOverResults[[4]]

   # Needs to be a matrix again.
   Psi_k<-lapply(Psi_k,function(x) x*diag(J))
   total_logl <- TakeOverResults[[5]]
   DMV_list<- TakeOverResults[[6]]

   # From list to matrix.
   saveDMV <- matrix(NA,nrow=n_state,ncol=n_sub)
   for(k in 1: n_state){
     saveDMV[k,] <- DMV_list[[k]]
   }
   C_k <- TakeOverResults[[7]]

   AllParameters <-list(pi_k,              #state proportions
                        nu_k,              #state-specific intercepts
                        Lambda_k,          #state-specific loading matrices
                        lapply(Psi_k,diag),#state-specific unique variances
                        total_logl,        #loglikelihood value
                        DMV_list,          #state-specific resp. probabilities
                        C_k,               #sample covariance matrix
                        list(estimation))

   M_step_parameters <- list(Lambda_k,          #state-specific loading matrices
                             lapply(Psi_k,diag))#state-specific unique variances

   #*******************************************************************************#
   # NOTE: Set the loglikelihood equal to the one from the start set, the number of
   # iterations equal to zero, and the difference between loglikelihood values to 1.
   # The values will be updated in the loop.
   #*******************************************************************************#
   LL <- total_logl
   iteration <- 0
   differenceLL <- 1

   while(iteration<n_initial_ite){
    #the other line (below) would be necessary for calculating the hitrate
    #while((sum(abs(differenceLL)>1e-3, iteration < max_iterations)==2)){

     iteration <- iteration+1

     #*******************************************************************************#
     # NOTE: In the following, I numbered the main 'update' steps,
     # not the 'computation' steps that are only required to obtain
     # updates of the relevant parameters
     #*******************************************************************************#

     #===============================================================================#
     # 1: Update: post. state-membership probabilities z_ig for
     #===============================================================================#

     # Self-created function (see '2UpdatePosteriorProb.R')
     z_ik <- updExpMem(Lambda_k, Psi_k, n_state, C_k, n_fact, DMV=saveDMV,n_sub,pi_k)

     #===============================================================================#
     # 2: Update: state-sp. samplesize N_k, state proportions pi_k,
     # and intercepts nu_k
     #===============================================================================#

     # Number of observations per state; kx1 vector.
     N_k <- lapply(z_ik,sum)

     # State proportions; kx1 vector.
     pi_k <- lapply(z_ik,mean)

     #-------------------------------------------------------------------------------#
     # Compute: Sample covariance matrix
     #-------------------------------------------------------------------------------#

     C_k <- rep(list(NA),n_state)
     for(sc in 1:n_state){
       #this is an existing function to obtain the weighted cov matrix
       SaveCov <- cov.wt(x,z_ik[[sc]],method = 'ML',center = T )
       C_k[[sc]] <- SaveCov$cov
       nu_k[[sc]] <- SaveCov$center #Here we also obtain nu_k
     }

     m_step <- 0
     SumParameterChange <- 100
     while(sum(m_step<n_m_step,SumParameterChange>m_step_tolerance)==2){
       m_step <- m_step+1
       #-------------------------------------------------------------------------------#
       # Compute: Regression Weights (Beta)
       #-------------------------------------------------------------------------------#

       # Self-created function (see '3ComputeBeta.R').
       Beta_k <- comBetas(Lambda_k, Psi_k, n_state,n_fact)

       #-------------------------------------------------------------------------------#
       # Compute: Exp. of crossproduct of factor scores given data (Theta)
       #-------------------------------------------------------------------------------#

       # Self-created function (see '4ComputeTheta.R')
       Theta_k <- comTheta(Beta_k, n_state, C_k, Lambda_k, n_fact)

       #===============================================================================#
       # 3: Update: Loadings Lambda_k and unique variances Psi_k
       #===============================================================================#

       # Self-created function (see '5UpdateLoadingUnique.R')
       LambPsi<-updLambPsi(n_state, C_k, n_fact, Beta_k, Theta_k,
                           residualVariance,J)
       Lambda_k <- LambPsi$Lambda_k
       Psi_k <- LambPsi$Psi_k


       #*******************************************************************************#
       # Evaluate intermediate parameter change to determine M-step iterations
       #*******************************************************************************#
       M_step_parametersN <- list(Lambda_k,          #state-specific loading matrices
                                  lapply(Psi_k,diag))#state-specific unique variances

       sumParChangeOld <- c(unlist(M_step_parameters))
       sumParChangeNew <- c(unlist(M_step_parametersN))

       #*******************************************************************************#
       # NOTE: The following ifelse replacement is done if a loading is ~zero. Otherwise
       # one would divide by 0. This is only relevant to calculate the sum of the
       # relative change. The parameters themselves are not constrained and can be zero.
       #*******************************************************************************#
       sumParChangeOld <- ifelse(abs(sumParChangeOld)<
                                   0.00000000001,0.00000000001,sumParChangeOld)
       sumParChangeNew <- ifelse(abs(sumParChangeNew)<
                                   0.00000000001,0.00000000001,sumParChangeNew)
       SumParameterChange <-sum(abs((sumParChangeOld-sumParChangeNew)/sumParChangeOld))
       M_step_parameters <- M_step_parametersN

     }

     #-------------------------------------------------------------------------------#
     # Compute: Observed-data loglikelihood
     #-------------------------------------------------------------------------------#

     #*******************************************************************************#
     # NOTE: FIRST calculate the DMV (response probabilities) to be able to re-use
     # them because this takes quite some time.
     # THEN calculate the likelihood while reusing the response probabilities DMV.
     #*******************************************************************************#

     # Self-created function (see '6ComputeResponseprobSaveDMV.R').
     saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,nu_k)

     # Obtain the observed-data loglikelihood.
     logli <- c(rep(0,n_sub))
     for(i in 1:n_sub){
       for(k in 1:n_state){
         logli[i] <-logli[i]+pi_k[[k]]*saveDMV[k,i]
       }
     }
     total_logl <-sum(log(logli))

     # Make DMV a list for convenient storage.
     DMV_list <- rep(list(NA),n_state)
     for(k in 1:n_state){ DMV_list[[k]] <- saveDMV[k,]}

     #-------------------------------------------------------------------------------#
     # Compute: Convergence
     #-------------------------------------------------------------------------------#
     differenceLL <- LL-total_logl
     LL <- total_logl
     estimation[iteration,] <- c(iteration,LL,LambPsi$act_constraints)
     AllParametersN <-list(pi_k,              #state proportions
                           nu_k,              #state-specific intercepts
                           Lambda_k,          #state-specific loading matrices
                           lapply(Psi_k,diag),#state-specific unique variances
                           total_logl,        #loglikelihood value
                           DMV_list,          #state-specific resp. probabilities
                           C_k,               #sample covariance matrix
                           list(estimation),
                           iteration)

     AllParameters <- AllParametersN
   }
   return(AllParameters) #end multistart loop after this one
  }
  stopCluster(cl)
  stopImplicitCluster()

  # Extract the likelihood values
  loglikMulti <- as.data.frame(matrix(unlist(lapply(MultistartResults2,
                                                    function(x) {x[[5]]})),ncol=1))
  row.names(loglikMulti) <- 1:(n_starts)

  # Obtain the number of the best starts
  best <- order(loglikMulti,decreasing = T)[1]

  equalLogli <- round(loglikMulti[,1])
  equalLogli <- unlist(equalLogli)
  hitrate <- (length(equalLogli[duplicated(equalLogli)]))/(length(equalLogli)-1)
  #hitrate <- 1- length(unique(equalLogli))/length(equalLogli)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Continue looping through the start set with the best loglikelihood value
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  # Extract parameter values belonging to the best loglikelihood values
  TakeOverResults <- MultistartResults2[[c(best)]]
  pi_k <- TakeOverResults[[1]]
  nu_k<- TakeOverResults[[2]]
  Lambda_k <- TakeOverResults[[3]]
  Psi_k<- TakeOverResults[[4]]
  total_logl <- TakeOverResults[[5]]
  DMV_list <- TakeOverResults[[6]]
  C_k <- TakeOverResults[[7]]
  estimation <- TakeOverResults[[8]][[1]]
  iteration <-TakeOverResults[[9]]

  # Needs to be a matrix again.
  Psi_k<-lapply(Psi_k,function(x) x*diag(J))


  # From list to matrix.
  saveDMV <- matrix(NA,nrow=n_state,ncol=n_sub)
  for(k in 1: n_state){
    saveDMV[k,] <- DMV_list[[k]]
  }


  AllResultsBestSet <- list()
  resultNumber <- 0
  AllParameters <-list(pi_k,              #state proportions
                       nu_k,              #state-specific intercepts
                       Lambda_k,          #state-specific loading matrices
                       lapply(Psi_k,diag),#state-specific unique variances
                       total_logl,        #loglikelihood value
                       DMV_list,          #state-specific resp. probabilities
                       C_k)               #sample covariance matrix

  M_step_parameters <- list(Lambda_k,          #state-specific loading matrices
                            lapply(Psi_k,diag))#state-specific unique variances

  #*******************************************************************************#
  # NOTE: Set the loglikelihood equal to the one from the start set, the number of
  # iterations equal to n_initial_ite, and the difference between loglikelihood
  # values to 1. The values will be updated in the loop.
  #*******************************************************************************#
  LL <- total_logl
  #iteration <- n_initial_ite
  differenceLL <- 1

  #*******************************************************************************#
  # NOTE: Give some value to the sum of parameter change to start the loop with
  # (will be updated in the loop).
  #*******************************************************************************#
  SumParameterChange <- 100
  cat("\n")
  cat(paste("3.Analyzing data with the best start set..."))

  LLcorrect <- TRUE

  while((sum(abs(differenceLL)>em_tolerance,SumParameterChange>em_tolerance,
  iteration < max_iterations)==3)){
    iteration <- iteration+1
    resultNumber <-resultNumber+1
    #if(iteration==max_iterations) stop("maximum number of iterations reached without convergence")

    #*******************************************************************************#
    # NOTE: In the following, I numbered the main 'update' steps,
    # not the 'computation' steps that are only required to obtain
    # updates of the relevant parameters
    #*******************************************************************************#

    #===============================================================================#
    # 1: Update: post. state-membership probabilities z_ig for
    #===============================================================================#

    # Self-created function (see '2UpdatePosteriorProb.R')
    z_ik <- updExpMem(Lambda_k, Psi_k, n_state, C_k, n_fact, DMV=saveDMV,n_sub,pi_k)

    #===============================================================================#
    # 2: Update: state-sp. samplesize N_k, state proportions pi_k,
    # and intercepts nu_k
    #===============================================================================#

    # Number of observations per state; kx1 vector.
    N_k <- lapply(z_ik,sum)

    # State proportions; kx1 vector.
    pi_k <- lapply(z_ik,mean)

    #-------------------------------------------------------------------------------#
    # Compute: Sample covariance matrix
    #-------------------------------------------------------------------------------#

    C_k <- rep(list(NA),n_state)
    for(sc in 1:n_state){
      #this is an existing function to obtain the weighted cov matrix
      SaveCov <- cov.wt(x,z_ik[[sc]],method = 'ML',center = T )
      C_k[[sc]] <- SaveCov$cov
      nu_k[[sc]] <- SaveCov$center #Here we also obtain nu_k
    }

    m_step <- 0
    SumParameterChange <- 100
    while(sum(m_step<n_m_step,SumParameterChange>m_step_tolerance)==2){
      m_step <- m_step+1

      #-------------------------------------------------------------------------------#
      # Compute: Regression Weights (Beta)
      #-------------------------------------------------------------------------------#

      # Self-created function (see '3ComputeBeta.R').
      Beta_k <- comBetas(Lambda_k, Psi_k, n_state,n_fact)

      #-------------------------------------------------------------------------------#
      # Compute: Exp. of crossproduct of factor scores given data (Theta)
      #-------------------------------------------------------------------------------#

      # Self-created function (see '4ComputeTheta.R')
      Theta_k <- comTheta(Beta_k, n_state, C_k, Lambda_k, n_fact)

      #===============================================================================#
      # 3: Update: Loadings Lambda_k and unique variances Psi_k
      #===============================================================================#

      # Self-created function (see '5UpdateLoadingUnique.R')
      LambPsi<-updLambPsi(n_state, C_k, n_fact, Beta_k, Theta_k,
                          residualVariance,J)
      Lambda_k <- LambPsi$Lambda_k
      Psi_k <- LambPsi$Psi_k


      #*******************************************************************************#
      # Evaluate intermediate parameter change to determine M-step iterations
      #*******************************************************************************#
      M_step_parametersN <- list(Lambda_k,          #state-specific loading matrices
                                 lapply(Psi_k,diag))#state-specific unique variances

      sumParChangeOld <- c(unlist(M_step_parameters))
      sumParChangeNew <- c(unlist(M_step_parametersN))

      #*******************************************************************************#
      # NOTE: The following ifelse replacement is done if a loading is ~zero. Otherwise
      # one would divide by 0. This is only relevant to calculate the sum of the
      # relative change. The parameters themselves are not constrained and can be zero.
      #*******************************************************************************#
      sumParChangeOld <- ifelse(abs(sumParChangeOld)<
                                  0.00000000001,0.00000000001,sumParChangeOld)
      sumParChangeNew <- ifelse(abs(sumParChangeNew)<
                                  0.00000000001,0.00000000001,sumParChangeNew)
      SumParameterChange <-sum(abs((sumParChangeOld-sumParChangeNew)/sumParChangeOld))
      M_step_parameters <- M_step_parametersN

    }
    #-------------------------------------------------------------------------------#
    # Compute: Observed-data loglikelihood
    #-------------------------------------------------------------------------------#

    #*******************************************************************************#
    # NOTE: FIRST calculate the DMV (response probabilities) to be able to re-use
    # them because this takes quite some time.
    # THEN calculate the likelihood while reusing the response probabilities DMV.
    #*******************************************************************************#

    # Self-created function (see '6ComputeResponseprobSaveDMV.R').
    saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,nu_k)

    # Obtain the observed-data loglikelihood.
    logli <- c(rep(0,n_sub))
    for(i in 1:n_sub){
      for(k in 1:n_state){
        logli[i] <-logli[i]+pi_k[[k]]*saveDMV[k,i]
      }
    }
    total_logl <-sum(log(logli))

    # Make DMV a list for convenient storage.
    DMV_list <- rep(list(NA),n_state)
    for(k in 1:n_state){ DMV_list[[k]] <- saveDMV[k,]}

    #-------------------------------------------------------------------------------#
    # Compute: Convergence
    #-------------------------------------------------------------------------------#
    differenceLL <- LL-total_logl
    LLcorrect <- total_logl>LL
    LL <- total_logl

    AllParametersN <- list(pi_k,              #state proportions
                           nu_k,              #state-specific intercepts
                           Lambda_k,          #state-specific loading matrices
                           lapply(Psi_k,diag),#state-specific unique variances
                           total_logl,        #loglikelihood value
                           DMV_list,          #state-specific resp. probabilities
                           C_k)               #sample covariance matrix


    sumParChangeOld <- c(unlist(AllParameters[1:4]))
    sumParChangeNew <- c(unlist(AllParametersN[1:4]))

    #*******************************************************************************#
    # NOTE: The following ifelse replacement is done if a loading is ~zero. Otherwise
    # one would divide by 0. This is only relevant to calculate the sum of the
    # relative change. The parameters themselves are not constrained and can be zero.
    #*******************************************************************************#
    sumParChangeOld <- ifelse(abs(sumParChangeOld)<
                                0.00000000001,0.00000000001,sumParChangeOld)
    sumParChangeNew <- ifelse(abs(sumParChangeNew)<
                                0.00000000001,0.00000000001,sumParChangeNew)

    SumParameterChange <-sum(abs((sumParChangeOld-sumParChangeNew)/sumParChangeOld))

    # Store the loglikelihood and activated constraints per iteration.
    estimation[iteration,] <- c(iteration,LL,LambPsi$act_constraints)

    AllParameters <- AllParametersN

    AllResultsBestSet[[resultNumber]] <-AllParameters

  }

  #
  if(iteration == max_iterations){
     warning("maximum number of iterations reached without convergence")
  }
  
  #-------------------------------------------------------------------------------#
  # Obtain the BIC values.
  #-------------------------------------------------------------------------------#

  # Calculate the number of parameters.
  R_k <- c(NA)
  for(k in 1:n_state){
    intercepts <- J
    variances <- J
    factorloadings <- J*(n_fact[k])
    covariances <- 0
    R_k[k] <- sum(intercepts,variances,factorloadings,covariances)
  }
  # Sum across states, subtract activated constraints, and add state proportions.
  R_T <- sum(R_k)-c(estimation[iteration,3])+n_state-1
  R_T <- as.numeric(R_T)

  #BIC_N <- -2*LL + R_T *log(n_cases)
  BIC_T <- -2*LL + R_T *log(n_sub)


  cat("\n")
  requiredTime <- as.numeric((proc.time() - ptm)[3])
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                       Step 2: Classification Error
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Obtain the final posterior state-membership probabilities.
  z_ik <- updExpMem(Lambda_k, Psi_k, n_state, C_k, n_fact, DMV=saveDMV,n_sub,pi_k)

  posterior_data <- matrix(c(NA),ncol=n_state,nrow = nrow(x))
  for(k in 1:n_state){
    posterior_data[,k] <- z_ik[[k]]
  }
  
  #--------------------------------------------------------------------------------#
  #                reorder parameters based on the size of the states
  #--------------------------------------------------------------------------------#
  orderStates <- c()
for(i in 1:n_state){
  orderStates <- c(orderStates,pi_k[[i]])
}

orderStates <- order(orderStates,decreasing = TRUE)

n_fact <- n_fact[orderStates]

pi_k_reorder <-pi_k
Lambda_k_reorder <- Lambda_k
Psi_k_reorder <- Psi_k
nu_k_reorder <- nu_k
posterior_data_reorder <- posterior_data
C_k_reorder <- C_k
for(i in 1:n_state){
  pi_k_reorder[[i]] <- pi_k[[orderStates[i]]]
  Lambda_k_reorder[[i]] <- Lambda_k[[orderStates[i]]]
  Psi_k_reorder[[i]] <- Psi_k[[orderStates[i]]]
  nu_k_reorder[[i]] <- nu_k[[orderStates[i]]]
  posterior_data_reorder[,i] <- posterior_data[,orderStates[i]]
  C_k_reorder[[i]] <-C_k[[orderStates[i]]]
}
pi_k <- pi_k_reorder
Lambda_k <- Lambda_k_reorder
Psi_k <- Psi_k_reorder
nu_k <- nu_k_reorder
posterior_data <- posterior_data_reorder
C_k <- C_k_reorder

  #--------------------------------------------------------------------------------#
  #                            continue with step 2
  #--------------------------------------------------------------------------------#

  # Obtain the modal state assignments.
  modal_data <- max.col(posterior_data)

  Posteriors <-cbind.data.frame(modal_data,posterior_data)
  colnames(Posteriors) <- c("Modal", paste("State",1:n_state,sep=""))

  ModalClassificationTable <- matrix(NA,ncol=n_state,nrow=n_state)
  for(i in 1:n_state){
    for(j in 1:n_state){
      ModalClassificationTable[j,i] <- sum((Posteriors[Posteriors$Modal==i,j+1]))
    }
  }

  # Calculate probabilities based on counts.

  W_mod <-ModalClassificationTable/rowSums(ModalClassificationTable)
  
  
  #-------------------------------------------------------------------------------#
  # Obtain R-squared entropy.
  #-------------------------------------------------------------------------------#
  if(n_state>1){
probVector <-c(NA)
  entropy <- function(p) sum(-p * log(p))
  for(i in 1:n_state){
    probVector[i] <- pi_k[[i]]
  }
  error_prior <- entropy(probVector) # Class proportions
  posteriors <-  Posteriors[,-1] #just for internal usage
  posteriors[(posteriors==0)] <- 1e-21
  error_post <- mean(apply(posteriors, 1, entropy))
  R2_entropy <- (error_prior - error_post) / error_prior
  }else{
    R2_entropy <-1
  }
  

  #-------------------------------------------------------------------------------#
  # Obtain standardized loadings and proportions of unique variances.
  #-------------------------------------------------------------------------------#
  x$State <- Posteriors[,1]
  #get items' standard deviation per state
  SDList <- list()
  for(i in 1:n_state){
    SDMatrix <- matrix(NA,J)
    rownames(SDMatrix) <- indicators
    for(item in 1:J){
      SDMatrix[item,]<- sd(unlist(x[x$State==i,indicators[item]]),na.rm = T)
    }
    SDList[[i]] <- SDMatrix
  }

  #get items' standard deviation across states
  SDList2 <- list()
  for(i in 1:n_state){
    SDMatrix <- matrix(NA,J)
    rownames(SDMatrix) <- indicators
    for(item in 1:J){
      SDMatrix[item,]<- sd(unlist(x[,indicators[item]]),na.rm = T)
    }
    SDList2[[i]] <- SDMatrix
  }

  #standardize loadings per state for better within-state comparison
  standLambda <- Lambda_k
  for(i in 1:n_state){
    for(fact in 1:n_fact[i]){
      for(j in 1:J){
        if(SDList[[i]][j]!=0){
          standLambda[[i]][j,fact] <- Lambda_k[[i]][j,fact]/SDList[[i]][j] 
        }
      }
    }
  }

  #standardize loadings across states for better between-state comparison
  standLambda2 <- Lambda_k
  for(i in 1:n_state){
    for(fact in 1:n_fact[i]){
      for(j in 1:J){
        if(SDList2[[i]][j]!=0){
          standLambda2[[i]][j,fact] <- Lambda_k[[i]][j,fact]/SDList2[[i]][j] 
        }
      }
    }
  }
  
  uniqueVariances <- lapply(Psi_k,diag)
  
  #proportions of unique variance per state for better within-state comparison
  standPsi <- uniqueVariances
  for(i in 1:n_state){
    for(j in 1:J){
        if(SDList[[i]][j]!=0){
          standPsi[[i]][j] <- uniqueVariances[[i]][j]/(SDList[[i]][j]^2) 
        }
      }
  }

  #proportions of unique variance across states for better between-state comparison
  
  standPsi2 <- uniqueVariances
  for(i in 1:n_state){
    for(j in 1:J){
        if(SDList2[[i]][j]!=0){
          standPsi2[[i]][j] <- uniqueVariances[[i]][j]/(SDList2[[i]][j]^2) 
        }
      }
  }

  #-------------------------------------------------------------------------------#
  # Obtain rotated solutions
  #-------------------------------------------------------------------------------#
  #Lambda_k_obli <- lapply(Lambda_k, function(x) GPFoblq(x, method = "oblimin", normalize = FALSE)$loadings[])
  #standLambda_obli <- lapply(standLambda, function(x) GPFoblq(x, method = "oblimin", normalize = FALSE)$loadings[])
  #standLambda2_obli <- lapply(standLambda2, function(x) GPFoblq(x, method = "oblimin", normalize = FALSE)$loadings[])
  #correlations_obli <- lapply(standLambda, function(x) GPFoblq(x, method = "oblimin", normalize = FALSE)$Phi[])
  Lambda_obli <- Lambda_k
  standLambda_obli <- standLambda
  standLambda2_obli <- standLambda2
  correlations_obli <- list()
  correlations_obli_unstandardized <- list()

  for(i in 1:n_state){
    if(n_fact[i]>1){
      rotationResults <- suppressWarnings(GPFoblq(standLambda2[[i]], method = "oblimin", normalize = FALSE))#between
      rotationResultsunstandardized <-suppressWarnings(GPFoblq(Lambda_k[[i]], method = "oblimin", normalize = TRUE))
      standLambda_obli[[i]] <- suppressWarnings(GPFoblq(standLambda[[i]], method = "oblimin", normalize = FALSE)$loadings[])#within
      correlations_obli[[i]] <- rotationResults$Phi[]#between
      correlations_obli_unstandardized[[i]] <- rotationResultsunstandardized$Phi[]#unstandardized; usually, for unstandardized loadings, only normalizing works, which can lead to small deviations in correlations
      standLambda2_obli[[i]] <- rotationResults$loadings[]#between
      Lambda_obli[[i]] <- rotationResultsunstandardized$loadings[]#unstandardized
    }else{
      correlations_obli[[i]] <- 1
      correlations_obli_unstandardized[[i]]<- 1 #usually, for unstandardized loadings, only normalizing works, which can lead to small deviations in correlations
    }
  }

  #sometimes the rotation results in warnings (thus, the results are not reliable)
  #for now we only check this for the between-state standardized loadings and unstandardized loadings because the within-state standardized loadings are not reported
  #standardized
  warning_stand_loadings <- c()

  for(i in 1:n_state){
    if(n_fact[i]>1){

          WarningStandardized <- suppressWarnings(tryCatch(GPFoblq(standLambda2[[i]], method = "oblimin", normalize = FALSE),
            warning = function(w) return(list(GPFoblq(standLambda2[[i]], method = "oblimin", normalize = FALSE),w))))


    if(length(grep("simpleWarning", as.character(WarningStandardized[[2]])))>0){
      warning_stand_loadings <- c(warning_stand_loadings,1)
    }else{
      warning_stand_loadings <- c(warning_stand_loadings,0)
    }
   }
  }

  if(sum(warning_stand_loadings)>0){
    warningRotationStandardized <- c("Warning message: convergence for rotating loadings in at least one state was not obtained")
  }else{
    warningRotationStandardized <- c("no warning")
  }

  #unstandardized
  warning_loadings <- c()

  for(i in 1:n_state){
    if(n_fact[i]>1){

    WarningUnstandardized <- suppressWarnings(tryCatch(GPFoblq(Lambda_k[[i]], method = "oblimin", normalize = TRUE),
                                                    warning = function(w) return(list(GPFoblq(Lambda_k[[i]], method = "oblimin", normalize = TRUE),w))))


    if(length(grep("simpleWarning", as.character(WarningUnstandardized[[2]])))>0){
      warning_loadings <- c(warning_loadings,1)
    }else{
      warning_loadings <- c(warning_loadings,0)
    }
   }
  }
    
  if(sum(warning_loadings)>0){
    warningRotationUnstandardized <- c("Warning message: convergence for rotating loadings in at least one state was not obtained")
  }else{
    warningRotationUnstandardized <- c("no warning")
  }


  #-------------------------------------------------------------------------------#
  # Obtain explained variance per state and in total.
  #-------------------------------------------------------------------------------#
  #sum of squared loadings per state
  SSL_l_k <- lapply(standLambda, function(x) sum(x^2))
  Percent_expl_var_k <- lapply(SSL_l_k, function(x) x/J)
  #weighted by state proportion
  Percent_expl_var_k_w <- mapply('*',Percent_expl_var_k, pi_k)
  #total explained variance
  Percent_expl_var <- Reduce("+", Percent_expl_var_k_w)

  #Number of variances that are equal to zero
  NumberZeroVariance <- Reduce("+",lapply(SDList,function(x) sum(x^2==0)))

  if(NumberZeroVariance>0) warning("one or more states contain zero item variences and therefore, the explained variance per state is not interpretable")
 
  Psi_k = lapply(lapply(Psi_k,diag),round,rounding)
  Psi_k_st_w = lapply(standPsi,round,rounding)
  Psi_k_st_b = lapply(standPsi2,round,rounding)
  for(i in 1:n_state){
    names(Psi_k[[i]])<- indicators
    names(Psi_k_st_w[[i]])<- indicators
    names(Psi_k_st_b[[i]])<- indicators
  }

  #-------------------------------------------------------------------------------#
  # Transpose unique variances and intercepts
  #-------------------------------------------------------------------------------#
  for(i in 1:n_state){
    Psi_k[[i]] <- t(t(Psi_k[[i]]))
    nu_k[[i]] <- t(t(nu_k[[i]]))
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                    On exit
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  x <- x[,-c(ncol(x))]

  on.exit(stopImplicitCluster(), add=TRUE)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                    flip loadings and correlations
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #copy loadings to flipped versions
  Lambda_k_fl <- Lambda_k
  standLambda_fl <- standLambda
  standLambda2_fl <- standLambda2
  standLambda_obli_fl <- standLambda_obli
  standLambda2_obli_fl <- standLambda2_obli
  correlations_obli_fl <- correlations_obli
  Lambda_obli_fl <- Lambda_obli
  correlations_obli_unstandardized_fl <- correlations_obli_unstandardized

  #normal unrotated unstandardized
  for(i in 1:n_state){
    for(j in 1:n_fact[i]){
      normal <- sum(abs(Lambda_k[[i]][,j][(Lambda_k[[i]][,j])>0]))
      flipped <- sum(abs(Lambda_k[[i]][,j][(Lambda_k[[i]][,j])* -1 > 0]))
      if(flipped>normal){
        Lambda_k_fl[[i]][,j] <- Lambda_k[[i]][,j]*-1
      }
    }
  }

  #within Standardized
  for(i in 1:n_state){
    for(j in 1:n_fact[i]){
      normal <- sum(abs(standLambda[[i]][,j][(standLambda[[i]][,j])>0]))
      flipped <- sum(abs(standLambda[[i]][,j][(standLambda[[i]][,j])* -1 > 0]))
      if(flipped>normal){
        standLambda_fl[[i]][,j] <- standLambda[[i]][,j]*-1
      }
    }
  }

  #between Standardized
  for(i in 1:n_state){
    for(j in 1:n_fact[i]){
      normal <- sum(abs(standLambda2[[i]][,j][(standLambda2[[i]][,j])>0]))
      flipped <- sum(abs(standLambda2[[i]][,j][(standLambda2[[i]][,j])* -1 > 0]))
      if(flipped>normal){
        standLambda2_fl[[i]][,j] <- standLambda2[[i]][,j]*-1
      }
    }
  }
  #within standardized rotated
  for(i in 1:n_state){
    for(j in 1:n_fact[i]){
      normal <- sum(abs(standLambda_obli[[i]][,j][(standLambda_obli[[i]][,j])>0]))
      flipped <- sum(abs(standLambda_obli[[i]][,j][(standLambda_obli[[i]][,j])* -1 > 0]))
      if(flipped>normal){
        standLambda_obli_fl[[i]][,j] <- standLambda_obli[[i]][,j]*-1
      }
    }
  }

  #between Standardized rotated
  for(i in 1:n_state){
    for(j in 1:n_fact[i]){
      normal <- sum(abs(standLambda2_obli[[i]][,j][(standLambda2_obli[[i]][,j])>0]))
      flipped <- sum(abs(standLambda2_obli[[i]][,j][(standLambda2_obli[[i]][,j])* -1 > 0]))
      if(flipped>normal){
        standLambda2_obli_fl[[i]][,j] <- standLambda2_obli[[i]][,j]*-1
        if(n_fact[i]>1){
          correlations_obli_fl[[i]][,j] <-  correlations_obli_fl[[i]][,j]*-1 
          correlations_obli_fl[[i]][j,] <-  correlations_obli_fl[[i]][j,]*-1 
        }
      }
    }
  }

  #Unstandardized unrotated
  for(i in 1:n_state){
    for(j in 1:n_fact[i]){
      normal <- sum(abs(Lambda_obli_fl[[i]][,j][(Lambda_obli_fl[[i]][,j])>0]))
      flipped <- sum(abs(Lambda_obli_fl[[i]][,j][(Lambda_obli_fl[[i]][,j])* -1 > 0]))
      if(flipped>normal){
        Lambda_obli_fl[[i]][,j] <- Lambda_obli_fl[[i]][,j]*-1
        if(n_fact[i]>1){
          correlations_obli_unstandardized_fl[[i]][,j] <-  correlations_obli_unstandardized_fl[[i]][,j]*-1 
          correlations_obli_unstandardized_fl[[i]][j,] <-  correlations_obli_unstandardized_fl[[i]][j,]*-1 
        }
      }
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                             get nice output
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
state_proportions <- table(Posteriors[,1])/n_sub
# for(i in 1:length(pi_k)){
#   state_proportions <- cbind(state_proportions,pi_k[[i]])
# }
# state_proportions <- as.data.frame(state_proportions)
# colnames(state_proportions) <- c(paste("S",rep(1:n_state),sep=""))
# state_proportions <- as.matrix(state_proportions)

factornames <- c()
for(i in 1:n_state){
  factornames <- c(factornames,1:n_fact[i])
}

loadings_w_obli_fl<- c()
for(i in 1:length(standLambda_obli_fl)){
  loadings_w_obli_fl<- cbind(loadings_w_obli_fl,standLambda_obli_fl[[i]])
}
loadings_w_obli_fl<- as.data.frame(loadings_w_obli_fl)
colnames(loadings_w_obli_fl) <- c(paste("S",rep(1:n_state,n_fact),"F", factornames,sep=""))


loadings_b_obli_fl<- c()
for(i in 1:length(standLambda2_obli_fl)){
  loadings_b_obli_fl<- cbind(loadings_b_obli_fl,standLambda2_obli_fl[[i]])
}
loadings_b_obli_fl<- as.data.frame(loadings_b_obli_fl)
colnames(loadings_b_obli_fl) <- c(paste("S",rep(1:n_state,n_fact),"F", factornames,sep=""))

intercepts <- c()
for(i in 1:length(nu_k)){
  intercepts <- cbind(intercepts,nu_k[[i]])
}
intercepts <- as.data.frame(intercepts)
colnames(intercepts) <- c(paste("S",rep(1:n_state),sep=""))

unique_variances <- c()
for(i in 1:length(Psi_k)){
  unique_variances <- cbind(unique_variances,Psi_k[[i]])
}
unique_variances <- as.data.frame(unique_variances)
colnames(unique_variances) <- c(paste("S",rep(1:n_state),sep=""))

names(correlations_obli) <- c(paste("S",rep(1:n_state),sep=""))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                    Return Step 1 and Step 2 Results
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  cat("\n")
  
  if(iteration<max_iterations){
    convergence <- 1
    cat("\n")
    cat(paste("Estimation converged after",iteration,"iterations."))
  }else{
    convergence <- 0
    cat("\n")
    cat("Maximum number of iterations reached without convergence.")
  }
  cat("\n")
  cat(paste("LL",round(LL,2),sep=" = "),"\n")
  cat("\n")
  cat("-------------------------------------------------------------")
  cat("\n")  
    output <- list(n_it = iteration,
              seconds = requiredTime,
              convergence = convergence,
              LL = round(LL,rounding),
              BIC = round(BIC_T,rounding),
              intercepts = round(intercepts, rounding),
              #loadings_w_obli = round(loadings_w_obli_fl, rounding),
              loadings_stand_obli = round(loadings_b_obli_fl, rounding),
              unique_variances = round(unique_variances, rounding),
              state_proportions = round(state_proportions, rounding),
              n_obs = n_sub,
              n_par = R_T,
              explained_variance = round(Percent_expl_var, rounding),
              n_state = n_state,
              n_fact = n_fact,
              #state_proportions_list = lapply(pi_k, round, rounding),
              intercepts_list = lapply(nu_k, round, rounding),
              loadings_list = lapply(Lambda_k_fl, round, rounding),
              #loadings_w_list = lapply(standLambda_fl, round, rounding),
              loadings_stand_list = lapply(standLambda2_fl, round, rounding),
              loadings_obli_list = lapply(Lambda_obli_fl,round,rounding),
              #loadings_w_obli_list = lapply(standLambda_obli_fl, round, rounding),
              loadings_stand_obli_list = lapply(standLambda2_obli_fl, round, rounding), 
              
              unique_variances_list = lapply(Psi_k, round, rounding),
              factor_correlations_stand_obli_list = lapply(correlations_obli_fl, round, rounding), 
              factor_correlations_obli_list = lapply(correlations_obli_unstandardized_fl, round, rounding), 
              #Psi_k_st_w = Psi_k_st_w,
              #Psi_k_st_b = Psi_k_st_b,
              activated_contraints = estimation[iteration,3],
              #standard_dev_k = SDList,
              #standard_dev = SDList2,
              classification_posteriors = Posteriors,
              classification_errors = round(ModalClassificationTable, rounding),
              classification_errors_prob = round(W_mod, rounding),
              R2_entropy = round(R2_entropy, rounding),
              warning_loadings = warningRotationUnstandardized,
              warning_loadings_stand = warningRotationStandardized,
              #sample_cov_matrix_list = C_k,
              raw_data = x
              #hitrate = hitrate,
              #loglies = loglikMulti,
              #estimation = estimation
              
              )
  class(output) = "lmfa_step1"

  if(modelselection==TRUE){
    all_estimated_models[[comparingmodels]] <-output
    names(all_estimated_models)[[comparingmodels]] <-paste("[",paste(n_fact,collapse = ""),"]",sep="")
  }
  } #closes comparing loop models
  
  if(modelselection==TRUE){
    class(all_estimated_models) <- "lmfa_modelselection"
    all_estimated_models
  }else{
    output
  }
  
}
