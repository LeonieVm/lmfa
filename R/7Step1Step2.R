#' Conducts steps 1 and 2 from the three-step estimation of CT-LMFA.
#'
#'
#'
#'
#'
#'
#' @param input_file The dataset (must be a dataframe and contain complete cases only).
#' @param variable_columns The variable names of the indicators (must be a vector of characters).
#' @param n_state The number of states that should be estimated (must be a single scalar).
#' @param n_fact The number of factors per state that should be estimated (must be a numeric vector of length n_state).
#' @param n_starts The number of random starts that should be used (must be a single scalar).
#' @param n_initial_ite The number of initial iterations for the best starts (must be a single scalar).
#' @param n_m_step The number of M-step iterations that should be used when parameters still change more than defined by the m_step_tolerance (must be a single scalar).
#' @param em_tolerance The convergence criterion for parameters and loglikelihood (must be a single scalar and smaller than m_step_tolerance).
#' @param m_step_tolerance The criterion for stopping the n_m_step M-step interations (must be a single scalar).
#' @param max_iterations The maximum number of iterations (must be a single scalar and larger than n_initial_ite).
#
#'
#' @return Returns the measurement model parameters, the proportional and
#' modal state assignments, and the classification errors.
#'
#' @examples
#' \dontrun{
#' fitStep1Step2 <- Step1Step2(input_file,variable_columns,n_state,
#'                       n_fact,n_starts=25,n_initial_ite=15,n_m_step=10,
#'                       em_tolerance=1e-6,m_step_tolerance=1e-3,max_iterations=500)
#' }
#' @export

Step1Step2 <- function(input_file,variable_columns,n_state,
                       n_fact,n_starts=25,n_initial_ite=15,n_m_step=10,
                       em_tolerance=1e-6,m_step_tolerance=1e-3,max_iterations=500){

  if(missing(input_file)) stop("argument input_file is missing, with no default")
  if(missing(variable_columns)) stop("argument variable_columns is missing, with no default")
  #if(missing(id_column)) stop("argument id_column is missing, with no default")
  if(missing(n_state)) stop("argument n_state is missing, with no default")
  if(missing(n_fact)) stop("argument n_fact is missing, with no default")

  if(!is.data.frame(input_file)) stop("input_file must be a dataframe")
  if(!is.character(variable_columns)) stop("variable_columns must be a vector of characters")
  #if(!is.character(id_column)) stop("id_column must be a single character")
  #if(length(id_column)>1) stop("id_column must be a single character")
  if(!is.numeric(n_state)) stop("n_state must be a single scalar")
  if(length(n_state)>1) stop("n_state must be a single scalar")
  if(!is.numeric(n_fact)) stop("n_state must be a numeric vector")
  if(length(n_fact)!=n_state) stop("n_fact must be of length n_state")

  if(!is.numeric(n_initial_ite)) stop("n_initial_ite must be a single scalar")
  if(length(n_initial_ite)>1) stop("n_initial_ite must be a single scalar")
  if(!is.numeric(n_m_step)) stop("n_m_step must be a single scalar")
  if(length(n_m_step)>1) stop("n_m_step must be a single scalar")
  if(!is.numeric(em_tolerance)) stop("em_tolerance must be a single scalar")
  if(length(em_tolerance)>1) stop("em_tolerance must be a single scalar")
  if(!is.numeric(m_step_tolerance)) stop("m_step_tolerance must be a single scalar")
  if(length(m_step_tolerance)>1) stop("m_step_tolerance must be a single scalar")
  if(!is.numeric(max_iterations)) stop("max_iterations must be a single scalar")
  if(length(max_iterations)>1) stop("max_iterations must be a single scalar")

  if(sum(complete.cases(input_file[,variable_columns])==FALSE)>0) stop("input_file must contain complete cases with regard to the indicators only")
  if(max_iterations <= n_initial_ite) stop("max_iterations must be larger than n_initial_ite")
  if(em_tolerance >= m_step_tolerance) stop("em_tolerance must be smaller than m_step_tolerance")

  ptm <- proc.time()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Input: b) defined in the package.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  # Perhaps set a seed in order to replicate results (e.g., same starting values).
  #set.seed(16)

  # Obtain the columns with the variables
  x <- input_file[,variable_columns]
  x <- as.data.frame(x)

  # Number of observations.
  #*******************************************************************************#
  # NOTE: Usualy, in mixture factor analysis, this would be the number
  # of subjects but here this number represents the independently treated total
  # number of observations.
  #*******************************************************************************#
  n_sub <- nrow(x)

  # Number of items.
  J <- ncol(x)

  # Number of cases.
  #n_cases <- length(unique(input_file[,id_column]))

  # List of multistart procedure results.
  MultistartResults1 <- rep(list(list(NA)),n_starts*10)
  MultistartResults2 <- rep(list(list(NA)),n_starts)

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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Initialize: a) n_starts*10 of random partitions (in parallel).
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

  ini_mclust <- Mclust(x, G =n_state,verbose=FALSE)
  ini_mclust <- ini_mclust$classification

  if(n_starts>0){ #otherwise only mclust is used
    MultistartResults1 <- foreach(multistart=1: (n_starts*10),
                                  .packages=c("doParallel",
                                              "mclust",
                                              "NPflow"),
                                  .export=c("RandVec",
                                            #"ini_mclust",
                                            "initializeStep1",
                                            "updExpMem",
                                            "comBetas",
                                            "comTheta",
                                            "updLambPsi",
                                            "DMV"),
                                  .verbose = FALSE)%dopar%{

    # Self-created function (see '1InitializeEM.R').
    InitialValues <- initializeStep1(x,n_sub,n_state,
                                     n_fact,J,startval="MCrandom",RandVec=RandVec,ini_mclust=ini_mclust)#random;

    # Extract all parameters that are going to be updated in the EM algorithm.
    z_ik<- InitialValues$z_ik         #expected state-membership-probabilities
    N_k<- InitialValues$N_k           #sample size per state
    pi_k<- InitialValues$pi_k         #state proportions
    mu_k<- InitialValues$mu_k         #state-specific intercepts
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
    saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,mu_k)

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
    # belong to the best startsets according to the loglikelihood.
    #*******************************************************************************#

    # Make DMV a list for convenient storage.
    DMV_list <- rep(list(NA),n_state)
    for(k in 1:n_state){ DMV_list[[k]] <- saveDMV[k,]}

    AllParameters <-list(pi_k,              #state proportions
                         mu_k,              #state-specific intercepts
                         Lambda_k,          #state-specific loading matrices
                         lapply(Psi_k,diag),#state-specific unique variances
                         total_logl,        #loglikelihood value
                         DMV_list,          #state-specific resp. probabilities
                         C_k)               #sample covariance matrix

    return(AllParameters)
   }
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Initialize: b) once based on mclust (in parallel).
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  # Self-created function (see '1InitializeEM.R').
  InitialValues <- initializeStep1(x,n_sub,n_state,n_fact,
                                   J,startval="mclust",
                                   RandVec=RandVec,
                                   ini_mclust=ini_mclust)#mclust;


  # Extract all parameters that are going to be updated in the EM algorithm.
  z_ik<- InitialValues$z_ik         #expected state-membership-probabilities
  N_k<- InitialValues$N_k           #sample size per state
  pi_k<- InitialValues$pi_k         #state proportions
  mu_k<- InitialValues$mu_k         #state-specific intercepts
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
  saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,mu_k)

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
  # belong to the best startsets according to the loglikelihood.
  #*******************************************************************************#

  # Make DMV a list for convenient storage.
  DMV_list <- rep(list(NA),n_state)
  for(k in 1:n_state){ DMV_list[[k]] <- saveDMV[k,]}

  AllParameters <-list(pi_k,              #state proportions
                       mu_k,              #state-specific intercepts
                       Lambda_k,          #state-specific loading matrices
                       lapply(Psi_k,diag),#state-specific unique variances
                       total_logl,        #loglikelihood value
                       DMV_list,          #state-specific resp. probabilities
                       C_k)               #sample covariance matrix

  # Add the mclust solution to the (existing multistart) results.
  MultistartResults1[[(length(lengths(MultistartResults1))+1)]] <- AllParameters

  # Extract the likelihood values.
  loglikMulti <- as.data.frame(matrix(unlist(lapply(MultistartResults1,
                                                    function(x) {x[[5]]})),ncol=1))
  row.names(loglikMulti) <- 1:((n_starts*10)+1)

  # Obtain the number of the best starts.
  best <- order(loglikMulti,decreasing = T)[1:(n_starts+1)]
  stopCluster(cl)
  stopImplicitCluster()
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Do n_initial_ite iterations for the best 10 percent of the startsets
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
            "of the startsets..."))
  MultistartResults2 <- foreach(multistart2=1:(n_starts+1),
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
   mu_k<- TakeOverResults[[2]]
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
                        mu_k,              #state-specific intercepts
                        Lambda_k,          #state-specific loading matrices
                        lapply(Psi_k,diag),#state-specific unique variances
                        total_logl,        #loglikelihood value
                        DMV_list,          #state-specific resp. probabilities
                        C_k,               #sample covariance matrix
                        list(estimation))

   M_step_parameters <- list(Lambda_k,          #state-specific loading matrices
                             lapply(Psi_k,diag))#state-specific unique variances

   #*******************************************************************************#
   # NOTE: Set the loglikelihood equal to the one from the startset, the number of
   # iterations equal to zero, and the difference between loglikelihood values to 1.
   # The values will be updated in the loop.
   #*******************************************************************************#
   LL <- total_logl
   iteration <- 0
   differenceLL <- 1

   while(iteration<n_initial_ite){
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
     # and intercepts mu_k
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
       mu_k[[sc]] <- SaveCov$center #Here we also obtain mu_k
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
     saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,mu_k)

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
                           mu_k,              #state-specific intercepts
                           Lambda_k,          #state-specific loading matrices
                           lapply(Psi_k,diag),#state-specific unique variances
                           total_logl,        #loglikelihood value
                           DMV_list,          #state-specific resp. probabilities
                           C_k,               #sample covariance matrix
                           list(estimation))

     AllParameters <- AllParametersN
   }
   return(AllParameters) #end multistart loop after this one
  }
  stopCluster(cl)
  stopImplicitCluster()

  # Extract the likelihood values
  loglikMulti <- as.data.frame(matrix(unlist(lapply(MultistartResults2,
                                                    function(x) {x[[5]]})),ncol=1))
  row.names(loglikMulti) <- 1:(n_starts+1)

  # Obtain the number of the best starts
  best <- order(loglikMulti,decreasing = T)[1]


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Continue looping through the startset with the best loglikelihood value
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  # Extract parameter values belonging to the best loglikelihood values
  TakeOverResults <- MultistartResults2[[c(best)]]
  pi_k <- TakeOverResults[[1]]
  mu_k<- TakeOverResults[[2]]
  Lambda_k <- TakeOverResults[[3]]
  Psi_k<- TakeOverResults[[4]]
  total_logl <- TakeOverResults[[5]]
  DMV_list <- TakeOverResults[[6]]
  C_k <- TakeOverResults[[7]]
  estimation <- TakeOverResults[[8]][[1]]

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
                       mu_k,              #state-specific intercepts
                       Lambda_k,          #state-specific loading matrices
                       lapply(Psi_k,diag),#state-specific unique variances
                       total_logl,        #loglikelihood value
                       DMV_list,          #state-specific resp. probabilities
                       C_k)               #sample covariance matrix

  M_step_parameters <- list(Lambda_k,          #state-specific loading matrices
                            lapply(Psi_k,diag))#state-specific unique variances

  #*******************************************************************************#
  # NOTE: Set the loglikelihood equal to the one from the startset, the number of
  # iterations equal to n_initial_ite, and the difference between loglikelihood
  # values to 1. The values will be updated in the loop.
  #*******************************************************************************#
  LL <- total_logl
  iteration <- n_initial_ite
  differenceLL <- 1

  #*******************************************************************************#
  # NOTE: Give some value to the sum of parameter change to start the loop with
  # (will be updated in the loop).
  #*******************************************************************************#
  SumParameterChange <- 100
  cat("\n")
  cat(paste("3.Analyzing data with the best startset..."))

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
    # and intercepts mu_k
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
      mu_k[[sc]] <- SaveCov$center #Here we also obtain mu_k
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
    saveDMV<- DMV(x, Lambda_k, Psi_k, n_state, J,n_sub,mu_k)

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
                           mu_k,              #state-specific intercepts
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

  #-------------------------------------------------------------------------------#
  # Obtain standardized loadings and proportions of unique variances.
  #-------------------------------------------------------------------------------#
  x$State <- Posteriors[,1]
  #get items' standard deviation per state
  SDList <- list()
  for(i in 1:n_state){
    SDMatrix <- matrix(NA,J)
    rownames(SDMatrix) <- variable_columns
    for(item in 1:J){
      SDMatrix[item,]<- sd(unlist(x[x$State==i,variable_columns[item]]),na.rm = T)
    }
    SDList[[i]] <- SDMatrix
  }

  #get items' standard deviation across states
  SDList2 <- list()
  for(i in 1:n_state){
    SDMatrix <- matrix(NA,J)
    rownames(SDMatrix) <- variable_columns
    for(item in 1:J){
      SDMatrix[item,]<- sd(unlist(x[,variable_columns[item]]),na.rm = T)
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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                    Return Step 1 and Step 2 Results
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  cat("\n")
  #print(list(pi_k=pi_k,mu_k=mu_k,Lambda_k=Lambda_k,Psi_k=lapply(Psi_k,diag),
             #classification_posterior=Posteriors,
             #classification_errors=ModalClassificationTable,
             #classification_errors_prob=W_mod,
             #act.contraints=estimation[iteration,3]),
             #R2_entropy=R2_entropy)
  if(iteration<max_iterations){
    cat("\n")
    cat(paste("Estimation converged after",iteration,"iterations."))
  }else{
    cat("\n")
    cat("Maximum number of iterations reached without convergence.")
  }
  cat("\n")
  cat(paste("LL",round(LL,4),sep="="),"\n")
  cat("\n")



  return(list(LL=LL,
              BIC_timepoints=BIC_T,
              #BIC_cases=BIC_N,
              #number_of_cases=n_cases,
              #id_column=id_column,
              number_of_timepoints=n_sub,
              number_of_parameters=R_T,
              pi_k=pi_k,
              mu_k=mu_k,
              Lambda_k=lapply(Lambda_k,round,16),
              Lambda_k_st_w=lapply(standLambda,round,16),
              Lambda_k_st_b=lapply(standLambda2,round,16),
              Psi_k=lapply(lapply(Psi_k,diag),round,16),
              Psi_k_st_w=lapply(standPsi,round,16),
              Psi_k_st_b=lapply(standPsi2,round,16),
              explained_var=Percent_expl_var,
              standard_dev_k=SDList,
              standard_dev=SDList2,
              classification_posterior=Posteriors,
              classification_errors=ModalClassificationTable,
              classification_errors_prob=W_mod,
              iterations=iteration,
              act.contraints=estimation[iteration,3],
              R2_entropy=R2_entropy,
              seconds=requiredTime))
}
