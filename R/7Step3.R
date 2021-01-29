#' Conducts step 3 from the three-step estimation of CT-LMFA.
#'
#' 
#'
#'
#'
#'
#' @param input_file The dataset. Should contain complete cases only 
#' (unequal intervals are automatically dealt with in the CT-LMM).
#' @param interval_column The column number or name with intervals.
#' @param id_column The column number or name with subject-id.
#' @param n_state The number of states that should be estimated 
#' (has to be a scalar).
#' @param fitStep1Step2 The output file created with Step1Step2().
#' @param transitionCovariates Vector with covariates for the transition intensities.
#' @param initialCovariates Vector with covariates for the initial state probabilities.
#' @param i.method The type of optimization method that should be used.
#' @param i.maxit The maximum number of iterations that should be used.
#' @param i.reltol The tolerance to evaluate convergend that should be used.
#' @param n_q The number of start values for the transition intensity parameters that should be used.
#' @param n_initial_ite The number of initial iterations for the different start sets that should be used.
#' @param previousCov Indicates whether the covariate at t or t-1 should be used.
#
#'
#' @return Returns .
#'
#' @examples
#' \dontrun{
#' step3Results <- Step3(input_file,
#' interval_column,
#' id_column,
#' n_state,
#' fitStep1Step2,
#' transitionCovariates = NULL,
#' initialCovariates = NULL,
#' i.method = "BFGS",
#' i.maxit = 10000,
#' i.reltol = 1e-8,
#' n_q = 10,
#' n_initial_ite = 10,
#' previousCov = FALSE)
#' }
#' @export



Step3 <- function(input_file,
                  interval_column,
                  id_column,
                  n_state,
                  fitStep1Step2,
                  transitionCovariates = NULL,
                  initialCovariates = NULL,
                  i.method = "BFGS",
                  i.maxit = 10000,
                  i.reltol = 1e-8,
                  n_q = 10,
                  n_initial_ite = 10,
                  previousCov = FALSE){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                                   Step 3
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ptm <- proc.time()
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Obtain all necessary elements (from user input orr from step 1 and 2).
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Obtain the time_column from interval_column
  n_cases <- length(unlist(unique(input_file[,id_column])))
  newData <- c()
  if(!is.null(interval_column)){
    for(i in 1:n_cases){
      datai <- subset(input_file,get(
        noquote(id_column))==unlist(unique(input_file[,id_column]))[i])
      int <- c(0)
      for(ti in 2:nrow(datai)){
        int[ti] <- int[ti-1]+unlist(datai[,interval_column])[ti]
      }
      datai$time <- int
      newData<-rbind(newData, datai)
    }
  }else{
    for(i in 1:n_cases){
      datai <- subset(input_file,get(
        noquote(id_column))==unique(input_file[,id_column])[i])
      datai$time <- seq(1:nrow(datai))
      newData<-rbind(newData, datai)
    }
  }
  
  # Define the time_column.
  time_column <- "time"

  # Add the classification from step 1 and 2 to the internally-used dataset.
  newData$State <- c(fitStep1Step2$classification_posterior[,"Modal"])

  # Get the classification from step 1 and 2.
  W_mod <-  fitStep1Step2$classification_errors_prob
  
  if(sum(W_mod<1e-7)>0){
    W_mod[which(W_mod<1e-7)] <- 1e-7
  }
  
  W_mod <- W_mod/rowSums(W_mod)
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
  
  # Obtain a list of initial values (considering the average time-interval)
  # Caclulate average time-interval.
  myid.uni <- unlist(unique(newData[,id_column]))
  myid.uni_length <-length(unlist(myid.uni))
  newData2 <- c()
  for (i in 1:myid.uni_length) {
    temp<-subset(newData, get(noquote(id_column))==myid.uni[i])
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
  probabilityVector <- runif(n_initial_ite,min = 0.5,max=1)
  for(i in 1:n_initial_ite){
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
  for(i in 1:n_q){
    step3Results <-  suppressWarnings(
      msm(as.formula(
        paste("State", "~",time_column, sep="")), 
        subject = get(noquote(id_column)), 
        data = as.data.frame(newData),
        qmatrix = initialQm[[i]],
        ematrix = W_mod,
        est.initprobs=TRUE,
        hessian = TRUE,
        fixedpars = c(fixed_responseprobabilities)+additionalCounts,
        method  = i.method,
        control=list(maxit = n_initial_ite,reltol = i.reltol),
        covariates = defineCovariates,
        initcovariates = defineInitialCovariates))
                                          
    q_bestloglik[[i]] <-step3Results$Qmatrices$baseline
    bestloglik[[i]] <-step3Results$minus2loglik*-2
  } 
  

  # Consider the transition intensities from the best start set.
  Qm <- q_bestloglik[[which.max(sapply(bestloglik, max))]]
 
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Final Analysis with best startset.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  step3Results <-msm(
    as.formula(paste("State", "~",time_column, sep="")), 
    subject = get(noquote(id_column)), 
    data = as.data.frame(newData),
    qmatrix = Qm,
    ematrix = W_mod,
    est.initprobs=TRUE,
    hessian = TRUE,
    fixedpars = c(fixed_responseprobabilities)+additionalCounts,
    method  = i.method,
    control=list(maxit = i.maxit,reltol = i.reltol),
    covariates = defineCovariates,
    initcovariates = defineInitialCovariates)
                    
  
  requiredTime <- as.numeric((proc.time() - ptm)[3])

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Extracting the results
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------------------------------------#
  #            Preparation
  #--------------------------------------#
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
    c(paste("initial state",
            1:(n_state-1)),
      #cov on initial state probabilities
      paste(rep(iniName,
                each=n_state-1),
            rep(1:(n_state-1),
                length(iniName))),
      
      #intensities
      paste("transition intercepts",
            1:length(fixed_responseprobabilities)),
      #cov on intensities
      paste(rep(transitionCovariates,
                each=length(fixed_responseprobabilities)),
            rep(1:length(fixed_responseprobabilities),
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

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #                  --------------------------------------
  #                         Return Step 3 Results
  #                  --------------------------------------
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  cat("\n")
  print(list(loglikelihood=step3Results$minus2loglik/-2,
             estimates=round(parameterEstimates,4), WaldTests=round(waldMatrix,4)))



  return(list(loglikelihood=step3Results$minus2loglik/-2,
              classification_posterior=viterbi.msm(step3Results),
              estimates=round(parameterEstimates,4), 
              WaldTests=round(waldMatrix,4),
              seconds=requiredTime))
}