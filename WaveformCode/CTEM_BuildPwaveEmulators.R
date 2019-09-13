
if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath,'EmulatingCriticalThresholdDataPwave.RData') ){
  load(paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataPwave.RData'))
  Im_Crit <- PwaveOEMTrainingSet[[1]]
  X <- PwaveOEMTrainingSet[[3]]
}else{
  
  PwaveOEMTrainingSet <- CTEM_CreateTrainingsetPwaves(XPwave,
                                                      PriorNonImplausibleSet,
                                                      E_Beta,
                                                      V_Beta,
                                                      numbertrainingpoints = 1100,
                                                      numberofrepitions = 50000,
                                                      clength = 0.2 , 
                                                      q = 0.99 )
  
  Im_Crit <- PwaveOEMTrainingSet[[1]]
  X <- PwaveOEMTrainingSet[[3]]
  save( PwaveOEMTrainingSet , file = paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataPwave.RData') )
}


Training_X = X[1:1000 , ]
Training_ImCrit = Im_Crit[1:1000 , 1:2]
Validation_X = X[1001:1100 ,]
Validation_ImCrit = Im_Crit[1001:1100 , 1:2]

{ 
  PwaveEmulatorParametersMean <- BE_CreateDefaultEmulationClass()
  PwaveEmulatorParametersMean$X <- Training_X
  PwaveEmulatorParametersMean$Y <- Training_ImCrit[,1]
  PwaveEmulatorParametersMean$w <- function(X){
    return(0.001*var(PwaveEmulatorParametersMean$Y)*diag(dim(as.matrix(X))[1]))
  }
  PwaveEmulatorParametersMean$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X ) 
    return(H)}
  PwaveEmulatorParametersMean$CorrelationLength <- function(X , n){
    return(0.2*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  { 
    Xstar <- Validation_X
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = PwaveEmulatorParametersMean  )
    
    V_std <- var((Validation_ImCrit[,1] - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) ))
    R_sqM <- (1- var(Validation_ImCrit[,1] - EmulatorOutput$E_D_MX)/var(Validation_ImCrit[,1]))
    R_sq <- (1- var(Validation_ImCrit[,1] - EmulatorOutput$E_D_fX)/var(Validation_ImCrit[,1]))
    
    
    BE_PlotStdResiduals(Validation_ImCrit[,1] , EmulatorOutput) + ggtitle('Standardised Residuals Mean Implausibility Threshold ')+ ggtitle(paste0('R Squared = ' , round(R_sq , 2) , ' Variance Std Resid = ' , round(V_std , 2) ) )
  }
}
{
  PwaveEmulatorParametersMax <- BE_CreateDefaultEmulationClass()
  PwaveEmulatorParametersMax$X <- Training_X
  PwaveEmulatorParametersMax$Y <- Training_ImCrit[,2]
  PwaveEmulatorParametersMax$w <- function(X){
    return(0.01*diag(dim(as.matrix(X))[1]))
  }
  #PwaveEmulatorParametersMax$MeanFunction <- function(X){
  #  X <- as.matrix(X)
  #  H = cbind(as.matrix(1 + 0*X[,1]) , X ) 
  #  return(H)}
  PwaveEmulatorParametersMax$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X ,
              X[,1]*X[ , 2:dim(X)[2]] , 
              X[,2]*X[ , 3:dim(X)[2]] , 
              X[,4]*X[ , 5:dim(X)[2]] ,
              X^2 )
    return(H)
  }
  PwaveEmulatorParametersMax$CorrelationLength <- function(X , n){
    return(0.3*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  {Xstar <- Validation_X
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = PwaveEmulatorParametersMax  )
    
    V_std <- var((Validation_ImCrit[,2] - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) ))
    R_sqM <- (1- var(Validation_ImCrit[,2] - EmulatorOutput$E_D_MX)/var(Validation_ImCrit[,2]))
    R_sq <- (1- var(Validation_ImCrit[,2] - EmulatorOutput$E_D_fX)/var(Validation_ImCrit[,2]))
    
    
    BE_PlotStdResiduals(Validation_ImCrit[,2] , EmulatorOutput) + ggtitle('Standardised Residuals Mean Implausibility Threshold ')+ ggtitle(paste0('R Squared = ' , round(R_sq , 2) , ' Variance Std Resid = ' , round(V_std , 2) ) )
  }
}

ImThresholdMeanPwave <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSet ,EmulatorSettings = PwaveEmulatorParametersMean  )
ImThresholdMaxPwave <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSet ,EmulatorSettings = PwaveEmulatorParametersMax  )

