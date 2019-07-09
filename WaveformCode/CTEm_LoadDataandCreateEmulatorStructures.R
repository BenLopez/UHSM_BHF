
numberofsimulations = 1100
numberofrepeats = 5000
precomputedfolderpath <- DP_SelectPrecomputedFolder()

if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath,'EmulatingCriticalThresholdData.RData') ){
  load(paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdData.RData'))
  Im_Crit <- EmulatingCriticalThresholdData[[1]]
  X <- EmulatingCriticalThresholdData[[2]]
}else{
  outputstruct <- CTEm_CreateTrainingSet( x = xTotal,
                                          PriorNonImplausibleSetTotal = PriorNonImplausibleSetTotal,
                                          F_total = F_total,f_total = f_total,MD_Total = MD_Total  ,
                                          numberofsimulations = numberofsimulations,
                                          c = c , 
                                          l = l , 
                                          N = n ,
                                          numberofrepeats = numberofrepeats )
  EmulatingCriticalThresholdData <- list( Im_crit = outputstruct , X = PriorNonImplausibleSetTotal[1:numberofsimulations , ])
  save( EmulatingCriticalThresholdData , file = paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdData.RData') )
}

Training_X = X[1:1000 , ]
Training_ImCrit = Im_Crit[1:1000 , 1:2]
Validation_X = X[1001:1100 ,]
Validation_ImCrit = Im_Crit[1001:1100 , 1:2]

{ EmulatorParametersMean <- BE_CreateDefaultEmulationClass()
  EmulatorParametersMean$X <- Training_X[ , -3]
  EmulatorParametersMean$Y <- Training_ImCrit[,1]
  EmulatorParametersMean$w <- function(X){
    return(0.001*diag(dim(as.matrix(X))[1]))
  }
  EmulatorParametersMean$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X ) 
    return(H)}
  EmulatorParametersMean$MeanFunction <- function(X){
   X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X ,
          X[,1]*X[ , 2:dim(X)[2]] , 
         X[,2]*X[ , 3:dim(X)[2]] , 
         X[,4]*X[ , 5:dim(X)[2]] , 
        X[,6]*X[ , 7:dim(X)[2]] ,
       X[,7]*X[ , 8:dim(X)[2]] ,
      X[,8]*X[ , 9:dim(X)[2]] ,
     X^2 )
    return(H)
  }
  EmulatorParametersMean$CorrelationLength <- function(X , n){
    return(0.35*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  {Xstar <- Validation_X[ , -3]
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = EmulatorParametersMean  )
    
    V_std <- var((Validation_ImCrit[,1] - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) ))
    R_sq <- (1- var(Validation_ImCrit[,1] - EmulatorOutput$E_D_fX)/var(Validation_ImCrit[,1]))
    
    
    BE_PlotStdResiduals(Validation_ImCrit[,1] , EmulatorOutput) + ggtitle('Standardised Residuals Mean Implausibility Threshold ')+ ggtitle(paste0('R Squared = ' , round(R_sq , 2) , ' Variance Std Resid = ' , round(V_std , 2) ) )
  }
}
{EmulatorParametersMax <- BE_CreateDefaultEmulationClass()
  EmulatorParametersMax$X <- Training_X[ , -3]
  EmulatorParametersMax$Y <- Training_ImCrit[,2]
  EmulatorParametersMax$w <- function(X){
    return(0.001*diag(dim(as.matrix(X))[1]))
  }
  #EmulatorParametersMax$MeanFunction <- function(X){
  #  X <- as.matrix(X)
  #  H = cbind(as.matrix(1 + 0*X[,1]) , X ) 
  #  return(H)}
  EmulatorParametersMax$MeanFunction <- function(X){
   X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X ,
          X[,1]*X[ , 2:dim(X)[2]] , 
         X[,2]*X[ , 3:dim(X)[2]] , 
         X[,4]*X[ , 5:dim(X)[2]] , 
        X[,6]*X[ , 7:dim(X)[2]] ,
       X[,7]*X[ , 8:dim(X)[2]] ,
     X[,8]*X[ , 9:dim(X)[2]] ,
     X^2 )
    return(H)
  }
  EmulatorParametersMax$CorrelationLength <- function(X , n){
    return(0.275*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  {Xstar <- Validation_X[ , -3]
   EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = EmulatorParametersMax  )
  
    V_std <- var((Validation_ImCrit[,2] - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) ))
    R_sq <- (1- var(Validation_ImCrit[,2] - EmulatorOutput$E_D_fX)/var(Validation_ImCrit[,2]))
  
  
    BE_PlotStdResiduals(Validation_ImCrit[,2] , EmulatorOutput) + ggtitle('Standardised Residuals Mean Implausibility Threshold ')+ ggtitle(paste0('R Squared = ' , round(R_sq , 2) , ' Variance Std Resid = ' , round(V_std , 2) ) )
  }
}



ImThresholdMeanRegular <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSetRegular[ , -3] ,EmulatorSettings = EmulatorParametersMean  )
ImThresholdMaxRegular <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSetRegular[ , -3] ,EmulatorSettings = EmulatorParametersMax  )
ImThresholdMeanIrregularlyIrregular <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSetRegularyIreRegular[ , -3] ,EmulatorSettings = EmulatorParametersMean  )
ImThresholdMaxIrregularlyIrregular <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSetRegularyIreRegular[ , -3] ,EmulatorSettings = EmulatorParametersMax  )
ImThresholdMeanTotal <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSetTotal[ , -3] ,EmulatorSettings = EmulatorParametersMean  )
ImThresholdMaxTotal <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar =  PriorNonImplausibleSetTotal[ , -3] ,EmulatorSettings = EmulatorParametersMax  )


ImThresholdMeanRegular <- ImThresholdMeanRegular$E_D_fX + 3*sqrt(ImThresholdMeanRegular$V_D_fX)
ImThresholdMaxRegular <- ImThresholdMaxRegular$E_D_fX + 3*sqrt(ImThresholdMaxRegular$V_D_fX)
ImThresholdMeanIrregularlyIrregular <- ImThresholdMeanIrregularlyIrregular$E_D_fX + 3*sqrt(ImThresholdMeanIrregularlyIrregular$V_D_fX)
ImThresholdMaxIrregularlyIrregular <- ImThresholdMaxIrregularlyIrregular$E_D_fX + 3*sqrt(ImThresholdMaxIrregularlyIrregular$V_D_fX)
ImThresholdMeanTotal <- ImThresholdMeanTotal$E_D_fX + 3*sqrt(ImThresholdMeanTotal$V_D_fX)
ImThresholdMaxTotal <- ImThresholdMaxTotal$E_D_fX + 3*sqrt(ImThresholdMaxTotal$V_D_fX)
