
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
                                                      clength = 0.03 , 
                                                      q = 0.99 )
  
  Im_Crit <- PwaveOEMTrainingSet[[1]]
  X <- PwaveOEMTrainingSet[[3]]
  save( PwaveOEMTrainingSet , file = paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataPwave.RData') )
}


{Training_X = X[1:1000 , ]
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
    return(0.12*(apply(X , 2 , max) - apply(X , 2 , min)) )
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
    return(0.001*diag(dim(as.matrix(X))[1]))
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
    return(0.15*(apply(X , 2 , max) - apply(X , 2 , min)) )
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
ImThresholdMaxPwave <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = PriorNonImplausibleSet ,EmulatorSettings = PwaveEmulatorParametersMax  )}

# Observed Implausibility Emulator 

{ObservedIm_rr <- seq(0,2,0.005)
  ObservedIm_rrmax <- seq(0,6,0.015)
  
X <- PwaveOEMTrainingSet[[3]]

Training_X = X[1:1000 , ]
Training_ImCrit = PwaveOEMTrainingSet[[2]][1:1000 , 1:2, ]
Validation_X = X[1001:1100 ,]
Validation_ImCrit = PwaveOEMTrainingSet[[2]][1001:1100 ,1:2 , ]

{ 
  PWaveEmulatorParametersCDFMean <- BE_CreateDefaultEmulationClass()
  PWaveEmulatorParametersCDFMean$X <- Training_X
  PWaveEmulatorParametersCDFMean$Y <- Training_ImCrit[ ,1 , ]
  PWaveEmulatorParametersCDFMean$w <- function(X){
    return(0.0001*diag(dim(as.matrix(X))[1]))
  }
  PWaveEmulatorParametersCDFMean$MeanFunction <- function(X){
     X <- as.matrix(X)
     H = cbind(as.matrix(1 + 0*X[,1]) , X , X^2) 
      return(H)}
 PWaveEmulatorParametersCDFMean$CorrelationLength <- function(X , n){
    return(0.12*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  {
    Xstar <- Validation_X
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimatesMO(xstar = Xstar ,EmulatorSettings =PWaveEmulatorParametersCDFMean  )
    
    StdResiduals <- matrix(0 , dim(EmulatorOutput$E_D_fX)[1] , dim(EmulatorOutput$E_D_fX)[2])
    
    RsqMean <- 1 - var( (matrix(EmulatorOutput$E_D_fX , length(Validation_ImCrit[ , 1 , ]) , 1 ) -  matrix(Validation_ImCrit[ , 1 , ] , length(Validation_ImCrit[ , 1 , ]) , 1) ) )/var(matrix(Validation_ImCrit[ , 1 , ] , length(Validation_ImCrit[ , 1 , ]) , 1))
    RsqMean2 <- 1 - var( (matrix(EmulatorOutput$E_D_MX , length(Validation_ImCrit[ ,1  , ]) , 1 ) -  matrix(Validation_ImCrit[ ,1  , ] , length(Validation_ImCrit[ ,1  , ]) , 1) ) )/var(matrix(Validation_ImCrit[ ,1  , ] , length(Validation_ImCrit[ ,1  , ]) , 1))
    
    
    for(i in 1:dim(Validation_ImCrit)[1]){
      StdResiduals[i,] <-  ( (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,1, ])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ))
    }  
    
    sampleofpoints <- sample(1:dim(StdResiduals)[1] , 100)
    
    x11()
    output <- ggplot( data.frame(Im = 1 - Validation_ImCrit[ sampleofpoints ,1, 100] , StdResids = StdResiduals[sampleofpoints , 100]) , aes(Im , StdResids) ) +
      geom_point( color = 'blue') + ggtitle('Max Implausibility Validation') + geom_hline( yintercept  = 3)+ geom_hline( yintercept  = -3) +
      ggtitle(paste0('R Squared = ', round(RsqMax , 3) ,' Variance Standardised Residuals = ', round(var(StdResiduals[sampleofpoints , 100]), 3 ) ) ) +
      ylab('Standardised Residuals') + xlab('P[Im > p]')
    
    
    print(output)
    
    #x11()
    #plot(xx,(EmulatorOutput$E_D_fX[1,] -  Validation_ImCrit[ 1 ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[1,1]*EmulatorOutput$SigmaHat) ) , type ='l'  , col = rgb(0,0,1,alpha =0.1) , ylim = c(0,5) , xlab = 'a' , ylab = 'StandardisedResiudals' )
    #for(i in 2:dim(Validation_ImCrit)[1]){
    #  lines(xx , (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ), col = rgb(0,0,1,alpha =0.1))
    #}
    
    
  }
}

{ 
  PWaveEmulatorParametersCDFMax <- BE_CreateDefaultEmulationClass()
  PWaveEmulatorParametersCDFMax$X <- Training_X[ ,]
  PWaveEmulatorParametersCDFMax$Y <- Training_ImCrit[ , 2 , ]
  PWaveEmulatorParametersCDFMax$w <- function(X){
    return(0.001*diag(dim(as.matrix(X))[1]))
  }
  PWaveEmulatorParametersCDFMax$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X , X^2 ) 
    return(H)}
  PWaveEmulatorParametersCDFMax$CorrelationLength <- function(X , n){
    return( 0.1*( apply(X , 2 , max) - apply(X , 2 , min) ) )
  }
  {
    Xstar <- Validation_X
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimatesMO(xstar = Xstar ,EmulatorSettings = PWaveEmulatorParametersCDFMax  )
    
    RsqMax <- 1 - var( (matrix(EmulatorOutput$E_D_fX , length(Validation_ImCrit[ ,2  , ]) , 1 ) -  matrix(Validation_ImCrit[ ,2  , ] , length(Validation_ImCrit[ ,2  , ]) , 1) ) )/var(matrix(Validation_ImCrit[ ,2  , ] , length(Validation_ImCrit[ ,1  , ]) , 1))
    RsqMax2 <- 1 - var( (matrix(EmulatorOutput$E_D_MX , length(Validation_ImCrit[ , 2 , ]) , 1 ) -  matrix(Validation_ImCrit[ , 2 , ] , length(Validation_ImCrit[ , 2 , ]) , 1) ) )/var(matrix(Validation_ImCrit[ ,2  , ] , length(Validation_ImCrit[ , 1 , ]) , 1))
    
    
    StdResiduals <- matrix(0 , dim(EmulatorOutput$E_D_fX)[1] , dim(EmulatorOutput$E_D_fX)[2])
    
    for(i in 1:dim(Validation_ImCrit)[1]){
      StdResiduals[i,] <-  ( (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,2, ])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ))
    }  
    
    sampleofpoints <- sample(1:dim(StdResiduals)[1] , 100)
    
    x11()
    output <- ggplot( data.frame(Im = 1 - Validation_ImCrit[ sampleofpoints ,2, 100] , StdResids = StdResiduals[sampleofpoints , 100]) , aes(Im , StdResids) ) +
      geom_point( color = 'blue') + ggtitle('Max Implausibility Validation') + geom_hline( yintercept  = 3)+ geom_hline( yintercept  = -3) +
      ggtitle(paste0('R Squared = ', round(RsqMax , 3) ,' Variance Standardised Residuals = ', round(var(StdResiduals[sampleofpoints , 100]), 3 ) ) ) +
      ylab('Standardised Residuals') + xlab('P[Im > p]')
    
    
    print(output)
    
  }
}
}
