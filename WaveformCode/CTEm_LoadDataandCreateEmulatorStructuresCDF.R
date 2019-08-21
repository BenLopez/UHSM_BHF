
numberofsimulations = 1100
numberofrepeats = 50000

if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath,'EmulatingCriticalThresholdDataCDF.RData') ){
  load(paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataCDF.RData'))
  Im_Crit <- EmulatingCriticalThresholdData[[1]]
  X <- EmulatingCriticalThresholdData[[2]]
}else{
  ObservedIm_xx <- seq(0,5,0.01)
  outputstruct <- CTEm_CreateTrainingSetCDF(ObservedIm_xx,xTotal,PriorNonImplausibleSetTotal,F_total,f_total,MD_Total  ,numberofsimulations = numberofsimulations, c = c , l = l , N = n ,numberofrepeats = numberofrepeats)
  EmulatingCriticalThresholdData <- list( Im_crit = outputstruct , X = PriorNonImplausibleSetTotal[1:numberofsimulations , ])
  save( EmulatingCriticalThresholdData , file = paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataCDF.RData') )
  Im_Crit <- EmulatingCriticalThresholdData[[1]]
  X <- EmulatingCriticalThresholdData[[2]]
}

ObservedIm_xx <- seq(0,5,0.01)

Training_X = X[1:1000 , ]
Training_ImCrit = Im_Crit[1:1000 , , 1:2]
Validation_X = X[1001:1100 ,]
Validation_ImCrit = Im_Crit[1001:1100 , , 1:2]

{ 
  EmulatorParametersCDFMean <- BE_CreateDefaultEmulationClass()
  EmulatorParametersCDFMean$X <- Training_X[ , -3]
  EmulatorParametersCDFMean$Y <- Training_ImCrit[ , , 1]
  EmulatorParametersCDFMean$w <- function(X){
    return(0.0001*diag(dim(as.matrix(X))[1]))
  }
  # EmulatorParametersCDFMean$MeanFunction <- function(X){
  #   X <- as.matrix(X)
  #   H = cbind(as.matrix(1 + 0*X[,1]) , X ) 
  #    return(H)}
  EmulatorParametersCDFMean$MeanFunction <- function(X){
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
  EmulatorParametersCDFMean$CorrelationLength <- function(X , n){
    return(0.25*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  {
    Xstar <- Validation_X[ , -3]
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimatesMO(xstar = Xstar ,EmulatorSettings = EmulatorParametersCDFMean  )
    
    StdResiduals <- matrix(0 , dim(EmulatorOutput$E_D_fX)[1] , dim(EmulatorOutput$E_D_fX)[2])
    
    RsqMean <- 1 - var( (matrix(EmulatorOutput$E_D_fX , length(Validation_ImCrit[ ,  , 1]) , 1 ) -  matrix(Validation_ImCrit[ ,  , 1] , length(Validation_ImCrit[ ,  , 1]) , 1) ) )/var(matrix(Validation_ImCrit[ ,  , 1] , length(Validation_ImCrit[ ,  , 1]) , 1))
    RsqMean2 <- 1 - var( (matrix(EmulatorOutput$E_D_MX , length(Validation_ImCrit[ ,  , 1]) , 1 ) -  matrix(Validation_ImCrit[ ,  , 1] , length(Validation_ImCrit[ ,  , 1]) , 1) ) )/var(matrix(Validation_ImCrit[ ,  , 1] , length(Validation_ImCrit[ ,  , 1]) , 1))
    
    
    for(i in 1:dim(Validation_ImCrit)[1]){
      StdResiduals[i,] <-  ( (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ))
    }  
    
    sampleofpoints <- sample(1:dim(StdResiduals)[1] , 100)
    
    x11()
    output <- ggplot( data.frame(x_2 = Validation_X[sampleofpoints , 8] , StdResids = StdResiduals[sampleofpoints , 116]) , aes(x_2 , StdResids) ) +
      geom_point( color = 'blue') + ggtitle('Mean Implausibility Validation') + geom_hline( yintercept  = 3)+ geom_hline( yintercept  = -3) +
    ggtitle(paste0('R Squared = ', round(RsqMean , 3) ,' Variance Std Residuals = ', round(var(StdResiduals[sampleofpoints , 116]), 3 ) ) ) +
      ylab('P[Im > p]') + xlab('p')
    
    print(output)
    
    #x11()
    #plot(xx,(EmulatorOutput$E_D_fX[1,] -  Validation_ImCrit[ 1 ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[1,1]*EmulatorOutput$SigmaHat) ) , type ='l'  , col = rgb(0,0,1,alpha =0.1) , ylim = c(0,5) , xlab = 'a' , ylab = 'StandardisedResiudals' )
    #for(i in 2:dim(Validation_ImCrit)[1]){
    #  lines(xx , (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ), col = rgb(0,0,1,alpha =0.1))
    #}
    

  }
}

{ 
  EmulatorParametersCDFMax <- BE_CreateDefaultEmulationClass()
  EmulatorParametersCDFMax$X <- Training_X[ , -3]
  EmulatorParametersCDFMax$Y <- Training_ImCrit[ , , 2]
  EmulatorParametersCDFMax$w <- function(X){
    return(0.001*diag(dim(as.matrix(X))[1]))
  }
  EmulatorParametersCDFMax$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X , X^2 ) 
    return(H)}
  EmulatorParametersCDFMax$CorrelationLength <- function(X , n){
    return( 0.35*( apply(X , 2 , max) - apply(X , 2 , min) ) )
  }
  {
    Xstar <- Validation_X[ , -3]
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimatesMO(xstar = Xstar ,EmulatorSettings = EmulatorParametersCDFMax  )
    
    
    
    RsqMax <- 1 - var( (matrix(EmulatorOutput$E_D_fX , length(Validation_ImCrit[ ,  , 2]) , 1 ) -  matrix(Validation_ImCrit[ ,  , 2] , length(Validation_ImCrit[ ,  , 2]) , 1) ) )/var(matrix(Validation_ImCrit[ ,  , 2] , length(Validation_ImCrit[ ,  , 1]) , 1))
    RsqMax2 <- 1 - var( (matrix(EmulatorOutput$E_D_MX , length(Validation_ImCrit[ ,  , 2]) , 1 ) -  matrix(Validation_ImCrit[ ,  , 2] , length(Validation_ImCrit[ ,  , 2]) , 1) ) )/var(matrix(Validation_ImCrit[ ,  , 2] , length(Validation_ImCrit[ ,  , 1]) , 1))
    
    
    StdResiduals <- matrix(0 , dim(EmulatorOutput$E_D_fX)[1] , dim(EmulatorOutput$E_D_fX)[2])
    
    for(i in 1:dim(Validation_ImCrit)[1]){
      StdResiduals[i,] <-  ( (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,, 2])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ))
    }  
    
    sampleofpoints <- sample(1:dim(StdResiduals)[1] , 100)
    
    x11()
    output <- ggplot( data.frame(Im = 1 - Validation_ImCrit[ sampleofpoints ,400, 2] , StdResids = StdResiduals[sampleofpoints , 400]) , aes(Im , StdResids) ) +
      geom_point( color = 'blue') + ggtitle('Max Implausibility Validation') + geom_hline( yintercept  = 3)+ geom_hline( yintercept  = -3) +
      ggtitle(paste0('R Squared = ', round(RsqMax , 3) ,' Variance Standardised Residuals = ', round(var(StdResiduals[sampleofpoints , 400]), 3 ) ) ) +
      ylab('P[Im > p]') + xlab('p')
      
    
    print(output)
    
  }
}

