# Set up design matrix
PWaveHM_CreateDesignMatrix <- function(t_observation , x , PsimulatorFunction){
  t_observation <- as.matrix(t_observation)
  return( cbind( matrix(1 , dim(t_observation)[1] , 1 ) , t_observation , as.matrix(PsimulatorFunction(x , t_observation))  ) )
}
PWaveHM_CalculateBetas <- function(H , z){
  z <- as.matrix(z)
  return(  solve(t(H)%*%H)%*%t(H)%*%z  )
}
PWaveHM_CalculateMSE <- function(H , Beta , z){
  return( mean(abs(z - H%*%Beta)) )
}
PWaveHM_PlotMatch <- function(t_observation , E_Z , z , v_me){
  x11()
  plot(  t_observation  ,  E_Z , col = 'blue'  ,  xlab = 'Normalised Time'  ,  ylab = 'Observed' , type = 'l'  , ylim = c(-40,40))
  points( t_observation  ,  z  ,  col = 'black' )
  points( t_observation  ,  z + 2*sqrt(v_me)  ,  col = 'red' )
  points( t_observation  ,  z - 2*sqrt(v_me)  ,  col = 'red' )
}
PWaveHM_CalculateImplausability <- function( t_observation , x ,  z ){
  H = PWaveHM_CreateDesignMatrix(t_observation , x , PsimulatorFunction)
  Beta = PWaveHM_CalculateBetas(H , z)
  MSE = PWaveHM_CalculateMSE(H , Beta , z)
  return( MSE )
} 
PWaveHM_CreateDefaultEmulationclass <- function(){ 
  EmulatorParameters <- BE_CreateDefaultEmulationClass()
  EmulatorParameters$w <- function(X){
    return( 0.5*diag(length(X) ) )
  }
  EmulatorParameters$CorrelationLengthfunction<- function(X , n){
    return(0.2)
  }
  return(EmulatorParameters)
}
PWaveHM_EmulateTQSegment <- function( QS_Struct , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar  , upper = 0.9 , lower = 0.1){
  
  tmp = 1:(min(dim(QS_Struct$Date)[2] , QS_Struct$numvalues[1 ]) -1)
  EmulatorParameters$X <- (1:length(QS_Struct$Date[1,tmp]))/length(QS_Struct$Date[1,tmp])
  EmulatorParameters$Y <- QS_Struct$Value[1,tmp] - mean(QS_Struct$Value[1,tmp])
  
  EmulatedQS <- matrix( 0 , dim(QS_Struct$Date)[1] , length(Xstar) )
  counter <- 1
  print( 'Emulating TQ Segments.' )
  for(i in 1:dim(QS_Struct$Date)[1]){
    if(QS_Struct$numvalues[i ] > as.numeric(quantile( QS_Struct$numvalues  , upper)) ){next}
    if(QS_Struct$numvalues[i ] < as.numeric(quantile( QS_Struct$numvalues  , lower)) ){next}
    tmp = 1:(min(dim(QS_Struct$Date)[2] , QS_Struct$numvalues[i ]) -1)
    EmulatorParameters$X <- (1:length(QS_Struct$Date[i,tmp]))/length(QS_Struct$Date[i,tmp])
    EmulatorParameters$Y <- QS_Struct$Value[i,tmp] - mean(QS_Struct$Value[i,tmp])
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(Xstar , EmulatorParameters)
    EmulatedQS[counter,] <- EmulatorOutput$E_D_fX
    counter <- counter +1
  }  
  print('TQ Segments Emulated.' )
  EmulatedQS <- EmulatedQS[-(counter:dim(EmulatedQS)[1] ), ]
  return(EmulatedQS)
}
PsimulatorFunction <- function( x , t_observation ){
  Pcen <- x[1]
  Pwidth <- x[2]
  return( ECGSim_Gaussian( x = t_observation , mu = Pcen , sigma = Pwidth ) )
}
PWaveHM_CalulateImplausabilityTQSegment <- function(Xstar , z ,   PriorNonImplausibleSet , ImplausabilityFunction = PWaveHM_CalculateImplausability){
  Implausability <- matrix(0 , dim(PriorNonImplausibleSet)[1] , 1)
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    Implausability[i,] = ImplausabilityFunction(Xstar , PriorNonImplausibleSet[i,] , z)
  }
  return(Implausability)
}
PWaveHM_CalulateBetasTQSegment <- function(Xstar , z ,   PriorNonImplausibleSet){
  H <- PWaveHM_CreateDesignMatrix(Xstar ,PriorNonImplausibleSet[1,] , PsimulatorFunction )
  Betas <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(H)[2])
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    Betas[i,] = PWaveHM_CalculateBetas(PWaveHM_CreateDesignMatrix(Xstar ,PriorNonImplausibleSet[i,] , PsimulatorFunction ), z)
  }
  return(Betas)
}
PWaveHM_EmulateEstimatePAmplitude <- function(QS_Struct , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar , PriorNonImplausibleSet){
  EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct , EmulatorParameters = EmulatorParameters , Xstar )
  mQS <- apply(EmulatedQS , 2 , median) 
  mQS <- mQS - mean( mQS )
  Implausability <- PWaveHM_CalulateImplausabilityTQSegment(Xstar , mQS , PriorNonImplausibleSet)
  XminIm <- PriorNonImplausibleSet[which.min(Implausability),]
  H = PWaveHM_CreateDesignMatrix(Xstar , XminIm , PsimulatorFunction)
  Beta = PWaveHM_CalculateBetas(H , mQS)
  E_Y = H%*%Beta
  P_Amplitude <- E_Y[ which.min(abs(Xstar - XminIm[1])) ]
  
  return(P_Amplitude)
}