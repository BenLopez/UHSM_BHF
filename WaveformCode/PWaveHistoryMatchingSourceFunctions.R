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
  lines( t_observation  ,  z  ,  col = 'black'  , type ='l')
  lines( t_observation  ,  z + 3*sqrt(v_me)  ,  col = 'red' , type ='l')
  lines( t_observation  ,  z - 3*sqrt(v_me)  ,  col = 'red' , type ='l')
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
  #EmulatedQS <- matrix( 0 , 100 , length(Xstar) )
  counter <- 1
  print( 'Emulating TQ Segments.' )
  for(i in 1:dim(QS_Struct$Date)[1]){
    if(QS_Struct$numvalues[i ] > as.numeric(quantile( QS_Struct$numvalues  , upper)) ){next}
    if(QS_Struct$numvalues[i ] < as.numeric(quantile( QS_Struct$numvalues  , lower)) ){next}
    tmp = 1:(min(dim(QS_Struct$Date)[2] , QS_Struct$numvalues[i ]) -1)
    EmulatorParameters$X <- (1:length(QS_Struct$Date[i,tmp]))/length(QS_Struct$Date[i,tmp])
    EmulatorParameters$Y <- QS_Struct$Value[i,tmp] - mean(QS_Struct$Value[i,tmp])
    if(length(  EmulatorParameters$X ) < 10){next}
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(Xstar , EmulatorParameters  , meanonly = 1)
    EmulatedQS[counter,] <- EmulatorOutput$E_D_fX
    counter <- counter +1
    #if(counter > 100){break}
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
PWaveHM_EmulateEstimatePAmplitude <- function(QS_Struct , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar , PriorNonImplausibleSet , Graphics = 0){
  EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct , EmulatorParameters = EmulatorParameters , Xstar )
  mQS <- apply(EmulatedQS , 2 , function(X){median(X[!is.na(X)])} ) 
  mQS <- mQS - mean( mQS )
  Implausability <- PWaveHM_CalulateImplausabilityTQSegment(Xstar , mQS , PriorNonImplausibleSet)
  XminIm <- PriorNonImplausibleSet[which.min(Implausability),]
  H = PWaveHM_CreateDesignMatrix(Xstar , XminIm , PsimulatorFunction)
  Beta = PWaveHM_CalculateBetas(H , mQS)
  E_Y = H%*%Beta
  P_Amplitude <- mQS[ which.min(abs(Xstar - XminIm[1])) ]
  
  if(Graphics == 1){
  x11()
  plot(Xstar , mQS , type ='l', xlab = 'Normalised time' , ylab = 'V')
  lines(Xstar , E_Y , type ='l' , col ='red')
  abline(P_Amplitude , 0)
  abline(v=Xstar[which.min(abs(Xstar - XminIm[1]))] )
  title(paste0('Im = ' , min(Implausability) , ' P-Amp = ' , P_Amplitude ))
  
  BC_PlotPairsFromTwoVariables(X = PriorNonImplausibleSet[which(Implausability > 1) , ] , Y = PriorNonImplausibleSet[which(Implausability < 1) , ]  , alpha = 0.1)
  
  }
  
  return(setNames(list(P_Amplitude , XminIm) ,  c('P_Amplitude' , 'XminIm')))
}
PWaveHM_HistoryMatchGroupofPwaves <- function(z , QS_Struct , PriorNonImplausibleSet , ModelDiscrepancyMatrix, Hinvstruct , Hstruct){
  # Precalculations
  
  if(!exists('Hinvstruct')){
    Hinvstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)
      return(solve(t(H)%*%H)%*%t(H))
    })}
  if(!exists('Hstruct')){
    Hstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction))
    })  
  }
  if(!exists('ModelDiscrepancyMatrix')){
    ModelDiscrepancyMatrix <- apply(PriorNonImplausibleSet , 1 , function(X){sqrt(ModelDiscrepancy(X , Xstar , PsimulatorFunction))} )
  }
  
  # History Match    
  Implausability <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    Betas = matrix(Hinvstruct[ , i] , 4  , 51)%*%t(z)
    
    Implausability[i , ] <- colMeans2( abs(t(z) - matrix(Hstruct[ , i] , 51  , 4)%*%Betas) / ModelDiscrepancyMatrix[,i])
    #apply(apply(abs(t(z) - H%*%Betas) , 2 ,  function(X){X / ModelDiscrepancyMatrix[,i]} ) , 2 , mean)
    #DP_WaitBar(i / dim(PriorNonImplausibleSet)[1])
  }
  
  # Extract Useful Data
  
  XminStruct <- PriorNonImplausibleSet[apply(Implausability , 2 , which.min) , ]
  XminStruct <- XminStruct[apply(Implausability  , 2 , min) < 2, ]

  if(length(XminStruct) < 25){
    return(rep(0,8))
  }else{
  T_start <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff(PsimulatorFunction(X , Xstar)> 0.075))  )==1 )[1]]} )
  T_end <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff(PsimulatorFunction(X , Xstar)> 0.075))  )==1 )[2]]} )
  T_end[is.na(T_end)] <- 1
  T_start[is.na(T_start)]<- mean(T_start , na.rm = T)
  PWaveDurations <- as.numeric((T_end - T_start)*(median(diff(QS_Struct$t_start) , na.rm = T) - (20*0.005) ))*1000
  PAmplitudes <- apply(XminStruct , 1 ,function(X){which.max( PsimulatorFunction(X , Xstar) ) } )
  for(i in 1:length(PAmplitudes)){
    PAmplitudes[i] <- z[i , PAmplitudes[i]]
  }
  PLocations <- apply(XminStruct , 1 ,function(X){which.max( PsimulatorFunction(X , Xstar) ) } )
  for(i in 1:length(PAmplitudes)){
    PLocations[i] <- Xstar[PLocations[i]]
  }
  
  output <- rep(0 , 8)
  output[1] <-  as.numeric(quantile(PWaveDurations , 0.98) - quantile(PWaveDurations , 0.02))
  output[2] <-  as.numeric(quantile(PWaveDurations , 0.98))
  output[3] <-  as.numeric(quantile(PWaveDurations , 0.02))
  output[4] <-  mean( PWaveDurations , na.rm = T)
  output[5] <-  var( PWaveDurations , na.rm = T) 
  output[6] <-  as.numeric((quantile(T_end , 0.99) - quantile(T_start , 0.01))*(median(diff(QS_Struct$t_start) , na.rm = T) - (20*0.005) )*1000)
  output[7] <-  mean(PAmplitudes, na.rm = T)
  output[8] <-  var(PAmplitudes , na.rm = T)
  
  return(output)
  }
}
