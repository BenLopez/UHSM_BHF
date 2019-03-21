{
  z <- EmulatedQS
  
  Implausability <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    Betas = matrix(Hinvstruct[ , i] , 4  , 51)%*%t(z)
    
    Implausability[i , ] <- colMeans2( abs(t(z) - matrix(Hstruct[ , i] , 51  , 4)%*%Betas) / ModelDiscrepancyMatrix[,i])
    #DP_WaitBar(i / dim(PriorNonImplausibleSet)[1])
  }
  
  XminStruct <- PriorNonImplausibleSet[apply(Implausability , 2 , which.min) , ]
  
  # Calculate Betas
  
  BetaStruct <- matrix(0,dim(XminStruct)[1] , 4)
  for(i in 1:dim(XminStruct)[1]){
    H = PWaveHM_CreateDesignMatrix(Xstar , XminStruct[i,] , PsimulatorFunction)
    alpha <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + diag(V_me))
    
    BetaStruct[i , ] <-  FMPWaveHM_BLUBetas( H , z = z[i,]  , E_Beta , V_Beta  , V_me  , alpha = 0)$E_z_Beta
  }
  
  MinImStruct <- apply(Implausability , 2 , min)
  
  # Calculate P-amplitude statistics
  E_PA <- mean( BetaStruct[,4] )
  V_PA <- var( BetaStruct[,4] )
  
  # Calculate Beta and drift Statistics 
  E_HM_Beta <- colMeans2( BetaStruct[ , 1:3 ] )
  V_HM_Beta <- cov( BetaStruct[ , 1:3 ] )
  
  # Calculate durations in ms
  T_start <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff( PsimulatorFunction( X , Xstar ) > 0.075) )  )==1 )[1]]} )
  T_end <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff( PsimulatorFunction( X , Xstar ) > 0.075) )  )==1 )[2]]} )
  
  PWaveDurations <- as.numeric((T_end - T_start)*(c(diff(QS_Struct$t_start) , median(diff(QS_Struct$t_start))   )))*1000
  E_PWaveDurations <- mean( PWaveDurations , na.rm = T )
  V_PWaveDurations <- var( PWaveDurations , na.rm = T)
  
  PRInterval <- as.numeric((1 - T_start)*(c(diff(QS_Struct$t_start) , median(diff(QS_Struct$t_start))   )))*1000
  E_PRInterval <- mean( PRInterval , na.rm = T)
  V_PRInterval <- var( PRInterval , na.rm = T)
}