{
  z <- EmulatedQS
  mQS <- apply( EmulatedQS , 2 , mean )
  VQS <- apply( EmulatedQS , 2 , var )
  
  {E_Beta = matrix(c(0 , 0 , 0 , 14) , 4 , 1)
    V_Beta = diag(diag(matrix(0 , 4 , 4)))
    V_Beta[1 , 1] <- 1000
    V_Beta[2 , 2] <- 1000
    V_Beta[3 , 3] <- 1000
    V_Beta[4 , 4] <- 450
    V_me = VQS
    
    C = diag(4)
    C[1 , 2] <- -0.95
    C[2 , 1] <- C[1 , 2]
    C[1 , 3] <- 0.95
    C[3 , 1] <- C[1 , 3]
    C[1 , 4] <- -0.8
    C[4 , 1] <- C[1 , 4]
    C[3 , 2] <- -0.7
    C[2 , 3] <- C[3 , 2] 
    C[4 , 2] <- 0.8
    C[2 , 4] <- C[4 , 2]
    C[4 , 3] <- -0.8
    C[3 , 4] <- C[4 , 3]
    
    V_Beta <- ((diag(sqrt(V_Beta) ))%*%t(diag(sqrt(V_Beta) )) )*C}
  
  
  Implausability <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    Betas = matrix(Hinvstruct[ , i] , 4  , 51)%*%t(z)
    
    Betas[4 , Betas[4 , ] < -30] <- -30
    Betas[4 , Betas[4 , ] > 30] <-  30
    
    Implausability[i , ] <- colMeans2( abs(t(z) - matrix(Hstruct[ , i] , 51  , 4)%*%Betas) / (2*ModelDiscrepancyMatrix[,i]) )
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
  E_PA <- mean( BetaStruct[apply(Implausability , 2 , min) < 0.6,4] )
  V_PA <- var( BetaStruct[apply(Implausability , 2 , min) < 0.6,4] )
  
  # Calculate Beta and drift Statistics 
  E_HM_Beta <- colMeans2( BetaStruct[apply(Implausability , 2 , min) < 0.6 , 1:3 ] )
  V_HM_Beta <- cov( BetaStruct[apply(Implausability , 2 , min) < 0.6 , 1:3 ] )
  
  # Calculate durations in ms
  T_start <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff( PsimulatorFunction( X , Xstar ) > 0.075) )  )==1 )[1]]} )
  T_end <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff( PsimulatorFunction( X , Xstar ) > 0.075) )  )==1 )[2]]} )
  
  PWaveDurations <- as.numeric((T_end - T_start)*(c(diff(QS_Struct$t_start) , median(diff(QS_Struct$t_start))   )))*1000
  E_PWaveDurations <- mean( PWaveDurations[apply(Implausability , 2 , min) < 0.6] , na.rm = T )
  V_PWaveDurations <- var( PWaveDurations[apply(Implausability , 2 , min) < 0.6] , na.rm = T)
  
  PRInterval <- as.numeric((1 - T_start)*(c(diff(QS_Struct$t_start) , median(diff(QS_Struct$t_start))   )))*1000
  E_PRInterval <- mean( PRInterval[apply(Implausability , 2 , min) < 0.6] , na.rm = T)
  V_PRInterval <- var( PRInterval[apply(Implausability , 2 , min) < 0.6] , na.rm = T)
}

