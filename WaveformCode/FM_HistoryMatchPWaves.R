{
  z <- EmulatedQS
  mQS <- apply( EmulatedQS , 2 , median )
  VQS <- apply( EmulatedQS , 2 , var )

  PwaveImThreshold <- 1
  
  Implausability <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  ImplausabilityAF <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  maxImplausability <-  matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  maxImplausabilityAF <-  matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  
  tz <- t(z)

  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    # Calculate Adjusted expectations
    
    Betas = as.vector(E_Beta) + (matrix(Hinvstruct[ , i] , 4  , 51)%*%( tz - EHstruct[,i] ))  
    
    # Old code
    # Betas = matrix(Hinvstruct[ , i] , 4  , 51)%*%( apply(t(z) , 2 , function(X){X - matrix(Hstruct[,i] , 51 , 4)%*%E_Beta}) ) + EBetaMat
    # Remove non-implausible points
     Betas[4 , Betas[4 , ] < -40] <- -40
     Betas[4 , Betas[4 , ] > 40] <-  40
     Betas[4 , ((Betas[4 , ] <0)*(Betas[4 , ] > -8) ) ==1  ] <- -8
     Betas[4 , ((Betas[4 , ] >0)*(Betas[4 , ] < 8) ) ==1  ] <-  8
    
    #Betas[4 , ((abs(Betas[4 , ])<10) + (abs(Betas[4 , ])>30)) > 0 ] <- 10000
    
    tmpmat <- abs(tz - matrix(Hstruct[ , i] , 51  , 4)%*%Betas) / sqrt(ModelDiscrepancyMatrix[,i])
    
    Implausability[i , ] <- colMeans2( tmpmat )
    maxImplausability[i , ] <- colMaxs( tmpmat  )
    
    Betas = matrix(HinvstructP[ , i] , 3  , 51)%*%tz
    tmpmat <- abs(tz - matrix(HstructP[ , i] , 51  , 3)%*%Betas) / sqrt(ModelDiscrepancyMatrix[,i])
    
    ImplausabilityAF[i , ] <- colMeans2( tmpmat )
    maxImplausabilityAF[i , ] <- colMaxs( tmpmat  )
    
    #DP_WaitBar(i / dim(PriorNonImplausibleSet)[1])
  }
  
  
  XminStruct <- PriorNonImplausibleSet[apply(Implausability , 2 , which.min) , ]
  
  # Calculate Betas
  
  BetaStruct <- matrix(0,dim(XminStruct)[1] , 4)
  for(i in 1:dim(XminStruct)[1]){
    H = PWaveHM_CreateDesignMatrix(Xstar , XminStruct[i,] , PsimulatorFunction)
    alpha <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + 10*diag(dim(H)[1]))
    
    BetaStruct[i , ] <-  FMPWaveHM_BLUBetas( H , z = z[i,]  , E_Beta , V_Beta  , 10*diag(dim(H)[1])  , alpha = alpha)$E_z_Beta
  }
  
  MinImStruct <- apply(Implausability , 2 , min)
  MinImStructAF <- apply(ImplausabilityAF , 2 , min)
  
  MinMaxImStruct <- apply(maxImplausability , 2 , min)
  MinMaxImStructAF <- apply(maxImplausabilityAF , 2 , min)
  
  BeatLogicalReg <- ((MinImStruct < PwaveImThreshold)*(MinMaxImStruct < 3)) == 1 
  BeatLogicalAF <- ((MinImStructAF < 1.6*PwaveImThreshold)*(MinMaxImStructAF < 3)) == 1 
  
  
  # Calculate P-amplitude statistics
  
  # Calculate Beta and drift Statistics 
  if(sum(BeatLogicalReg ) < 2 ){
    E_HM_Beta <- 0
    V_HM_Beta <- 0
    E_PA <- 0
    V_PA <- 0
    E_PWaveDurations <- NA
    V_PWaveDurations <- NA
    E_PRInterval <- NA
    V_PRInterval <- NA
  }else{
    E_PA <- median( BetaStruct[BeatLogicalReg, 4 ] )
    V_PA <- IQR( BetaStruct[BeatLogicalReg , 4 ] )
    
  E_HM_Beta <- colMeans2( BetaStruct[BeatLogicalReg , 1:3 ] )
  V_HM_Beta <- cov( BetaStruct[BeatLogicalReg , 1:3 ] )

  # Calculate durations in ms
  T_start <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff( PsimulatorFunction( X , Xstar ) > 0.075) )  )==1 )[1]]} )
  T_end <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff( PsimulatorFunction( X , Xstar ) > 0.075) )  )==1 )[2]]} )
  
  PWaveDurations <- as.numeric((T_end - T_start)*(c(diff(QS_Struct$t_start) , median(diff(QS_Struct$t_start))   )))*1000
  E_PWaveDurations <- median( PWaveDurations[BeatLogicalReg] , na.rm = T )
  V_PWaveDurations <- IQR( PWaveDurations[BeatLogicalReg] , na.rm = T)
  
  PRInterval <- as.numeric((1 - T_start)*(c(diff(QS_Struct$t_start) , median(diff(QS_Struct$t_start))   )))*1000
  E_PRInterval <- median( PRInterval[BeatLogicalReg] , na.rm = T)
  V_PRInterval <- IQR( PRInterval[BeatLogicalReg] , na.rm = T)
}
  }

