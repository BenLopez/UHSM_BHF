{
  
  QSwidth = 12

  QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
  if(is.null(QS_Struct$Date) ){
    QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGI , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
    #StartBeat <- StartBeat + numberofBeats
    #next
    }
  if((length(QS_Struct$Date) == 1)){
    QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGI , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
    # StartBeat <- StartBeat + numberofBeats
    #next
      }
  if((length(QS_Struct$Date) == 0)){
    #QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGI , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
    StartBeat <- StartBeat + numberofBeats
    next
    }
  
  if(is.null(QS_Struct$Date) ){
    #QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGI , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
    StartBeat <- StartBeat + numberofBeats
    next
  }
  if((length(QS_Struct$Date) == 1)){
    #QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGI , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
     StartBeat <- StartBeat + numberofBeats
    next
  }
  if((length(QS_Struct$Date) == 0)){
    #QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGI , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
    StartBeat <- StartBeat + numberofBeats
    next
  }
  
  EmulatedQS <- FMPWaveHM_EmulateTQSegment( QS_Struct = QS_Struct , EmulatorParameters = EmulatorParameters , Xstar = seq(0.5 , 1 , 0.5/51) )
  if(is.null(dim(EmulatedQS))){
    ECGAbscence[i] <- 1
    StartBeat <- StartBeat + numberofBeats
    next}
  if(nrow(EmulatedQS) < 100 ){
    ECGAbscence[i] <- 1
    StartBeat <- StartBeat + numberofBeats
    next}
  
 
  z <- apply(EmulatedQS , 2, median)
  V_me <- apply(EmulatedQS , 2, var)/dim(EmulatedQS)[1]
  
  EmulatorParameters2 <- EmulatorParameters
  EmulatorParameters2$X <- seq(0.5 , 1 , 0.5/51)
  EmulatorParameters2$Y <- z
  EmulatorParameters2$w<-function(X){
    return(0.01*diag(length(X)))
  }
  
  SetforLookup <- seq(0.5, 1 , 0.5/10000)
  LookupTable <- BE_BayesLinearEmulatorLSEstimatesBatchMode(as.matrix(SetforLookup),EmulatorSettings = EmulatorParameters2 )$E_D_fX
  
  z <- t(matrix(FM_LookupPoints(Lookupinputs = SetforLookup,
                             LookupValues = LookupTable,
                             LookupPoints = matrix(XPwave ,length(XPwave) , 1)) 
             , dim(XPwave)[1] , dim(XPwave)[2]))
  
  #EmulatorParameters3 <- EmulatorParameters
  #EmulatorParameters3$X <- seq(0.5 , 1 , 0.5/51)
  #EmulatorParameters3$Y <- V_me
  #EmulatorParameters3$w<-function(X){
  #  return(0.00001*diag(length(X)))
  #}
  
  #SetforLookup <- seq(0.5, 1 , 0.5/10000)
  #LookupTable <- BE_BayesLinearEmulatorLSEstimatesBatchMode(as.matrix(SetforLookup),EmulatorSettings = EmulatorParameters3 )$E_D_fX
  
  #V_me <- t(matrix(FM_LookupPoints(Lookupinputs = SetforLookup,
  #                              LookupValues = LookupTable,
  #                              LookupPoints = matrix(XPwave ,length(XPwave) , 1)) 
  #              , dim(XPwave)[1] , dim(XPwave)[2]))
  #V_me[V_me <0 ] <- 0
  
  #H = PWaveHM_CreateDesignMatrix(Xstar , x=PriorNonImplausibleSet[1,] , PsimulatorFunction)
  #Hinvstruct <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + 10*diag(dim(H)[1])  ) 
  
  if( !exists('HHVBetaStruct') || !exists('VHBeta')){
  
  VHBeta <- matrix(0, dim(FStruct)[1] , dim(FStruct)[2])
  for(iii in 1:dim(Betas)[1]){
    HHInv <- matrix(Hinvstruct[ , iii] , 3 , 51)
    HH <- matrix(Hstruct[ , iii] , 51 , 3)
    VHBeta[,iii] <- diag( HH%*%( V_Beta - HHInv%*%(HH%*%V_Beta))%*%t(HH) )
  }
  }
  
  Betas <- matrix(0 , dim(z)[2] , length(E_Beta))
  HBeta <- matrix(0, dim(FStruct)[1] , dim(FStruct)[2])
  
  for(iii in 1:dim(Betas)[1]){
    # Might be able to vectorise this
    Betas[iii,] <- E_Beta +  matrix(Hinvstruct[ , iii] , 3 , 51)%*%( z[,iii] - FStruct[,iii])
    HBeta[,iii] <-  matrix(Hstruct[ , iii] , 51 , 3)%*%Betas[iii,]
    #Sigma <- var(z[,iii] - (FStruct[,iii]+ (HBeta[,iii])) )
  }
  
  {
  EZ <- HBeta + FStruct
  
  ImMatrix <- ( abs(EZ - z) / sqrt(ModelDiscrepancyMatrix + VHBeta + 0.1) )
  ImMatrix[is.na(ImMatrix)] <- 100  
  
  MeanIm <- colMeans(ImMatrix)
  MaxIm <- colMaxs(ImMatrix)
  
  MaxIm[is.infinite(MaxIm)] <- 100
  MaxIm[is.na(MaxIm)] <- 100
  MeanIm[is.infinite(MeanIm)] <- 100
  MeanIm[is.na(MeanIm)] <- 100
  }
  

  tmplog <- ((MaxIm < (ImThresholdMaxPwave$E_D_fX + 3*sqrt(ImThresholdMaxPwave$V_D_fX))*(MeanIm < (ImThresholdMeanPwave$E_D_fX + 3*sqrt(ImThresholdMeanPwave$V_D_fX)) )) ) == 1
   
  NonImplausibleX <- PriorNonImplausibleSet[ tmplog , ]
  ImplausibleX <- PriorNonImplausibleSet[ tmplog == 0 , ]
  
  if(length(NonImplausibleX) == 5){
    NonImplausibleX = t(as.matrix(NonImplausibleX))
  }
  
  
  if(length(NonImplausibleX) < 12 ){
    tmplog <- ( (MaxIm < (1.2*ImThresholdMaxPwave$E_D_fX ) )*(MeanIm < (1.2*ImThresholdMeanPwave$E_D_fX) ) ) == 1 
    tmplog2 <- ( (MaxIm < (1.2*ImThresholdMaxPwave$E_D_fX ) )*(MeanIm < (1.2*ImThresholdMeanPwave$E_D_fX)) ) == 0
    
    NonImplausibleX <- PriorNonImplausibleSet[ tmplog , ]
    ImplausibleX <- PriorNonImplausibleSet[ tmplog2 , ]
  
    }
  
  if(length(NonImplausibleX) == 6){
    NonImplausibleX = t(as.matrix(NonImplausibleX))
  }
  
  Wave2Pwave <- FALSE
  if(sum(tmplog) > 2 & tmplog[length(tmplog)] ==1  ){
  Wave2Pwave <-   FM_PwaveWave2Logic(NonImplausibleX,PriorNonImplausibleSet,PWaveEmulatorParametersCDFMean,PWaveEmulatorParametersCDFMax,ObservedIm_rr,ObservedIm_rrmax)
  if(Wave2Pwave == TRUE ){
    disp('Pwave abscent in wave two.')
  }else{
    disp('Pwave present in wave two.')
  }
  
  }
  
if(dovalidationplots == 1 ){
    
    
if(nrow(as.matrix(NonImplausibleX) ) != 0 ){
  BC_PlotPairsFromThreeVariables( X =ImplausibleX[which(ImplausibleX[,5]==0)[1:2000],] , Y = ImplausibleX[which(ImplausibleX[,5]!=0)[1:2000],] , Z = NonImplausibleX[1:min(2000 , dim(NonImplausibleX)[1]),]  , alpha = c(0.1 , 0.1 , 0.1), labels = c('X1' ,'X2'  , 'X3','X4','X5' ) , main = TeX('Acceptable Matches Average P-wave') )
}
 
  x11()
  MeanIm[tmplog ==F ] = 100
  index <- which.min(MeanIm[1:(length(MeanIm)-2)])
  #index <- dim(PriorNonImplausibleSet)[1]
  p1 <- ggplot(data = data.frame( t = XPwave[  index, ] , z = EZ[  , index] ) , aes(t , z)) + 
    geom_line( col = 'blue') +
    geom_line(data = data.frame( t =  XPwave[  index, ] , z = EZ[  , index] + 3*sqrt(ModelDiscrepancyMatrix[ , index] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = EZ[  , index] - 3*sqrt(ModelDiscrepancyMatrix[ , index] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = z[  ,index]  ) , aes(t , z) , col ='black') + 
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = matrix(Hstruct[,index] , 51 , 3)%*%Betas[index,]) , aes(t , z) , col = 'green')+
    ggtitle(TeX(paste0('Acceptable Match' )  ))
  print(p1)
  x11()
  index <- length(MeanIm)
  #index <- dim(PriorNonImplausibleSet)[1]
  p1 <- ggplot(data = data.frame( t = XPwave[  index, ] , z = EZ[  , index] ) , aes(t , z)) + 
    geom_line( col = 'blue') +
    geom_line(data = data.frame( t =  XPwave[  index, ] , z = EZ[  , index] + 3*sqrt(ModelDiscrepancyMatrix[ , index] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = EZ[  , index] - 3*sqrt(ModelDiscrepancyMatrix[ , index] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = z[  ,index]  ) , aes(t , z) , col ='black') + 
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = matrix(Hstruct[,index] , 51 , 3)%*%Betas[index,]) , aes(t , z) , col = 'green')+
    ggtitle(TeX(paste0('Acceptable Match' )  ))
  print(p1)
  
}
 
}



