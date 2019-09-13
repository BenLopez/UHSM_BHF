{
  
  QSwidth = 11

  QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = QSwidth)
  if(is.null(QS_Struct$Date) ){
    StartBeat <- StartBeat + numberofBeats
    next}
  if((length(QS_Struct$Date) == 1)){
    StartBeat <- StartBeat + numberofBeats
    next}
  if((length(QS_Struct$Date) == 0)){
    StartBeat <- StartBeat + numberofBeats
    next}
  
  EmulatedQS <- FMPWaveHM_EmulateTQSegment( QS_Struct = QS_Struct , EmulatorParameters = EmulatorParameters , Xstar = seq(0.5 , 1 , 0.5/51) )
  if(is.null(dim(EmulatedQS)) ){
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
  
  EmulatorParameters3 <- EmulatorParameters
  EmulatorParameters3$X <- seq(0.5 , 1 , 0.5/51)
  EmulatorParameters3$Y <- V_me
  EmulatorParameters3$w<-function(X){
    return(0.00001*diag(length(X)))
  }
  
  SetforLookup <- seq(0.5, 1 , 0.5/10000)
  LookupTable <- BE_BayesLinearEmulatorLSEstimatesBatchMode(as.matrix(SetforLookup),EmulatorSettings = EmulatorParameters3 )$E_D_fX
  
  V_me <- t(matrix(FM_LookupPoints(Lookupinputs = SetforLookup,
                                LookupValues = LookupTable,
                                LookupPoints = matrix(XPwave ,length(XPwave) , 1)) 
                , dim(XPwave)[1] , dim(XPwave)[2]))
  V_me[V_me <0 ] <- 0
  
  #H = PWaveHM_CreateDesignMatrix(Xstar , x=PriorNonImplausibleSet[1,] , PsimulatorFunction)
  #Hinvstruct <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + 10*diag(dim(H)[1])  ) 
  
  Betas <- matrix(0 , dim(z)[2] , length(E_Beta))
  HBeta <- matrix(0, dim(FStruct)[1] , dim(FStruct)[2])
  for(i in 1:dim(Betas)[1]){
    Betas[i,] <- E_Beta + t(matrix(Hinvstruct[ , i] , 51 , 3))%*%as.matrix(( z[,i] - FStruct[,i]))
    HBeta[,i] <- matrix(Hstruct[ , i] , 51 , 3)%*%Betas[i,]
  }
  
  EZ <- HBeta + FStruct
  
  ImMatrix <- ( abs(EZ - z) / sqrt(ModelDiscrepancyMatrix + V_me) )
  #ImMatrix <- (abs(EZ - z)  ) 
  
  MeanIm <- apply(ImMatrix , 2 , function(X){mean(X , na.rm = T)})
  MaxIm <- apply(ImMatrix , 2 , function(X){max(X , na.rm = T)})
  
  MaxIm[is.infinite(MaxIm)] <- 100
  MaxIm[is.na(MaxIm)] <- 100
  MeanIm[is.infinite(MeanIm)] <- 100
  MeanIm[is.na(MeanIm)] <- 100
  

  tmplog <- ((MaxIm < (ImThresholdMaxPwave$E_D_fX + 3*sqrt(ImThresholdMaxPwave$V_D_fX))*(MeanIm < (ImThresholdMeanPwave$E_D_fX + 3*sqrt(ImThresholdMeanPwave$V_D_fX)) )) ) == 1
   
  NonImplausibleX <- PriorNonImplausibleSet[ tmplog , ]
  ImplausibleX <- PriorNonImplausibleSet[ tmplog == 0 , ]
  
  if(length(NonImplausibleX) == 5){
    NonImplausibleX = t(as.matrix(NonImplausibleX))
  }
  
  
  if(nrow(NonImplausibleX) < 2 ){
    tmplog <- ( (MaxIm < (MaxImThrehsold + 1) )*(MeanIm < (MeanImThreshold + 0.4) ) ) == 1 
    tmplog2 <- ( (MaxIm < (MaxImThrehsold + 1) )*(MeanIm < (MeanImThreshold + 0.4)) ) == 0
    
    NonImplausibleX <- PriorNonImplausibleSet[ tmplog , ]
    ImplausibleX <- PriorNonImplausibleSet[ tmplog2 , ]
  
    }
  
  if(length(NonImplausibleX) == 5){
    NonImplausibleX = t(as.matrix(NonImplausibleX))
  }
  
  if(dovalidationplots == 1 ){
    
    
  if(nrow(as.matrix(NonImplausibleX) ) != 0 ){
  BC_PlotPairsFromThreeVariables( X =ImplausibleX[which(ImplausibleX[,5]==0)[1:2000],] , Y = ImplausibleX[which(ImplausibleX[,5]!=0)[1:2000],] , Z = NonImplausibleX[1:min(2000 , dim(NonImplausibleX)[1]),]  , alpha = c(0.1 , 0.1 , 0.1), labels = c('X1' ,'X2'  , 'X3','X4','X5' ) , main = TeX('Acceptable Matches Average P-wave') )
    }
 
  x11()
  MeanIm[tmplog ==F ] = 100
  index <- which.min(MeanIm)
  p1 <- ggplot(data = data.frame( t = XPwave[  index, ] , z = EZ[  , index] ) , aes(t , z)) + 
    geom_line( col = 'blue') +
    geom_line(data = data.frame( t =  XPwave[  index, ] , z = EZ[  , index] + 3*sqrt(ModelDiscrepancyMatrix[ , index] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = EZ[  , index] - 3*sqrt(ModelDiscrepancyMatrix[ , index] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = z[  ,index]  ) , aes(t , z) , col ='black') + 
    geom_line(data = data.frame( t =  XPwave[ index, ] , z = matrix(Hstruct[,index] , 51 , 3)%*%Betas[index,]) , aes(t , z) , col = 'green')+
    ggtitle(TeX(paste0('Acceptable Match' )  ))
  print(p1)
  }

  if(sum(NonImplausibleX[,5] ==0) > 0 && sum(NonImplausibleX[,5] !=0) > 0){
    ImDiff <- sum(c( min(MaxIm[which(((PriorNonImplausibleSet[ , 5] ==0)*(tmplog))==1 )])/3 , 
                 min(MeanIm[which(((PriorNonImplausibleSet[ , 5] ==0)*(tmplog))==1 )])/0.8 ) - 
      c( min(MaxIm[which(((PriorNonImplausibleSet[ , 5] !=0)*(tmplog))==1 )])/3 , 
            min(MeanIm[which(((PriorNonImplausibleSet[ , 5] !=0)*(tmplog))==1 )])/0.8 )) 
    
  }
  
  if(sum(NonImplausibleX[,5] ==0) > 0 && sum(NonImplausibleX[,5] !=0) > 0 & dovalidationplots == 1){
  x11()
  p1 <- ggplot(data = data.frame( t = Xstar , z = EZ[  ,which(((PriorNonImplausibleSet[ , 5] !=0)*(MaxIm < MaxImThrehsold)*(MeanIm < MeanImThreshold))==1 )[1]] ) , aes(t , z)) + 
    geom_line( col = 'blue') +
    geom_line(data = data.frame( t = Xstar , z = EZ[  , which(((PriorNonImplausibleSet[ , 5] !=0)*(MaxIm < MaxImThrehsold)*(MeanIm < MeanImThreshold))==1 )[1]] + MaxImThrehsold*sqrt(ModelDiscrepancyMatrix[ ,which(((PriorNonImplausibleSet[ , 5] !=0)*(MaxIm < 3)*(MeanIm < 0.8))==1 )[1]]) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t = Xstar , z = EZ[  , which(((PriorNonImplausibleSet[ , 5] !=0)*(MaxIm < MaxImThrehsold)*(MeanIm < MeanImThreshold))==1 )[1]] - MaxImThrehsold*sqrt(ModelDiscrepancyMatrix[ , which(((PriorNonImplausibleSet[ , 5] !=0)*(MaxIm < 3)*(MeanIm < 0.8))==1 )[1]]) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t = Xstar , z = z ) , aes(t , z) , col ='black') + 
    geom_line(data = data.frame( t = Xstar , z = H%*%Betas[,which.min(MeanIm)]) , aes(t , z) , col = 'green')
  ggtitle(TeX(paste0('Acceptable Match' )  ))
  print(p1)
  }  
 
  
}

#plot( Xstar ,  z , type ='l')
#for(i in 1 : dim(NonImplausibleX)[1] ){
#  lines(Xstar , EZ[,which(( (MaxIm < MaxImThrehsold)*(MeanIm < MeanImThreshold) ) == 1)[i]] , col = rgb(0,0,1,alpha = 0.1))
#}