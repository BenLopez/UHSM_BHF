{
  
  if(!exists('QSwidth')){
    QSwidth = 12
  }
  
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
  
  EmulatedQS <- FMPWaveHM_EmulateTQSegment( QS_Struct = QS_Struct , EmulatorParameters = EmulatorParameters , Xstar = Xstar )
  if(is.null(dim(EmulatedQS)) ){
    StartBeat <- StartBeat + numberofBeats
    next}
  
  MeanImThreshold <- 0.8
  MaxImThrehsold <- 3
  
  z <- apply(EmulatedQS , 2, median)
  V_me <- apply(EmulatedQS , 2, var)/dim(EmulatedQS)[1]
  
  H = PWaveHM_CreateDesignMatrix(Xstar , x=PriorNonImplausibleSet[1,] , PsimulatorFunction)
  Hinvstruct <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + 10*diag(dim(H)[1])  ) 
  
  Betas = as.vector(E_Beta) + (Hinvstruct%*%( z - FStruct ))
  
  EZ <- H%*%Betas + FStruct
  
  ImMatrix <- (abs(EZ - z) / sqrt(ModelDiscrepancyMatrix  ) )
  
  MeanIm <- apply(ImMatrix , 2 , function(X){mean(X , na.rm = T)})
  MaxIm <- apply(ImMatrix , 2 , function(X){max(X , na.rm = T)})
  
  MaxIm[is.infinite(MaxIm)] <- 100
  MaxIm[is.na(MaxIm)] <- 100
  MeanIm[is.infinite(MeanIm)] <- 100
  MeanIm[is.na(MeanIm)] <- 100
  

  tmplog <- ( (MaxIm < MaxImThrehsold)*(MeanIm < MeanImThreshold) ) == 1
  
  NonImplausibleX <- PriorNonImplausibleSet[ tmplog , ]
  ImplausibleX <- PriorNonImplausibleSet[ ( (MaxIm < MaxImThrehsold)*(MeanIm < MeanImThreshold) ) == 0 , ]
  
  
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
  p1 <- ggplot(data = data.frame( t = Xstar , z = EZ[  , which.min(MeanIm)] ) , aes(t , z)) + 
    geom_line( col = 'blue') +
    geom_line(data = data.frame( t = Xstar , z = EZ[  , which.min(MeanIm)] + MaxImThrehsold*sqrt(ModelDiscrepancyMatrix[ , which.min(MeanIm)] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t = Xstar , z = EZ[  , which.min(MeanIm)] - MaxImThrehsold*sqrt(ModelDiscrepancyMatrix[ , which.min(MeanIm)] ) ) , aes(t , z) , col ='red') +
    geom_line(data = data.frame( t = Xstar , z = z ) , aes(t , z) , col ='black') + 
    geom_line(data = data.frame( t = Xstar , z = H%*%Betas[,which.min(MeanIm)]) , aes(t , z) , col = 'green')
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