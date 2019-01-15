AFD_ExtractIHVAFScore <- function( RWaveExtractedData , binlims = c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3) , n = 250 ){
  binmatrix <- CalulateBinMatrix(RWaveExtractedData , binlims , n )
  t <- RWaveExtractedData[!is.na(binmatrix[ , 1]) ,1] 
  binmatrix <-  binmatrix[!is.na(binmatrix[ , 1]),]
  output <- setNames(data.frame(t , 1/apply(binmatrix , 1 , var)) , c('t' , 'IHAVFScore'))
  return( output )
}
AFD_CalulateBinMatrix <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.8  , 0.05  )), n = 250){
  binmatrix <- matrix(0 , length( RWaveExtractedData$RR ) , length(binlims) - 1 )
  
  for(i in 1:(length(binlims) - 1))
  {
    binmatrix[ , i] <- smth( (RWaveExtractedData$RR > binlims[i])*(RWaveExtractedData$RR <= binlims[i+1])  , method = 'sma'  , n = n)
  }
  return( binmatrix )
}
AFD_CalulateBinMatrixKernelDensityEstimated <- function(RWaveExtractedData , binlims=  c(0, seq(from = 0  , to = 1.8  , 0.025  )) , n = 250 , window = 0.1){
  binMatrix <- AFD_CalulateBinMatrix(RWaveExtractedData , binlims , n )
  #binMatrix <- binMatrix[!is.na(binMatrix[ , 1]) , ] 
  binMatrix <- t(apply(binMatrix  , 1 , function(X){ smth(X , method ='gaussian'  , window)} ))
  #binMatrix <- binMatrix[, !is.na(binMatrix[ 1 , ])] 
  return(binMatrix)
}
AFD_ExtractIHVAFScore <- function( RWaveExtractedData , binlims = c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3) , n = 250 ){
  binmatrix <- CalulateBinMatrix(RWaveExtractedData , binlims , n )
  t <- RWaveExtractedData[!is.na(binmatrix[ , 1]) ,1] 
  binmatrix <-  binmatrix[!is.na(binmatrix[ , 1]),]
  output <- setNames(data.frame(t , 1/apply(binmatrix , 1 , function(X){var(X)})) , c('t' , 'IHAVFScore'))
  return( output )
}
AFD_ExtractNumberofModes <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.5  , 0.025  )), n = 250 , densitythresh = 0.025){
  # Function to calulate the number of modes in a local region
  
  binmatrix <- AFD_CalulateBinMatrixKernelDensityEstimated(RWaveExtractedData , binlims = binlims ,  n = n) 
  tbin <- RWaveExtractedData[!is.na(binmatrix[ , round(dim(binmatrix)[2])/2]) ,1] 
  binmatrix <-  binmatrix[!is.na(binmatrix[ , round(dim(binmatrix)[2])/2]),]
  binmatrix <- binmatrix[ , !is.na(binmatrix[1 , ]) ]
  NumberModes <- apply(binmatrix , 1 , function(X){ sum( FindLocalTurningPoints( X > densitythresh ,  1:length( binmatrix[1 , ]) ) )})
  
  return(data.frame(t = tbin, NumModes = NumberModes))
  
}
AFD_Calculatemodalmode <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.5  , 0.025  )), n = 250 , densitythresh = 0.025 , nn = 500){
  
  output <- AFD_ExtractNumberofModes(RWaveExtractedData , binlims, n  , densitythresh )
  output$NumModes[output$NumModes > 3] <- 1
  output$NumModes <- round( smth( output$NumModes , method = 'sma' , n = nn ) )
  output$NumModes[1:nn] <-output$NumModes[(nn+1)]
  output$NumModes[(length(output$NumModes) - nn) : length(output$NumModes)] <- output$NumModes[length(output$NumModes) - nn -1 ]
  output$NumModes[output$NumModes < 1] = 1;
  return(output)
}
AFD_CalulateMeanNormal <- function(AFScore , NumberModes , initialtime = 1 , priormean = 50){
  
  logicalts <- (abs(difftime( NumberModes$t , NumberModes$t[1] , units = c('hours') )) < initialtime)*(NumberModes$NumModes < 1);   
  
  if(sum(logicalts) < 1000){
    output <- priormean
  }
  if(sum(logicalts)){  
    output <- mean(AFScore$IHAVFScore[logicalts == 1] )   
  }
  return(output)
}
AFD_CreateDefaultSettings <- function(){
  SettingsAFDetection<-list()
  SettingsAFDetection[[1]] <- 100
  SettingsAFDetection[[2]] <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3)
  SettingsAFDetection[[3]] <- 100
  SettingsAFDetection[[4]] <- c(0, seq(from = 0  , to = 1.8  , 0.025  ))
  SettingsAFDetection[[5]] <- 0.02
  SettingsAFDetection[[6]] <- 250
  SettingsAFDetection[[7]] <- 1
  SettingsAFDetection[[8]] <- 10
  SettingsAFDetection[[9]] <- 5
  SettingsAFDetection[[10]] <- 60
  SettingsAFDetection[[11]] <- 2
  SettingsAFDetection[[12]] <- 200
  SettingsAFDetection[[13]] <- 50
  
  
  SettingsAFDetection <- setNames(SettingsAFDetection ,
                                  c('TimeGapThreshold',
                                    'BinlimsScore',
                                    'BandWidthScore' ,
                                    'BinlimsMM',
                                    'DensityThresholdMM',
                                    'BadnWidthMM',
                                    'InitialTimetoCalulateGlobalComponent',
                                    'TimeinAFThresh',
                                    'TimeinMMThresh',
                                    'AFScoreThresh',
                                    'ModeThresh',
                                    'AFScoreUpperThresh',
                                    'PriorBaselineScore'))  
  
  return(SettingsAFDetection)
}
AFD_DetectionWrapper <- function(RWaveExtractedData , SettingsAFDetection = AFD_CreateDefaultSettings() ){
  
  AFScore <- AFD_ExtractIHVAFScore(RWaveExtractedData ,  binlims = SettingsAFDetection[['BinlimsScore']] , n = SettingsAFDetection[['BandWidthScore']] )
  TimeGaps <- AFScore$t[c(0,abs(diff(as.numeric(AFScore$t)))) > SettingsAFDetection[['TimeGapThreshold']]]
  
  NumberModes <- AFD_Calculatemodalmode( RWaveExtractedData , 
                                         binlims = SettingsAFDetection[['BinlimsMM']] , 
                                         n = SettingsAFDetection[['BandWidthScore']] ,
                                         densitythresh = SettingsAFDetection[['DensityThresholdMM']],
                                         nn = SettingsAFDetection[['BadnWidthMM']])
  
  m <- AFD_CalulateMeanNormal(AFScore , NumberModes , initialtime = SettingsAFDetection[['InitialTimetoCalulateGlobalComponent']] , priormean = SettingsAFDetection[['PriorBaselineScore']])
  AFScore$IHAVFScore <-  AFScore$IHAVFScore - m
  
  TimeIndexofGaps <- lapply(TimeGaps , function(X){ which.min( abs(X - AFScore$t)  ) } )
  
  AFScore$IHAVFScore[as.numeric(unique(as.matrix(TimeIndexofGaps)))] <- 0
  AFScore$IHAVFScore[as.numeric(unique(as.matrix(TimeIndexofGaps))) - 1] <- 0
  
  StartEndTimesMM <- AFD_GetStartEndAF(t = AFScore$t , logicaltimeseries =  ((NumberModes$NumModes >= (SettingsAFDetection[['ModeThresh']]))*(AFScore$IHAVFScore < SettingsAFDetection[['AFScoreUpperThresh']]) ) == 1 , minutethreshold = SettingsAFDetection[['TimeinMMThresh']] )
  AFScore <- AFD_zeroMultiModal(AFScore , StartEndTimesMM)
  StartEndTimesAF <- AFD_GetStartEndAF(t = AFScore$t , logicaltimeseries =  (AFScore$IHAVFScore > SettingsAFDetection[['AFScoreThresh']]) , minutethreshold = SettingsAFDetection[['TimeinMMThresh']] )
  
  return(setNames(list(AFScore , StartEndTimesAF , StartEndTimesMM , m , NumberModes  ) , c('AFScore' , 'StartEndTimesAF' , 'StartEndTimesMM' , 'GlobalParameter' , 'NumberModes')))
}
AFD_ExtractModeStatistics <- function( f_t , binlimits = SettingsAFDetection$BinlimsMM[2:length(SettingsAFDetection$BinlimsMM)] , thresh1 = 0.025){
  
  d1 <- c(0,diff(f_t))
  d2 <- c(0,0,diff(diff(f_t)))
  peakslogical<- ((d1<thresh1)*(d2<0)*(f_t>0.01)) ==1
  peakslogical[is.na(peakslogical)]=FALSE
  peakslogical <- FindLocalTurningPoints(peakslogical ,f_t)
  
  return( data.frame(locations = binlimits[peakslogical], densities = f_t[peakslogical])) 
  
}
AFD_Checkformissingdata <- function(StartEndTimesAF , AFScore , ECGI , ECGII , ECGIII , beatlims = c(45 , 300) , readinglim = 0.65) {
  output <- StartEndTimesAF
  for(i in 1:length(StartEndTimesAF[ , 1])){
    timeinAF <- difftime( StartEndTimesAF$End[i] , StartEndTimesAF$Start[i] , units = 'secs' )
    
    minbeats <- as.numeric((beatlims[1]/60)*timeinAF)
    maxbeats <- as.numeric((beatlims[2]/60)*timeinAF)
    
    beats <- sum(((AFScore$t > StartEndTimesAF$Start[i])*(AFScore$t < StartEndTimesAF$End[i])) == 1 ) 
    
    if(beats < minbeats || beats > maxbeats){
      output <- output[-i , ]
      warning('Implausible Heart Rate')
      next
    }
    
    expectedmeasurments <- as.numeric(timeinAF)/as.numeric(abs(ECGI$Date[1] - ECGI$Date[2]))
    numberofmeasurements = c(0,0,0)
    numberofmeasurements[1] <- sum(((ECGI$Date > StartEndTimesAF$Start[i])*(ECGI$Date < StartEndTimesAF$End[i])) == 1 )
    numberofmeasurements[2] <- sum(((ECGII$Date > StartEndTimesAF$Start[i])*(ECGII$Date < StartEndTimesAF$End[i])) == 1 )
    numberofmeasurements[3] <- sum(((ECGIII$Date > StartEndTimesAF$Start[i])*(ECGIII$Date < StartEndTimesAF$End[i])) == 1 )
    
    if(sum((numberofmeasurements / expectedmeasurments) < readinglim) > 1 ){
      output <- output[-i , ]
      warning('Missing Data')
      next
    }
    
    
  }
  
  return(output)
}
AFD_CreateDistributionSummaryNames <- function(){
  return(c('Median' , 
    'Mean' , 
    'Variance' , 
    'Skewness' , 
    'Kurtosis' , 
    'Variance densities' ,
    'Max Densitiy' , 
    'Num Modes',
    'Var Modes',
    'var Mode Density',
    'IQR',
    'time' ))
}
AFD_ExtractDistributionSummaries <- function(RRStruct , n = 251 , SettingsAFDetection = AFD_CreateDefaultSettings() , MetaData){
  
if(exists('SettingsAFDetection') == FALSE){
    SettingsAFDetection <- AFD_CreateDefaultSettings()
}
if(isodd(n)){n <- n+1}

  SummaryStats <- list()
  SummaryStats[[1]] <- rollmedian( RRStruct$RR , n , na.pad = TRUE)
  SummaryStats[[2]] <- rollmean(   RRStruct$RR , n , na.pad = TRUE)
  SummaryStats[[3]] <- rollapply(  data = RRStruct$RR , width = n , FUN = var , na.pad = TRUE , align = 'center')
  SummaryStats[[4]] <- rollapply(  data = RRStruct$RR , width = n , FUN = skewness , na.pad = TRUE, align = 'center')
  SummaryStats[[5]] <- rollapply(  data = RRStruct$RR , width = n , FUN = kurtosis , na.pad = TRUE, align = 'center')
  binMatrix <- AFD_CalulateBinMatrixKernelDensityEstimated(RRStruct , n =n)
  SummaryStats[[6]] <- apply(binMatrix , 1 , function(X){var(X[X>0] , na.rm = TRUE)} )
  SummaryStats[[7]] <- apply(binMatrix , 1 , function(X){max(X , na.rm = TRUE)} )
  SummaryStats[[7]][is.infinite(SummaryStats[[7]])] <- mean(SummaryStats[[7]][(2*n):10000] , rm.na = TRUE)
  SummaryStats[[8]]  <- apply(binMatrix , 1 , function(X){ length(AFD_ExtractModeStatistics(X)$densities) } )
  SummaryStats[[8]][SummaryStats[[8]] == 0] <- 1
  SummaryStats[[8]][is.na(SummaryStats[[8]])] <- 1
  SummaryStats[[9]]  <- apply(binMatrix , 1 , function(X){ var(AFD_ExtractModeStatistics(X)$locations) } )
  SummaryStats[[9]][is.na(SummaryStats[[9]])] <- mean(SummaryStats[[9]][1:min(10000,length(SummaryStats[[9]]))] , rm.na = TRUE)
  SummaryStats[[10]] <- apply(binMatrix , 1 , function(X){ var(AFD_ExtractModeStatistics(X)$densities) } )
  SummaryStats[[10]][is.na(SummaryStats[[10]])] <- mean(SummaryStats[[10]][1:min(10000,length(SummaryStats[[10]]))] , rm.na = TRUE)
  SummaryStats[[11]] <- rollapply(RRStruct$RR , width = n, na.pad = TRUE , FUN = function(X){IQR(X, na.rm= TRUE)} )
  SummaryStats[[11]][is.na(SummaryStats[[11]])] <- mean(SummaryStats[[11]][1:min(10000,length(SummaryStats[[11]]))] , rm.na = TRUE)
  SummaryStats[[12]] <- RRStruct$t 
  
  SummaryStats <- setNames( SummaryStats ,  AFD_CreateDistributionSummaryNames())
  return(SummaryStats)
}
AFD_GetStartEndAF <- function( t , logicaltimeseries , minutethreshold = 10){
  logicaltimeseries[length(logicaltimeseries)] = FALSE
  logicaltimeseries[1] = FALSE
  d_logicaltimeseries <- diff(logicaltimeseries)
  
  if(sum(d_logicaltimeseries) == 1)
  {
    d_logicaltimeseries[length(d_logicaltimeseries)] = -1
  }
  
  if(sum(d_logicaltimeseries) == -1)
  {
    d_logicaltimeseries[1] = 1
  }
  
  output <- setNames(data.frame( t[c(d_logicaltimeseries , 0) == 1]  , t[c(d_logicaltimeseries , 0) == -1]) , c("Start" , "End"))
  output <- output[ difftime(output$End , output$Start , units = 'secs') > (minutethreshold*60) , ]
  return(output)
}
AFD_ExtractSQ<-function(ECG , RPeaks , QSwidth = 8 , index = 1){
  logicalvector <- ((ECG$Date > (RPeaks$t[index] + (QSwidth*0.005))) )*(ECG$Date <= RPeaks$t[index + 1] - (QSwidth*0.005)) == 1
  t <- ECG[logicalvector , 1]
  Values <- ECG[logicalvector , 2]
  return(setNames(list( RPeaks$t[index] , t , Values ) , c('t_start' , 'Date' , 'Value')))
}
AFD_ExtractAllSQ <- function(ECG , RPeaks , QSwidth = 8){
  
  t_start <-  RPeaks$t[1:( length(RPeaks$t) -1 )] 
  Date <- matrix(NA , length(RPeaks$t) -1 , 600)
  Value <- matrix(NA , length(RPeaks$t) -1 , 600)
  numvalues <- matrix(NA,length(t_start),1)
  ECG <- ECG[((ECG$Date <= RPeaks$t[max(which(!is.na(RPeaks$t) ))])*(ECG$Date >= RPeaks$t[1])) ==1, ]
  
  if(nrow(ECG) == 0){
    return(setNames(list(1, 1 , 1 , 1) , c('t_start' , 'Date' , 'Value' , 'numvalues')))
  }
  
  for( index in 1:(length(RPeaks$t) -1) ){
    logicalvector <- ((ECG$Date >= (RPeaks$t[index] + (QSwidth*0.005))) )*(ECG$Date <= RPeaks$t[index + 1] - (QSwidth*0.005)) == 1
    Sum_LV <- sum(logicalvector)
    if(Sum_LV == 0){next}
    if(Sum_LV > 600){next}
    Date[index ,1:Sum_LV ] <- ECG[logicalvector , 1]
    Date[index ,1:Sum_LV ] <- (Date[index ,1:Sum_LV ] - Date[index ,1 ])
    Value[index ,1:Sum_LV ] <- (ECG[logicalvector , 2])
    numvalues[index , 1] <- Sum_LV
  }  
  
  Log_1 <- apply(Date , 1 , function(X){sum(!is.na(X))}) != 0
  Log_2 <- ((apply(Date , 2 , function(X){sum(!is.na(X))})/(length(RPeaks$t) -1) )*100) > 20
  
  Date <- as.matrix(Date[ , Log_2])
  Value <- as.matrix( Value[ , Log_2])
  Date <- Date[Log_1 , ]
  Value <- Value[Log_1 , ]
  numvalues <- numvalues[Log_1 , ]
  return(setNames(list(t_start, Date , Value , numvalues) , c('t_start' , 'Date' , 'Value' , 'numvalues')))
}
AFD_zeroMultiModal <- function(AFScore , StartEndTimesMM){
  output <- AFScore
  if(nrow(StartEndTimesMM) > 0){
    for( i in 1:nrow(StartEndTimesMM) ){
      output$IHAVFScore[(output$t > StartEndTimesMM$Start[i])*(output$t < StartEndTimesMM$End[i]) == 1] <- 0
    }
  }
  return(output)
}
AFD_CalulateCDFBinMatrix <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.8  , 0.05  )), n = 250){
  binmatrix <- matrix(0 , length( RWaveExtractedData$RR ) , length(binlims) - 1 )
  
  for(i in 1:(length(binlims) - 1))
  {
    binmatrix[ , i] <- smth( (RWaveExtractedData$RR <= binlims[i])  , method = 'sma'  , n = n)
  }
  return( binmatrix )
}
AFD_CalulateCDFBinMatrixKernelDensityEstimated <- function(RWaveExtractedData , binlims=  c(0, seq(from = 0  , to = 1.8  , 0.025  )) , n = 250 , window = 0.1){
  binMatrix <- AFD_CalulateCDFBinMatrix(RWaveExtractedData , binlims , n )
  #binMatrix <- binMatrix[!is.na(binMatrix[ , 1]) , ] 
  binMatrix <- t(apply(binMatrix  , 1 , function(X){ smth(X , method ='gaussian'  , window)} ))
  rownames(binMatrix) <- RWaveExtractedData$t
  colnames(binMatrix) <- binlims[-1]
  #binMatrix <- binMatrix[, !is.na(binMatrix[ 1 , ])] 
  return(setNames(list(binMatrix , RWaveExtractedData$t) , c('CDFs' , 'time')))
}
AFD_CalulateMedianfromBinMatrix <- function(Binlims , BinMatrix){
  return(Binlims[which.min(abs(BinMatrix - 0.5))])
}
AFD_CDFCalulateBaseline <- function(BinMatix , n = 2000){
  return( as.numeric(apply(BinMatix[1:n , ] , 2 , function(X){mean(X[!is.na(X)])}) ))
}