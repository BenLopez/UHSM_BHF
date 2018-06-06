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
  output <- setNames(data.frame(t , 1/apply(binmatrix , 1 , var)) , c('t' , 'IHAVFScore'))
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

AFD_CalulateMeanNormal <- function(AFScore , NumberModes , initialtime = 1){
  
logicalts <- (abs(difftime( NumberModes$t , NumberModes$t[1] , units = c('hours') )) < initialtime)*(NumberModes$NumModes < 2);   

if(sum(logicalts) < 1000){
output <- 50
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
  SettingsAFDetection[[8]] <- 9
  SettingsAFDetection[[9]] <- 9
  SettingsAFDetection[[10]] <- 60
  SettingsAFDetection[[11]] <- 2
  
  
  
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
                                    'ModeThresh'))  
  
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
  
  m <- AFD_CalulateMeanNormal(AFScore , NumberModes , initialtime = SettingsAFDetection[['InitialTimetoCalulateGlobalComponent']] )
  AFScore$IHAVFScore <-  AFScore$IHAVFScore - m
  
  TimeIndexofGaps <- lapply(TimeGaps , function(X){ which.min( abs(X - AFScore$t)  ) } )
  
  AFScore$IHAVFScore[as.numeric(unique(as.matrix(TimeIndexofGaps)))] <- 0
  AFScore$IHAVFScore[as.numeric(unique(as.matrix(TimeIndexofGaps))) - 1] <- 0
  
  StartEndTimesAF <- ASWF_GetStartEndAF(t = AFScore$t , logicaltimeseries = ( (AFScore$IHAVFScore > SettingsAFDetection[['AFScoreThresh']])*(NumberModes$NumModes < SettingsAFDetection[['ModeThresh']] )) == 1 , minutethreshold = SettingsAFDetection[['TimeinAFThresh']] )
  StartEndTimesMM <- ASWF_GetStartEndAF(t = AFScore$t , logicaltimeseries = ( (NumberModes$NumModes > (SettingsAFDetection[['ModeThresh']] -1) )) == 1 , minutethreshold = SettingsAFDetection[['TimeinMMThresh']] )

  return(setNames(list(AFScore , StartEndTimesAF , StartEndTimesMM , m , NumberModes  ) , c('AFScore' , 'StartEndTimesAF' , 'StartEndTimesMM' , 'GlobalParameter' , 'NumberModes')))
}