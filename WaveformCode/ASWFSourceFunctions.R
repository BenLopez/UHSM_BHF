ASWF_SegmentChange <- function(regionofinterest  , jump){
  
jumpschoices <- ASWF_CreateJumpChoices()
if( jump == jumpschoices[1] ){ regionofinterest <- regionofinterest + length(regionofinterest) }
if( jump == jumpschoices[2] ){ regionofinterest <- regionofinterest + 10*length(regionofinterest) }
if( jump == jumpschoices[3] ){ regionofinterest <- regionofinterest + 100*length(regionofinterest) }
if( jump == jumpschoices[4] ){ regionofinterest <- regionofinterest + 500*length(regionofinterest) }
if( jump == jumpschoices[5] ){ regionofinterest <- regionofinterest + 1000*length(regionofinterest) }
if( jump == jumpschoices[6] ){ regionofinterest <- regionofinterest - length(regionofinterest) }
if( jump == jumpschoices[7] ){ regionofinterest <- regionofinterest - 10*length(regionofinterest) }
if( jump == jumpschoices[8] ){ regionofinterest <- regionofinterest - 100*length(regionofinterest) }
if( jump == jumpschoices[9] ){ regionofinterest <- regionofinterest - 500*length(regionofinterest) }
if( jump == jumpschoices[10] ){ regionofinterest <- regionofinterest - 1000*length(regionofinterest) }
if( jump == jumpschoices[11] ){ regionofinterest <- regionofinterest + as.numeric(winDialogString(message = 'Enter step and direction. For example, default is back 20.' , default = '-20'))*length(regionofinterest)}
return(regionofinterest)

}
ASWF_CreateJumpChoices <- function(){
  return(c('next Segment' , 'next 10' , 'next 100' , 'next 500' , 'next 1000' , 'previous Segment' , 'previous 10' , 'previous 100' , 'previous 500' ,'previous 1000'  , 'custom'))
}
ASWF_Truncatetoregionwithdata <- function(regionofinterest , ECG){
  if( regionofinterest[1] <= 0 ){ regionofinterest <-  c(1:length(regionofinterest)) }
  if( regionofinterest[length(regionofinterest)]  >=  length(ECG[ , 1]) ){regionofinterest  <- c(length(ECG[ , 1]) - length(regionofinterest)):length(ECG[ , 1])}
return(regionofinterest)
}
ASWF_AlignRegionofInterests <- function(Waveform1 , Waveform2 , regionofinterest){
  startindex <- which.min( abs(Waveform2[ , 1] - Waveform1[ regionofinterest[1] , 1]) )
  endindex <- startindex + length(regionofinterest)
  regionofinterest2 <- startindex:endindex
  return(regionofinterest2)
}
ASWF_GetStartEndAF <- function( t , logicaltimeseries , minutethreshold = 10){
  return( AFD_GetStartEndAF( t, logicaltimeseries , minutethreshold ) )
}
ASWF_Choosetimeandreturnindex <- function(WaveData , PatientRecord){
  tmp <- WaveData$Date[seq(1,length(WaveData$Date),1000)]
  interestingtimepoint <- which.min( abs( difftime(
    as.POSIXct( as.vector(as.character(round.POSIXt(tmp , units = c('hours')))) ),
    as.POSIXct(select.list(as.vector(as.character(unique(round.POSIXt(
      tmp , units = c('hours')))))
      , preselect = PatientRecord$FirstNewAF[1]
      , multiple = FALSE
      , title = 'Select Interest Time Point'
      , graphics = TRUE ) ), units ='hours') ))
  timeindex <- which.min( abs(difftime( WaveData$Date ,  tmp[interestingtimepoint[1]] , units = 'secs')) )
return(timeindex)
}
ASWF_LoadandCropECG <- function(path , subList, numberrep ,  timeindex,  Filetoprocess , ECGI){
  if(DP_checkfilesprocessed(path , subList , Filetoprocess) == 1){
    WaveData <- DP_LoadECG( path , subList , numberrep , ECGNum = 3 )
  }
  if(DP_checkfilesprocessed(path , subList , Filetoprocess) == 0){
    WaveData <- ECGI}  
  if(DP_checkfilesprocessed(path , subList , Filetoprocess) == 1 ){
    { 
      numberhoursbefore = 0
      numberhoursafter = abs( difftime( ECGI$Date[length(ECGI$Date)] , ECGI$Date[1]  , units ='hours' ))
      HoursBeforeAndAFter = data.frame(numberhoursbefore , numberhoursafter)
    }
    
    WaveData <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , c(ECGI$Date[1] ,ECGI$Date[length(ECGI$Date)] ) , HoursBeforeAndAFter))
    print('ECGIII loaded.')
  }
  
  if(DP_checkfilesprocessed(path , subList , Filetoprocess) == 0 ){
    WaveData = ECGI
    warning('No ECGIII processed.')
  }
  return(WaveData)
}