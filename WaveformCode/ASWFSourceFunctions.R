ASWF_SegmentChange <- function(regionofinterest  , jump)
{
  
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

ASWF_Truncatetoregionwithdata <- function(regionofinterest , ECG)
{
  if( regionofinterest[1] <= 0 ){ regionofinterest <-  c(1:length(regionofinterest)) }
  if( regionofinterest[length(regionofinterest)]  >=  length(ECG[ , 1]) ){regionofinterest  <- c(length(ECG[ , 1]) - length(regionofinterest)):length(ECG[ , 1])}
return(regionofinterest)
}
  

ASWF_AlignRegionofInterests <- function(Waveform1 , Waveform2 , regionofinterest)
{
  startindex <- which.min( abs(Waveform2[ , 1] - Waveform1[ regionofinterest[1] , 1]) )
  endindex <- startindex + length(regionofinterest)
  regionofinterest2 <- startindex:endindex
  return(regionofinterest2)
}

ASWF_GetStartEndAF <- function( t , logicaltimeseries , minutethreshold = 10)
{
  logicaltimeseries[length(logicaltimeseries)] = FALSE
  
  
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