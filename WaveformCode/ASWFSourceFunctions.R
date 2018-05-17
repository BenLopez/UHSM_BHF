ASWF_SegmentChange <- function(regionofinterest  , jump)
{
  
jumpschoices <- c('next Segment' , 'next 10' , 'next 100' , 'next 1000' , 'previous Segment' , 'previous 10' , 'previous 100' , 'previous 1000'  )
if( jump == jumpschoices[1] ){ regionofinterest <- regionofinterest + length(regionofinterest) }
if( jump == jumpschoices[2] ){ regionofinterest <- regionofinterest + 10*length(regionofinterest) }
if( jump == jumpschoices[3] ){ regionofinterest <- regionofinterest + 100*length(regionofinterest) }
if( jump == jumpschoices[4] ){ regionofinterest <- regionofinterest + 1000*length(regionofinterest) }
if( jump == jumpschoices[5] ){ regionofinterest <- regionofinterest - length(regionofinterest) }
if( jump == jumpschoices[6] ){ regionofinterest <- regionofinterest - 10*length(regionofinterest) }
if( jump == jumpschoices[7] ){ regionofinterest <- regionofinterest - 100*length(regionofinterest) }
if( jump == jumpschoices[8] ){ regionofinterest <- regionofinterest - 1000*length(regionofinterest) }
return(regionofinterest)
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