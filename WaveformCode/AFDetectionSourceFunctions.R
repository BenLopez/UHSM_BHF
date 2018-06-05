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

AFD_ExtractIHVAFScore <- function( RWaveExtractedData , binlims = c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3) , n = 250 ){
  binmatrix <- CalulateBinMatrix(RWaveExtractedData , binlims , n )
  t <- RWaveExtractedData[!is.na(binmatrix[ , 1]) ,1] 
  binmatrix <-  binmatrix[!is.na(binmatrix[ , 1]),]
  output <- setNames(data.frame(t , 1/apply(binmatrix , 1 , var)) , c('t' , 'IHAVFScore'))
  return( output )
}

AFD_ExtractNumberofModes <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.5  , 0.025  )), n = 250 , densitythresh = 0.025){
  # Function to calulate the number of modes in a local region
  
  binmatrix <- CalulateBinMatrix(RWaveExtractedData , binlims = binlims ,  n = n) 
  tbin <- RWaveExtractedData[!is.na(binmatrix[ , 1]) ,1] 
  binmatrix <-  binmatrix[!is.na(binmatrix[ , 1]),]
  NumberModes <- apply(binmatrix , 1 , function(X){ sum( FindLocalTurningPoints( X > densitythresh ,  1:length( binmatrix[1 , ]) ) )})
  
  return(data.frame(t = tbin, NumModes = NumberModes))
  
}

AFD_Calculatemodalmode <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.5  , 0.025  )), n = 250 , densitythresh = 0.025 , nn = 500){
  
  output <- ExtractNumberofModes(RWaveExtractedData , binlims, n  , densitythresh )
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

