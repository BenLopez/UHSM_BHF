# Script to hold function for peak extraction

RPeakExtraction <- function( Times , ECGWaveFormData ){
  # Function to extract location of R peak times and amplitudes from ECG data.
  # The idea behind this approach is a heuristic attempt to finding the closest point
  # to where the second derivative is zero. "similar toR-peak detection algorithm for ECG using double difference and RR
  #  interval processing"
  
  # Inputs: time series of ECG data  
  # Ouputs: dataframe[peak times , peak amplitudes , peak R-R times]  
  
  # Take the approximate first derivtaive with respect to time.
  d_ECGWaveFormData_dt <- c(0,diff(ECGWaveFormData))
  # Find locations where the derivative is in the tail of the derivative distribution. 
  Rlogical <- abs(d_ECGWaveFormData_dt) > (mean(abs(d_ECGWaveFormData_dt)) + 3*sqrt(var(abs(d_ECGWaveFormData_dt)))) 
  
  # find locations where the sign changes
  diffRlogical <- c(c(diff(d_ECGWaveFormData_dt[Rlogical]),0) < 0)
  
  # Extract locations from data
  Rtimes <- Times[Rlogical]
  RAmplitudes <- ECGWaveFormData[Rlogical]
  
  # Extract locations where dervative changes sign
  Rtimes2 <- Rtimes[diffRlogical]
  RAmplitudes2 <- RAmplitudes[diffRlogical]
  
  # Clean data by removing locations where the difference between Rpeaks is too small. 
  diffRtimes2 = c(diff(Rtimes2), 1)
  Rtimes2 <- Rtimes2[diffRtimes2 > 0.1] 
  RAmplitudes2 <- RAmplitudes2[diffRtimes2 > 0.1]
  
  t  <- Rtimes2
  RA <- RAmplitudes2
  RR <- c(diff(Rtimes2) , mean(diff(Rtimes2)))
  
  output <- data.frame(t , RA , RR)
  
  return(output)
}
RPeakExtractionWavelet <- function(WaveData , Filter , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5 ){
  # Function to detect R-peaks and RR wave times using a wavelet decomposition. 
  
  # Inputs : WaveData[date , wavefirm], Filter structure for wavelet package  
  # Ouputs: dataframe[peak times , peak amplitudes , peak R-R times]  
  
  # Extract data from dataframe  
  f_t <- WaveData$Value
  t <- WaveData$Date
  
  imodoutput <- waveletremovecompenentsandreconstruct(f_t , Filter , nlevels , ComponetsToKeep )
  
  # Find local maxima 
  stdresid = imodoutput/sqrt(var(imodoutput))
 
  # Only look for positive peaks
  stdresid[ stdresid < 0] = 0 
  Peakslogical = stdresid > stdthresh
  
 #Peakslogical[f_t < 0] = FALSE
  Peakslogical <- FindLocalTurningPoints(Peakslogical , f_t , maxormin = 1)
  
  RPeakLocations <- t[Peakslogical == TRUE]
  RAmplitudes <- f_t[Peakslogical == TRUE ]
  
  t  <- RPeakLocations
  RA <- RAmplitudes
  RR <- c(diff(RPeakLocations) , mean(diff(RPeakLocations)))
  
  output <- data.frame(t , RA , RR)
  return(output)
}
FindLocalTurningPoints <- function(Peakslogical , f_t , maxormin = 1){
  # Function to find max in sets of logicals in a times series
  # Inputs: Peakslogical, a series of logicals definning the local regions with turning points.
  # f_t, the series to to find turing point.  maxormin a logical to dictate whether max or min is found.
  for (i in 1:(length(Peakslogical) -1))
  {
    if(Peakslogical[i] == FALSE){ next }
    
    if(Peakslogical[i] == TRUE && Peakslogical[i + 1] == TRUE )
    {# Count number of conected logicals
      j <- 1
      while(Peakslogical[min(i + j , length(Peakslogical))] ==  TRUE & ((i + j) <  length(Peakslogical)) ){ j <- (j+1) }
    }
    if(Peakslogical[i] == TRUE && Peakslogical[i + 1] == FALSE){j <- 0}
    # Find location of max of f_t in set of connected logicals
    
    if(maxormin == 1){maxindex <- which.max(f_t[i:(i+j)])}
    if(maxormin == 0){maxindex <- which.min(f_t[i:(i+j)])}
    Peakslogical[i:(i+j)] <- FALSE
    Peakslogical[i + (maxindex -1)] <- TRUE
  }
return(Peakslogical)
}
FindQRSComplex <- function( WaveData , Filter , nlevels = 12 , ComponetsToRemove = c(3,4) , Filter2 = rep(1 , 41) ){
# Function to find QRS complex from an ECG.   
# Inputs :   
  
  tt <- WaveData$Date
  f_tt <- WaveData$Value
  
  imodoutput <- waveletremovecompenentsandreconstruct(f_tt , Filter)
  
  stdresid <- imodoutput/sqrt(var(imodoutput))
  stdresid2 <- stdresid
  
  stdresid[ stdresid < 0] <- 0
  Rlogical <- FindLocalTurningPoints( stdresid > 2.8 , f_tt )
  RPeakLocations <- tt[ Rlogical == TRUE ]
  RAmplitudes <- f_tt[ Rlogical == TRUE ]  
  
  stdresid2[stdresid2 > 0] = 0
  stdresid2 <- abs( stdresid2 )
  stdresid2 <- (stdresid2 > ( mean(stdresid2) + 2*sqrt( var(stdresid2)) ))
  
  # Use a filter to create a local clique for every R peak.
  Rregion <- ((imfilter1D( Rlogical , Filter2)) < 0.5)
  stdresid2[Rregion] = 0
  
  QSlogical <- FindLocalTurningPoints(stdresid2 , f_tt , 0)
  QSlocations <- tt[QSlogical == TRUE]
  QSValues <- f_tt[QSlogical == TRUE]
  
  Temp <- get.knnx(QSlocations , RPeakLocations , k = 2)
  Temp$nn.index <- t(apply(Temp$nn.index, 1 , sort))
  QLocations <- QSlocations[Temp$nn.index[ , 1]]
  SLocations <- QSlocations[Temp$nn.index[ , 2]]
  
  QAmplitudes <- QSValues[Temp$nn.index[ , 1]]
  SAmplitudes <- QSValues[Temp$nn.index[ , 2]]
  QStime <- t(apply(Temp$nn.dist , 1 , sum))
  rm(Temp)
  
  Rt <- RPeakLocations
  RA <- RAmplitudes
  RR <- c(diff(RPeakLocations) , mean(diff(RPeakLocations)))
  Qt <- QLocations
  St <- SLocations
  QA <- QAmplitudes
  SA <- SAmplitudes
  QS <- QStime
  
  output <- data.frame(Rt ,RA , RR , Qt , St , QA , SA ,  QS ) 
  return(output)
  
}
waveletremovecompenentsandreconstruct <- function(f_tt , Filter , nlevels = 12 , ComponetsToKeep = c(3,4)){
  # Function to perform a wavelet decompostion remove components and reconstruct.
  # Inputs: signal, filter structure, number of levels, components to remove.
  # otuputs: Reconstructed signal
  a <- c(1:nlevels)
  a <- a[ -ComponetsToKeep ]
  modoutput <-  modwt( f_tt  , Filter , nlevels)
  W <- slot(modoutput , 'W')
  V <- slot(modoutput , 'V')
  W <- SetElementsOfListoToZero(W , a )
  V <- SetElementsOfListoToZero(V , a )
  slot(modoutput , 'W')<- W
  slot(modoutput , 'V')<- V
  imodoutput <- imodwt(modoutput, fast=TRUE)
  return(imodoutput)
}
SeparateWaveQRSandPTWaveforms <- function(WaveData , QRSoutput , bandincrement = (10*0.005)){
# Function to separate QRS and PT Wave Forms  
# Inputs: WaveData QRS complex locations and amplitudes, band incriment is a logical cliques round QRS complexes
  
  tt <- WaveData$Date
  f_tt <- WaveData$Value
  
  # Preallocate matrix
  QRSLogical <- matrix(0 , length(tt) , 1)
  for(i in 1:length( QRSoutput$Qt ))
  {
    QRSLogical[ ((tt >= ( QRSoutput$Qt[i] - bandincrement ) )*((tt <= (QRSoutput$St[i] + bandincrement)  ))) == 1 ] <- 
      1 +  QRSLogical[ ((tt >= (QRSoutput$Qt[i] - bandincrement)  )*((tt <= (QRSoutput$St[i] + bandincrement)  ))) == 1 ]
  }
  
  QRSWaveForm <- matrix(0 , length(tt) , 1)
  QRSWaveForm[QRSLogical == 1] <- f_tt[QRSLogical == 1]
  QRSWaveForm[QRSLogical == 0] <- NA
  QRSWaveForm[1] <- 0
  QRSWaveForm[length(QRSWaveForm)] <- 0
  QRSWaveForm <- na.approx(QRSWaveForm)
  
  PTWaveFrom <- matrix(0 , length(tt) , 1)
  PTWaveFrom[QRSLogical == 0] <- f_tt[QRSLogical == 0]
  PTWaveFrom[QRSLogical == 1] <- NA
  PTWaveFrom[1] <- 0
  PTWaveFrom[length(PTWaveFrom)] <- 0
  PTWaveFrom <- na.approx(PTWaveFrom)
  
  output <- data.frame(QRSWaveForm , PTWaveFrom)
  return(output)
    
}
ExtractPQRST <- function( WaveData , Filter , nlevels = 12 , ComponetsToKeep = c(3,4) , ComponetsToKeep2 = 7 , ComponetsToKeep3 = 6 ){
# A function to extract PQRST peaks for ECG waveform data.   
# Inputs: WaveData, a filter struture, nlevels for wavelet decomposition, compents to keep for WRS detection, compnets to keep for P detection Components to keep for PT detection 
  
tt <- WaveData$Date
f_tt <- WaveData$Value 
   
QRSoutput <- FindQRSComplex(WaveData , Filter)
SeparatedWaveForms <- SeparateWaveQRSandPTWaveforms( WaveData , QRSoutput)

# Find Pt locations
output <- waveletremovecompenentsandreconstruct(SeparatedWaveForms$PTWaveFrom , Filter  , 12 , ComponetsToKeep3)
stdresid <- output/sqrt(var(output))
stdresid[stdresid<0] = 0
PTLocations  <-  FindLocalTurningPoints(stdresid>mean(stdresid) , SeparatedWaveForms$PTWaveFrom , 1)
PTAmplitudes <-  f_tt[PTLocations]
PTLocations  <-  tt[PTLocations]

# Find T locations
output <- waveletremovecompenentsandreconstruct(SeparatedWaveForms$PTWaveFrom , Filter  , nlevels , ComponetsToKeep2)
stdresid <- output/sqrt(var(output))
stdresid[stdresid<0] = 0
TLocations <- FindLocalTurningPoints(stdresid>mean(stdresid) , SeparatedWaveForms$PTWaveFrom , 1)
TAmplitudes <- f_tt[TLocations]
TLocations <- tt[TLocations]

# Remove T locations to 
Plocations <- PTLocations
PAmpltidues <- PTAmplitudes

for (i in 1:length(TLocations)){
index <- which.min( abs(as.numeric(Plocations) - as.numeric(TLocations[i]) ) )
Plocations[index] <- NA
PAmpltidues[index] <- NA
}

Plocations <- Plocations[is.na(Plocations) == FALSE]
PAmpltidues <- PAmpltidues[is.na(PAmpltidues) == FALSE]

Rt <- QRSoutput$Rt
RA <- QRSoutput$RA
RR <- QRSoutput$RR
Qt <- QRSoutput$Qt
St <- QRSoutput$St
QA <- QRSoutput$QA
SA <- QRSoutput$SA
Tt <- TLocations
TA <- TAmplitudes
Pt <- Plocations
PA <- PAmpltidues

output <- list(Rt , RA , RR , Qt , St , QA , SA , Tt , TA , Pt , PA  )
outputnames <- c('Rt' , 'RA' , 'RR' , 'Qt' , 'St' , 'QA' , 'SA' , 'Tt' , 'TA' , 'Pt' , 'PA')
output <- setNames( output , outputnames) 
  
return(output)
}
ReturnWaveformwithPositiveOrientation <- function(WaveData){  
qus <- abs(quantile(WaveData$Value, probs = c(0.01 , 0.99) , na.rm = TRUE ))
if( as.numeric(qus[1]) > as.numeric(qus[2]) )
{
  Date   <-   WaveData$Date 
  Value  <-  -WaveData$Value
  WaveData <- data.frame(Date , Value)
}
output<-WaveData
return(output)
}
CleanRpeaks <- function( RWaveExtractedData , threshold = 2 ){
  
  # Function to clean R peaks by removing outliers caused by data gaps and missed beats.
  RWaveExtractedData$RR[RWaveExtractedData$RR > threshold] <- median(RWaveExtractedData$RR)
  RWaveExtractedData$RR[RWaveExtractedData$RR < 0] <- median(RWaveExtractedData$RR)
  
  m <- smth(RWaveExtractedData$RR , method = 'sma' , n = 100)
  m[is.na(m)] <- mean( m  , rm.na = TRUE)
  mm <- abs(RWaveExtractedData$RR - m)
  
  m <- smth(mm , method = 'sma' , n = 100)
  m[is.na(m)] <- mean( m  , rm.na = TRUE)
  
  v <- smth(mm , method = 'sma' , n = 100) - m^2
  v[is.na(m)] <- mean( v , rm.na = TRUE)
  stdresid <- abs((mm - m)/sqrt(v)) 
    
  RWaveExtractedData$RR[stdresid > 2] <- median(RWaveExtractedData$RR)
  
  return(RWaveExtractedData)
}
ExtractIHVAFScore <- function( RWaveExtractedData , binlims = c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3) , n = 250 ){
  binmatrix <- CalulateBinMatrix(RWaveExtractedData , binlims , n )
  t <- RWaveExtractedData[!is.na(binmatrix[ , 1]) ,1] 
  binmatrix <-  binmatrix[!is.na(binmatrix[ , 1]),]
  output <- setNames(data.frame(t , 1/apply(binmatrix , 1 , var)) , c('t' , 'IHAVFScore'))
  return( output )
}
CalulateBinMatrix <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.8  , 0.05  )), n = 250){
  binmatrix <- matrix(0 , length( RWaveExtractedData$RR ) , length(binlims) - 1 )
  
  for(i in 1:(length(binlims) - 1))
  {
    binmatrix[ , i] <- smth( (RWaveExtractedData$RR > binlims[i])*(RWaveExtractedData$RR <= binlims[i+1])  , method = 'sma'  , n = n)
  }
  return( binmatrix )
}
ExtractNumberofModes <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.5  , 0.025  )), n = 250 , densitythresh = 0.025){
  # Function to calulate the number of modes in a local region
  
  binmatrix <- CalulateBinMatrix(RWaveExtractedData , binlims = binlims ,  n = n) 
  tbin <- RWaveExtractedData[!is.na(binmatrix[ , 1]) ,1] 
  binmatrix <-  binmatrix[!is.na(binmatrix[ , 1]),]
  NumberModes <- apply(binmatrix , 1 , function(X){ sum( FindLocalTurningPoints( X > densitythresh ,  1:length( binmatrix[1 , ]) ) )})
  
  return(data.frame(t = tbin, NumModes = NumberModes))

}
Calculatemodalmode <- function(RWaveExtractedData , binlims= c(0, seq(from = 0.25  , to = 1.5  , 0.025  )), n = 250 , densitythresh = 0.025 , nn = 500){

  output <- ExtractNumberofModes(RWaveExtractedData , binlims, n  , densitythresh )
  output$NumModes <- round( smth( output$NumModes , method = 'sma' , n = nn ) )
  output$NumModes[1:nn] <-output$NumModes[(nn+1)]
  output$NumModes[(length(output$NumModes) - nn) : length(output$NumModes)] <- output$NumModes[length(output$NumModes) - nn -1 ]
  return(output)
}
PE_MultipleECGRPeaks <- function( outputdata , thresh = 0.01 ){
    # Calulate distance to closest peak 
  mindistancesI <- apply( as.matrix( as.numeric(outputdata$ECGIII$t) ) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGI$t)))} )
  mindistancesII <- apply( as.matrix( as.numeric(outputdata$ECGIII$t) ) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGII$t)))} )
  
  # Find good peaks
  Goodlogical <- (mindistancesI <= thresh) + (mindistancesII <= thresh)
  
  # Extract good peaks
  tt_good <- outputdata$ECGIII$t[Goodlogical > 0]
  tt_notgood <- outputdata$ECGIII$t[Goodlogical == 0]
  
  # Set of beats in ECGII which are classified as good.
  goodlogical_2 <- apply( as.matrix(as.numeric(outputdata$ECGII$t)) , 1 , function(X){min(abs(X - as.numeric(tt_good))) <= thresh} )
  
  ttmp <- outputdata$ECGII$t[ goodlogical_2 == 0 ]
  
  mindistancesI <- apply( as.matrix( as.numeric(ttmp) ) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGI$t)))} )
  mindistancesII <- apply( as.matrix( as.numeric(ttmp) ) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGIII$t)))} )
  
  Goodlogical <- (mindistancesI <= thresh) + (mindistancesII <= thresh)
  
  tt_good <- c( tt_good , outputdata$ECGII$t[Goodlogical > 0]  )
  tt_notgood <- c( tt_notgood , outputdata$ECGII$t[Goodlogical == 0] )
  
  goodlogical_1 <- apply( as.matrix(as.numeric(outputdata$ECGI$t)) , 1 , function(X){min(abs(X - as.numeric(tt_good))) <= thresh} )
  ttmp <- outputdata$ECGI$t[goodlogical_1 == 0 ]
  
  mindistancesI <- apply( as.matrix( as.numeric(ttmp) ) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGII$t)))} )
  mindistancesII <- apply( as.matrix( as.numeric(ttmp) ) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGIII$t)))} )
  
  Goodlogical <- (mindistancesI <= thresh) + (mindistancesII <= thresh)
  
  tt_good <- sort(c(tt_good , outputdata$ECGI$t[Goodlogical > 0] ))
  tt_notgood <- c(tt_notgood , outputdata$ECGI$t[Goodlogical == 0] )
  
  RR <- c(0 , diff(tt_good))
  tt_good <- tt_good[((RR > 0.1)*(RR < 4)) == 1]
  RR <- RR[((RR > 0.1)*(RR < 4)) == 1]
  
  ECG <- data.frame(t = tt_good , RA = 150 + 0*RR , RR = RR)
}
  


