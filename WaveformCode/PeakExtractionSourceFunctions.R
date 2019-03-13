# Script to hold function for peak extraction

PE_FindLocalTurningPoints <- function(Peakslogical , f_t , maxormin = 1){
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
    if(maxormin == 3){maxindex <- which.max(abs(f_t[i:(i+j)]))}
    Peakslogical[i:(i+j)] <- FALSE
    Peakslogical[i + (maxindex -1)] <- TRUE
  }
  return(Peakslogical)
}
PE_ReformulateUniqueTimeDataFrame <- function(TableofUniqueTimes){
  UniqueDays    <- unique( round.POSIXt(TableofUniqueTimes$Var1 , units = 'days') )
  UniqueHours   <- unique( format(TableofUniqueTimes$Var1 , format = '%H') )
  UniqueMinutes <- unique( format(TableofUniqueTimes$Var1 , format = '%M') )
  UniqueSeconds <- unique( format(TableofUniqueTimes$Var1 , format = '%S') )
  
  SecsStruct <- list()
  for(i in 1:length(UniqueSeconds)){
    SecsStruct[[i]]  <- 1
  }
  names(SecsStruct) <- UniqueSeconds
  MinsStruct <- list()
  for(i in 1:length(UniqueMinutes)){
    MinsStruct[[i]]  <- SecsStruct
  }
  names(MinsStruct) <- UniqueMinutes
  HoursStruct <- list()
  for(i in 1:length(UniqueHours)){
    HoursStruct[[i]] <- MinsStruct
  }
  names(HoursStruct) <- UniqueHours
  daysstruct <- list()
  for(i in 1:length(UniqueDays)){
    NewStruct[[i]] <- HoursStruct
  }
  names(NewStruct) <- UniqueDays
  
  for (i in 1:length( UniqueDays ) ){
    days <- TableofUniqueTimes[ round.POSIXt(TableofUniqueTimes$Var1 , units = 'days') == UniqueDays[i], ]
    #NewStruct[[i]][[j]]
    for(j in 1:length( UniqueHours )){
      hours <- days[ format(days$Var1 , format = '%H') == UniqueHours[j] ,]
      #NewStruct[[i]][[j]]
      for( k in 1:length(UniqueMinutes) ){
        #NewStruct[[i]][[j]][[k]]
        minutes <- hours[ format(hours$Var1 , format = '%M') == UniqueMinutes[k], ]
        for(l in 1:length(UniqueSeconds)){
          if(nrow(minutes[ format(minutes$Var1 , format = '%S') == UniqueSeconds[l], ]) > 0){
            NewStruct[[i]][[j]][[k]][[l]] <- minutes[ format(minutes$Var1 , format = '%S') == UniqueSeconds[l], ]}
          if(nrow(minutes[ format(minutes$Var1 , format = '%S') == UniqueSeconds[l], ]) == 0){
            NewStruct[[i]][[j]][[k]][[l]] <- data.frame(Var1=NA , Freq= NA)  
          }
        }
      }
    }
  }
  
  return(NewStruct)
}
PE_ReformulateUniqueTimeDataFrame2 <- function(TableofUniqueTimes){
  UnClassed <- unclass(TableofUniqueTimes$Var1)
  uniquemonths <- unique(UnClassed$mon)
  uniquedays <- unique(UnClassed$mday)
  uniquehours <- unique(UnClassed$hour)
  uniqueminutes <- unique(UnClassed$min)
  uniqueseconds <- unique(UnClassed$sec)
  
  NewStruct <- list()
  for(ii in 1:length(uniquemonths)){
    NewStruct[[ii]] <- list()
    for(jj in 1:length(uniquedays)){
      NewStruct[[ii]][[jj]]<-list()
      for(kk in 1:length(uniquehours)){
        NewStruct[[ii]][[jj]][[kk]]<-list()
        for(ll in 1:length(uniqueminutes) ){
          NewStruct[[ii]][[jj]][[kk]][[ll]] <- TableofUniqueTimes[
            (UnClassed$mon ==    uniquemonths[ii])*
              (UnClassed$mday ==   uniquedays[jj])*
              (UnClassed$hour ==   uniquehours[kk])*
              (UnClassed$min ==    uniqueminutes[ll]) == 1 ,]
        }
        names(NewStruct[[ii]][[jj]][[kk]]) <- uniqueminutes
      }
      names(NewStruct[[ii]][[jj]]) <- uniquehours
    }
    names(NewStruct[[ii]]) <- uniquedays
  }
  
  names(NewStruct) <- uniquemonths
  
  
  return(NewStruct)
  
}
PE_ExractECGActiveLogical <- function(RpeaksFile , ECG1 , ECG2){
  output <- matrix(0 , length(RpeaksFile$t) , 2)
  output[( RpeaksFile$t %in% ECG1 ) , 1] <- 1 
  output[( RpeaksFile$t %in% ECG2 ) , 2] <- 1
  return(output)  
}
PE_ReturnWaveformwithPositiveOrientation <- function(WaveData){  
  mm <- rollmedian(WaveData$Value ,  k = 21 , na.pad = TRUE)
  mm[is.na(mm)] <- 0
  qus <- abs(quantile(WaveData$Value - mm, probs = c(0.01 , 0.99) , na.rm = TRUE ))
  if( as.numeric(qus[1]) > as.numeric(qus[2]) )
  {
    Date   <-   WaveData$Date 
    Value  <-  -WaveData$Value
    WaveData <- data.frame(Date , Value)
  }
  output <- WaveData
  return(output)
}
PE_ExtractOrientationLogical <- function(WaveData , theshold = 10 , k1 = 21 , k2 = 41){
  tmp <- PE_ReturnWaveformwithPositiveOrientation(WaveData)$Value
  mm <-rollmedian(tmp ,  k = k1 , na.pad = TRUE)
  mm[is.na(mm)] <- 0
  tmp <- tmp - mm
  ECGMax <- rollmax(tmp , align = 'center' , na.pad = TRUE, k = k2)
  ECGMax[is.na(ECGMax)] <- 0
  ECGmin <- -rollmax(-(tmp) ,  align = 'center' , na.pad = TRUE ,  k = k2)
  ECGmin[is.na(ECGmin)] <- 0
  Oreintationlogical <- ((abs(ECGMax) + theshold) > (abs(ECGmin)))
  Oreintationlogical2 <- ((abs(ECGmin) + theshold ) > (abs(ECGMax)))
  
  return(data.frame(Oreintationlogical1 = Oreintationlogical , Oreintationlogical2 = Oreintationlogical2 ))  
}
PE_RPeakExtractionWavelet <- function(WaveData , Filter = wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5 ){
  # Function to detect R-peaks and RR wave times using a wavelet decomposition. 
  
  # Inputs : WaveData[date , wavefirm], Filter structure for wavelet package  
  # Ouputs: dataframe[peak times , peak amplitudes , peak R-R times]  
  #WaveData <- PE_ReturnWaveformwithPositiveOrientation(WaveData)
  # Extract data from dataframe  
  f_t <- WaveData$Value
  t <- WaveData$Date
  
  imodoutput <- waveletremovecompenentsandreconstruct(f_t , Filter , nlevels , ComponetsToKeep )
  
  # Find local maxima 
  vv <- VO_rollvar(imodoutput)
  vv[(is.na(vv))] <- mean(vv[!is.na(vv)])
  vv[vv == 0] <- mean(vv[!is.na(vv)])
  
  stdresid = imodoutput/sqrt(vv)
  
  
  Orientationlogical <- PE_ExtractOrientationLogical(WaveData = WaveData )
  
  # Only look for positive peaks
  stdresid[Orientationlogical$Oreintationlogical1 == 0] <- -stdresid[Orientationlogical$Oreintationlogical1 == 0]
  stdresid[stdresid < 0 ] = 0
  Peakslogical = abs(stdresid) > stdthresh
  
  #Peakslogical[f_t < 0] = FALSE
  Peakslogical <- PE_FindLocalTurningPoints(Peakslogical , f_t , maxormin = 3)
  
  RPeakLocations <- t[Peakslogical == TRUE]
  RAmplitudes <- f_t[Peakslogical == TRUE ]
  
  t  <- RPeakLocations
  RA <- RAmplitudes
  RR <- c(diff(RPeakLocations) , mean(diff(RPeakLocations)))
  
  output <- data.frame(t , RA , RR)
  return(output)
}
PE_MultipleECGRPeaks <- function( outputdata , ECGs ,  thresh = 0.02  ){
  
  ActiveMatrix <- rbind(1 + as.matrix(apply(PE_ExractECGActiveLogical(outputdata$ECGI , ECGs$ECGII$Date, ECGs$ECGIII$Date) , 1 , sum)) ,
                        1 + as.matrix(apply(PE_ExractECGActiveLogical(outputdata$ECGII , ECGs$ECGI$Date, ECGs$ECGIII$Date) , 1 , sum)) ,
                        1 + as.matrix(apply(PE_ExractECGActiveLogical(outputdata$ECGIII , ECGs$ECGII$Date, ECGs$ECGI$Date) , 1 , sum)))
  PeakMatrix <- c(outputdata$ECGI$t , outputdata$ECGII$t , outputdata$ECGIII$t )
  
  ActiveMatrix <- ActiveMatrix[order(PeakMatrix)]
  PeakMatrix <- sort(PeakMatrix)
  
  MinDIsMatrix <- cbind( as.matrix(apply(as.matrix(PeakMatrix) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGI$t)) )} )) ,
                         as.matrix(apply(as.matrix(PeakMatrix) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGII$t)) )} )) ,
                         as.matrix(apply(as.matrix(PeakMatrix) , 1 , function(X){min(abs(X - as.numeric(outputdata$ECGIII$t)) )} ))  )
  
  NumberofECGsPeaksAppearin <- apply( MinDIsMatrix <= thresh , 1 , sum  )
  
  # deal with active matrix
  NumberofECGsPeaksAppearin[ActiveMatrix == 1] = 3 
  Goodlogical = (NumberofECGsPeaksAppearin >= 2)
  Goodlogical[duplicated(PeakMatrix) ] = F
  Goodlogical[c(0,diff(PeakMatrix) ) <= thresh] = F
  
  t <- PeakMatrix[Goodlogical]
  RR <- c( 0 , diff(PeakMatrix[Goodlogical]))
  t <- t[((RR > 0.1)*(RR < 4)) == 1]
  RR <- RR[((RR > 0.1)*(RR < 4)) == 1]
          
          return( data.frame(t = t , RA = 150 + 0*RR , RR = RR) )        
          
}
PE_CleanRpeaks <- function( RWaveExtractedData , threshold = 2 ){
  
  # Function to clean R peaks by removing outliers caused by data gaps and missed beats.
  RWaveExtractedData$RR[RWaveExtractedData$RR > threshold] <- median(RWaveExtractedData$RR)
  RWaveExtractedData$RR[RWaveExtractedData$RR < 0.25] <- median(RWaveExtractedData$RR)
  
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


##### Mischananious functions and legacy functions #####

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
RPeakExtractionWavelet <- function(WaveData , Filter , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5){
  return(PE_RPeakExtractionWavelet(WaveData , Filter , nlevels , ComponetsToKeep  , stdthresh ))
}
FindLocalTurningPoints <- function(Peakslogical , f_t , maxormin = 1){
  return(PE_FindLocalTurningPoints(Peakslogical , f_t , maxormin))
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
  return(PE_ReturnWaveformwithPositiveOrientation(WaveData))
}
CleanRpeaks <- function( RWaveExtractedData , threshold = 2 ){
  return(PE_CleanRpeaks(RWaveExtractedData , threshold))
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
