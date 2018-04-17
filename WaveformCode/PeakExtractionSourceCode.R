# Script to hold function for peak extraction

RPeakExtraction <- function( Times , ECGWaveFormData )
{
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

RPeakExtractionWavelet <- function(WaveData , Filter , nlevels = 12 , ComponetsToRemove = c(3,4) )
{
  # Function to detect R-peaks and RR wave times using a wavelet decomposition. 
  
  # Inputs : WaveData[date , wavefirm], Filter structure for wavelet package  
  # Ouputs: dataframe[peak times , peak amplitudes , peak R-R times]  
  
  # Extract data from dataframe  
  f_t <- WaveData$Value
  t <- WaveData$Date
  
  a <- c(1:nlevels)
  a <- a[ -ComponetsToRemove ]
  
  # Perform discrete wavelet transform on data.
  modoutput <-  modwt( f_t , Filter , nlevels)
  
  # set components to zero to supress unwanted data features.
  W <- slot(modoutput , 'W')
  V <- slot(modoutput , 'V')
  W <- SetElementsOfListoToZero( W , a )
  V <- SetElementsOfListoToZero( V , a )
  slot( modoutput , 'W' ) <- W
  slot( modoutput , 'V' ) <- V
  
  # Perform inverse discrete wavelet tramsform to recronstruct cleaned signal.
  imodoutput <- imodwt(modoutput, fast=TRUE)
  
  # Find local maxima 
  stdresid = imodoutput/sqrt(var(imodoutput))
  
  # Only look for positive peaks
  stdresid[ stdresid < 0] = 0 
  Peakslogical = stdresid>2.8
  
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
    maxindex <- which.max(f_t[i:(i+j)])
    Peakslogical[i:(i+j)] <- FALSE
    Peakslogical[i + (maxindex -1)] <- TRUE
  }
  
  RPeakLocations <- t[Peakslogical == TRUE]
  RAmplitudes <- f_t[Peakslogical == TRUE ]
  
  t  <- RPeakLocations
  RA <- RAmplitudes
  RR <- c(diff(RPeakLocations) , mean(diff(RPeakLocations)))
  
  output <- data.frame(t , RA , RR)
  return(output)
}

FindLocalTurningPoints <- function(Peakslogical , f_t , maxormin = 1)
{
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
