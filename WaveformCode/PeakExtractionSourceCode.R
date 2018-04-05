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