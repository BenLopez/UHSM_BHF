
ECGSim_LaplaceDisFunction <- function(x , mu , b){
  (1/(2*b))*exp(-abs(x-mu)/b)
}

ECGSim_Rpeak <- function(x , mu , b){
  (2*b)*ECGSim_LaplaceDisFunction(x , mu , b)
}

ECGSim_Gaussian <- function(x , mu , sigma){
  exp( -0.5*((x - mu)^2)/(sigma^2)   )
}  

ECGSim_SkewGaussian <- function(x , mu , sigma , slant)
{
  tmp <- dsn(x ,xi = mu, omega = sigma , alpha = slant)
  index <- which.max(tmp) 
  if(tmp[index]  == 0){return(tmp)}  
  if(tmp[index] != 0){
  tmp <- tmp/tmp[index]
  index2 <- which.min(abs(x-mu))
  l <- index - index2
  
  if(sign(l) == -1 ){ output <- c( tmp[(length(x) - l):length(x)]  , tmp[1:(length(x) - abs(l) -1)]) }
  if(sign(l) == 1  ){ output <- c( tmp[(l+1):length(x)] , tmp[1:l] ) }
  if(sign(l) == 0  ){ output <- tmp }
  return( output )
   }
}

ECGSim_WrapperSingleBeat <- function(x , t_observation){
  
  Baseline <- x[1]
  Rcen <- x[2]
  Rwidth <- x[3]
  RA <- x[4]
  Qcen <- x[5]
  Qwidth <- x[6]
  QA <- x[7]
  Scen <-   x[8]
  Swidth <- x[9]
  SA <- x[10]
  Pcen <- x[11]
  Pwidth <- x[12]
  PA <- x[13]
  Pslant <- x[14]
  Tcen <- x[15]
  Twidth <- x[16]
  TA <- x[17]
  Tslant <- x[18]
  
  return(Baseline + RA*ECGSim_Rpeak( t_observation , Rcen , Rwidth ) +
           QA*ECGSim_Rpeak( t_observation , Qcen , Qwidth  ) +
           SA*ECGSim_Rpeak( t_observation , Scen , Swidth  ) +
           PA*ECGSim_SkewGaussian(x=t_observation , mu=Pcen , sigma=Pwidth , slant =Pslant ) +
           TA*ECGSim_SkewGaussian(t_observation , Tcen , Twidth , Tslant) )
}

