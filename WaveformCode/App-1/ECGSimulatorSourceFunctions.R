ECGSim_LaplaceDisFunction <- function(x , mu , b){
  (1/(2*b))*exp(-abs(x-mu)/b)
}
ECGSim_Rpeak <- function(x , mu , b){
  (2*b)*ECGSim_LaplaceDisFunction(x , mu , b)
}
ECGSim_Gaussian <- function(x , mu , sigma){
  exp( -0.5*((x - mu)^2)/(sigma^2)   )
}  
ECGSim_DiffGaussian <- function(x , mu , sigma){
  -((x - mu)/(sigma))*exp( -0.5*((x - mu)^2)/(sigma^2)   )
}
ECGSim_SkewGaussian <- function(x , mu , sigma , slant){
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
  
  return(Baseline + RA*ECGSim_Gaussian( t_observation , Rcen , Rwidth ) +
           QA*ECGSim_Gaussian( t_observation , Qcen , Qwidth  ) +
           SA*ECGSim_Gaussian( t_observation , Scen , Swidth  ) +
           PA*ECGSim_SkewGaussian(x=t_observation , mu=Pcen , sigma=Pwidth , slant =Pslant ) +
           TA*ECGSim_SkewGaussian(t_observation , Tcen , Twidth , Tslant) )
}
ECGSim_WrapperSingleBeat_m2 <- function(x , t_observation){
  
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
  PcenL <- x[11]
  PwidthL <- x[12]
  PAL <- x[13]
  PcenR <- x[14]
  PwidthR <- x[15]
  PAR <- x[16]
  TcenL <- x[17]
  TwidthL <- x[18]
  TAL <- x[19]
  TcenR <- x[20]
  TwidthR <- x[21]
  TAR <- x[22]
  
  
  return(Baseline + RA*ECGSim_Gaussian( t_observation , Rcen , Rwidth ) +
           QA*ECGSim_Gaussian( t_observation , Qcen , Qwidth  ) +
           SA*ECGSim_Gaussian( t_observation , Scen , Swidth  ) +
           PAL*ECGSim_Gaussian(t_observation , PcenL , PwidthL ) +
           PAR*ECGSim_Gaussian(t_observation , PcenR , PwidthR ) +
           TAR*ECGSim_Gaussian(t_observation , TcenR , TwidthR ) +
           TAL*ECGSim_Gaussian(t_observation , TcenL , TwidthL ) )
}  
ECGSim_LeftPwave <- function(x , t_observation){
  
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
  PcenL <- x[11]
  PwidthL <- x[12]
  PAL <- x[13]
  PcenR <- x[14]
  PwidthR <- x[15]
  PAR <- x[16]
  TcenL <- x[17]
  TwidthL <- x[18]
  TAL <- x[19]
  TcenR <- x[20]
  TwidthR <- x[21]
  TAR <- x[22]

  return(PAL*ECGSim_Gaussian(t_observation , PcenL , PwidthL ))
} 
ECGSim_RightPwave <- function(x , t_observation){
  
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
  PcenL <- x[11]
  PwidthL <- x[12]
  PAL <- x[13]
  PcenR <- x[14]
  PwidthR <- x[15]
  PAR <- x[16]
  TcenL <- x[17]
  TwidthL <- x[18]
  TAL <- x[19]
  TcenR <- x[20]
  TwidthR <- x[21]
  TAR <- x[22]
  
  return(PAR*ECGSim_Gaussian(t_observation , PcenR , PwidthR ))
}

ECGSim_QSSimulator <-function( X , t ){
  Baseline <- X[1]
  P_Cen    <- X[2]
  P_Amp    <- X[3] - x[1]  
  P_Width  <- X[4]
  T_Cen    <- X[5]
  T_Amp    <- X[6] - x[1]
  T_Width  <- X[7]
  
return(Baseline + 
       P_Amp*ECGSim_Gaussian( t , P_Cen , (2*P_Width)^2  ) +
       T_Amp*ECGSim_Gaussian( t , T_Cen , (2*T_Width)^2  )  )  
}
ECGSim_DiffQSSimulator <-function( X , t ){
  Baseline <- X[1]
  P_Cen    <- X[2]
  P_Amp    <- X[3] - x[1]  
  P_Width  <- X[4]
  T_Cen    <- X[5]
  T_Amp    <- X[6] - x[1]
  T_Width  <- X[7]
  
  return(P_Amp*ECGSim_DiffGaussian( t , P_Cen , (2*P_Width)^2  ) +
         T_Amp*ECGSim_DiffGaussian( t , T_Cen , (2*T_Width)^2  )  )  
}

ECGSim_QSModelDiscrepancy <- function(X , t ,  a =0 , b=1 , d =0 ,e = 0){
X[1]<-0
  tmp <- ECGSim_QSSimulator(X , t)
  return( a + b*(abs(tmp)^(-1)) + d*abs(ECGSim_DiffQSSimulator(X , t)) + e*((ECGSim_Gaussian(t, mean(t , na.rm = TRUE) , 0.4*sqrt(var(t, na.rm = TRUE))))^(-1)))
}
