CreateDefaultX <- function(mu1 = 0.6 , mu2= 0.9 , sigma1=0.01 , sigma2=0.01){
  X <- rep(0 , 10)
  
  X[1] <- 0.5
  X[2] <- 0
  X[3] <- 0.5
  X[4] <- mu1
  X[5] <- 0.0000001
  X[6] <- mu2
  X[7] <- sigma1
  X[8] <- 0.0000001
  X[9] <- sigma2
  X[10] <- 1
  return(X)
}
CreateDefaultX_PACS <- function(X1 =0.1 , X2 =0.8 , X3 = 0.1 , X4 =0.5 , X5 =0.8 , X6 =1 , X7 =0.01 , X8=0.01 , X9=0.01 , X10 =1 ){
  
  X <- rep(0 , 10)
  
  X[1] <- X1
  X[2] <- X2
  X[3] <- X3
  X[4] <- X4
  X[5] <- X5
  X[6] <- X6
  X[7] <- X7
  X[8] <- X8
  X[9] <- X9
  X[10] <- X10
  
  X[1:3] <- X[1:3]/sum(X[1:3])
  
  return(X)
  
}
PER_SampleGMMPACs <- function(X , N){
  
  # sample multinomial distribution for mixtures model
  SampleMultinomial <- t(rmultinom(N  , size = 1 , c(X[2] , X[1] + X[3]) ))
  
  output <- matrix(0 , sum(SampleMultinomial[,1]) + 2*(sum(SampleMultinomial[,2])) , 2)
  counter <- 1
  for(i in 1:dim(SampleMultinomial)[1]){
  if(SampleMultinomial[i,1] == 1){
    output[counter , 1] <- rnorm(1 , X[5] , X[8])  
    output[counter , 2] <- 1
    counter = counter +1
  }
  if(SampleMultinomial[i,1] == 0){
    output[counter , 1] <- rnorm(1 , X[4] , X[7])  
    output[counter , 2] <- 2
    output[counter + 1 , 1] <- rnorm(1 , X[6] , X[9])  
    output[counter + 1 , 2] <- 3
    counter = counter + 2
  }
  }
  return(output)
}
PER_CreateECG <- function( t , t_observation , RRTimes ){

output <- matrix(0 , length(t_observation) , 1)  

for(i in 1:length(RRTimes)){  
X_Sim <- rep(0 , 18)
X_Sim[1] <- 0 # Baseline
X_Sim[2] <- t[i] #Rceb
X_Sim[3] <- 0.01 #Rwidth
X_Sim[4] <- 164 #RA
X_Sim[5] <- t[i] - 0.02 #Qcen
X_Sim[6] <- 0.01 #Qwidth
X_Sim[7] <- -23 #QA
X_Sim[8] <- t[i] + 0.02 #Scen
X_Sim[9] <- 0.012 #Swidth
X_Sim[10] <- -30 #SA
X_Sim[11] <- t[i] - 0.155 # PcenL
X_Sim[12] <- 0.01 # PWidthL

if(mod(i , 2) == 1){
X_Sim[13] <- 0 # PAL
}else{
X_Sim[13] <- 8
}
X_Sim[14] <- t[i] - 0.105 #PcenR
X_Sim[15] <- 0.02 #PenR

if(mod(i , 2) == 1){
X_Sim[16] <- 0 # PAR
}else{
  X_Sim[16] <- 12
}
X_Sim[17] <- t[i] + 0.23 #TcenL
X_Sim[18] <- 0.05 #TwidthL
X_Sim[19] <- 20 #TwidthL
X_Sim[20] <- t[i] + 0.24 #TcenR
X_Sim[21] <- 0.04 #TwidthR
X_Sim[22] <- 18 #TAR


output <- output + ECGSim_WrapperSingleBeat_m2(X_Sim , t_observation )
}
return(output)
}
PER_CreateECGPACS <- function( t , t_observation , RRTimes ){
  
  output <- matrix(0 , length(t_observation) , 1)  
  
  for(i in 1:dim(RRTimes)[1]){  
    X_Sim <- rep(0 , 18)
    X_Sim[1] <- 0 # Baseline
    X_Sim[2] <- t[i] #Rceb
    X_Sim[3] <- 0.01 #Rwidth
    X_Sim[4] <- 164 #RA
    X_Sim[5] <- t[i] - 0.02 #Qcen
    X_Sim[6] <- 0.01 #Qwidth
    X_Sim[7] <- -23 #QA
    X_Sim[8] <- t[i] + 0.02 #Scen
    X_Sim[9] <- 0.012 #Swidth
    X_Sim[10] <- -30 #SA
    X_Sim[11] <- t[i] - 0.155 # PcenL
    X_Sim[12] <- 0.01 # PWidthL
    if(RRTimes[i,2] == 2){
    X_Sim[13] <- 0  
    }else{
    X_Sim[13] <- 8}
    X_Sim[14] <- t[i] - 0.105 #PcenR
    X_Sim[15] <- 0.02 #PenR
    if(RRTimes[i,2] == 2){
      X_Sim[16] <- 0  
    }else{
      X_Sim[16] <- 12}
    X_Sim[17] <- t[i] + 0.23 #TcenL
    X_Sim[18] <- 0.05 #TwidthL
    X_Sim[19] <- 20 #TwidthL
    X_Sim[20] <- t[i] + 0.24 #TcenR
    X_Sim[21] <- 0.04 #TwidthR
    X_Sim[22] <- 18 #TAR
    
    
    output <- output + ECGSim_WrapperSingleBeat_m2(X_Sim , t_observation )
  }
  return(output)
}
PER_CreateECGAFib <- function( t , t_observation , RRTimes ){
  
  output <- matrix(0 , length(t_observation) , 1)  
  
  for(i in 1:length(RRTimes)){  
    X_Sim <- rep(0 , 18)
    X_Sim[1] <- 0 # Baseline
    X_Sim[2] <- t[i] #Rceb
    X_Sim[3] <- 0.01 #Rwidth
    X_Sim[4] <- 164 #RA
    X_Sim[5] <- t[i] - 0.02 #Qcen
    X_Sim[6] <- 0.01 #Qwidth
    X_Sim[7] <- -23 #QA
    X_Sim[8] <- t[i] + 0.02 #Scen
    X_Sim[9] <- 0.012 #Swidth
    X_Sim[10] <- -30 #SA
    X_Sim[11] <- t[i] - 0.155 # PcenL
    X_Sim[12] <- 0.01 # PWidthL
    X_Sim[13] <- 0  
    X_Sim[14] <- t[i] - 0.105 #PcenR
    X_Sim[15] <- 0.02 #PenR
    X_Sim[16] <- 0  
    X_Sim[17] <- t[i] + 0.23 #TcenL
    X_Sim[18] <- 0.05 #TwidthL
    X_Sim[19] <- 20 #TwidthL
    X_Sim[20] <- t[i] + 0.24 #TcenR
    X_Sim[21] <- 0.04 #TwidthR
    X_Sim[22] <- 18 #TAR
    output <- output + ECGSim_WrapperSingleBeat_m2(X_Sim , t_observation )
  }
  return(output)
}
PER_CreateECGReg <- function( t , t_observation , RRTimes ){
  
  output <- matrix(0 , length(t_observation) , 1)  

  X_Sim <- rep(0 , 22)
  for(i in 1:length(RRTimes)){  
    X_Sim[1] <- 0 # Baseline
    X_Sim[2] <- t[i] #Rceb
    X_Sim[3] <- 0.01 #Rwidth
    X_Sim[4] <- 164 #RA
    X_Sim[5] <- t[i] - 0.02 #Qcen
    X_Sim[6] <- 0.01 #Qwidth
    X_Sim[7] <- -23 #QA
    X_Sim[8] <- t[i] + 0.02 #Scen
    X_Sim[9] <- 0.012 #Swidth
    X_Sim[10] <- -30 #SA
    X_Sim[11] <- t[i] - 0.125 # PcenL
    X_Sim[12] <- 0.02 # PWidthL
    X_Sim[13] <- 6  
    X_Sim[14] <- t[i] - 0.125 #PcenR
    X_Sim[15] <- 0.02 #PenR
    X_Sim[16] <- 6  
    X_Sim[17] <- t[i] + 0.23 #TcenL
    X_Sim[18] <- 0.05 #TwidthL
    X_Sim[19] <- 20 #TwidthL
    X_Sim[20] <- t[i] + 0.24 #TcenR
    X_Sim[21] <- 0.04 #TwidthR
    X_Sim[22] <- 18 #TAR
    output <- output + ECGSim_WrapperSingleBeat_m2(X_Sim , t_observation )
  }
  return(output)
}
PER_CreatePwaveECGXReg <- function(t = 1 , X1 =6 , X2 =6 , X3 = 0.125 , X4 = 0.125 , X5 =0.02 , X6 =0.02){
  X_Sim <- rep(0 , 22)
  X_Sim[1] <- 0 # Baseline
  X_Sim[2] <- t #Rceb
  X_Sim[3] <- 0.01 #Rwidth
  X_Sim[4] <- 164 #RA
  X_Sim[5] <- t - 0.02 #Qcen
  X_Sim[6] <- 0.01 #Qwidth
  X_Sim[7] <- -23 #QA
  X_Sim[8] <- t + 0.02 #Scen
  X_Sim[9] <- 0.012 #Swidth
  X_Sim[10] <- -30 #SA
  X_Sim[11] <- t - X3 # PcenL
  X_Sim[12] <- X5 # PWidthL
  X_Sim[13] <- X1  #PAL
  X_Sim[14] <- t - X4 #PcenR
  X_Sim[15] <- X6 #PenR
  X_Sim[16] <- X2  #PAR
  X_Sim[17] <- t + 0.23 #TcenL
  X_Sim[18] <- 0.05 #TwidthL
  X_Sim[19] <- 20 #TwidthL
  X_Sim[20] <- t + 0.24 #TcenR
  X_Sim[21] <- 0.04 #TwidthR
  X_Sim[22] <- 18 #TAR
  return(X_Sim)
}

