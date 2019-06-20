PER <- CreateDefaultX <- function(mu1 = 0.6 , mu2= 0.9 , sigma1=0.01 , sigma2=0.01){
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