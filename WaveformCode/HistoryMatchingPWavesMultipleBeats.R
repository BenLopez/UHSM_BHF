

x11()
plot(Xstar, EmulatedQS[3,] ,type ='l' , col =rgb(0,0,1,alpha = 0.1) , ylim = c(-20,20) )
for(i in 1:dim(EmulatedQS)[1]){
  lines(Xstar, EmulatedQS[i,] ,type ='l' , col =rgb(0,0,1,alpha = 0.1))
}
abline(v = PWaveSummariesOutputs[ 3 , kk])
abline(v = PWaveSummariesOutputs[ 4 , kk])
lines(Xstar , mQS)
lines(Xstar , mQS + 2*sqrt(vQS) , col ='red')
lines(Xstar , mQS - 2*sqrt(vQS) , col ='red')



Implausability <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])

ModelDiscrepancyMatrix <- apply(PriorNonImplausibleSet , 1 , function(X){sqrt(ModelDiscrepancy(X , Xstar , PsimulatorFunction))} )

z = EmulatedQS

for(i in 1:dim(PriorNonImplausibleSet)[1]){
Betas = matrix(Hinvstruct[ , i] , 4  , 51)%*%t(z)
  
Implausability[i , ] <- colMeans2( abs(t(z) - matrix(Hstruct[ , i] , 51  , 4)%*%Betas) / ModelDiscrepancyMatrix[,i])
#apply(apply(abs(t(z) - H%*%Betas) , 2 ,  function(X){X / ModelDiscrepancyMatrix[,i]} ) , 2 , mean)
#DP_WaitBar(i / dim(PriorNonImplausibleSet)[1])
}


XminStruct <- PriorNonImplausibleSet[apply(Implausability , 2 , which.min) , ]
XminStruct <-XminStruct[apply(Implausability  , 2 , min) < 2, ]

x11()
plot(Xstar , PsimulatorFunction( XminStruct[1 , ] , Xstar) , type = 'l', col =rgb(0,0,1,alpha = 0.1))
for(i in 1:dim(EmulatedQS)[1]){
  lines(Xstar,  PsimulatorFunction( XminStruct[i , ] , Xstar) ,type ='l' , col =rgb(0,0,1,alpha = 0.1))
}

T_start <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff(PsimulatorFunction(X , Xstar)> 0.075))  )==1 )[1]]} )
T_end <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff(PsimulatorFunction(X , Xstar)> 0.075))  )==1 )[2]]} )
T_end[is.na(T_end)] <- 1

PWaveDurations <- as.numeric((T_end - T_start)*(median(diff(QS_Struct$t_start) , na.rm = T) - (20*0.005) ))*1000
PAmplitudes <- apply(XminStruct , 1 ,function(X){which.max( PsimulatorFunction(X , Xstar) ) } )
for(i in 1:length(PAmplitudes)){
  PAmplitudes[i] <- z[i , PAmplitudes[i]]
}
#PLocations <- apply(XminStruct , 1 ,function(X){which.max( PsimulatorFunction(X , Xstar) ) } )
#for(i in 1:length(PAmplitudes)){
#  PLocations[i] <- Xstar[PLocations[i]]
#}

PWaveDispertion <- as.numeric(quantile(PWaveDurations , 0.98) - quantile(PWaveDurations , 0.02))
PMaxDis <-  as.numeric(quantile(PWaveDurations , 0.98))
PMinDis <-  as.numeric(quantile(PWaveDurations , 0.02))
PEDies <- mean( PWaveDurations , na.rm = T)
PVarDis <-  var( PWaveDurations , na.rm = T) 
PTotal <- as.numeric((quantile(T_end , 0.99) - quantile(T_start , 0.01))*(median(diff(QS_Struct$t_start) , na.rm = T) - (20*0.005) )*1000)
PEAmp <- mean(PAmplitudes, na.rm = T)
PVAmp <- var(PAmplitudes , na.rm = T)

output <- PWaveHM_HistoryMatchGroupofPwaves(z , QS_Struct , PriorNonImplausibleSet  , ModelDiscrepancyMatrix, Hinvstruct , Hstruct)

PWaveHM_HistoryMatchGroupofPwaves <- function(z , QS_Struct , PriorNonImplausibleSet , ModelDiscrepancyMatrix, Hinvstruct , Hstruct){
# Precalculations
  
if(!exists('Hinvstruct')){
Hinvstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
    H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)
    return(solve(t(H)%*%H)%*%t(H))
  })}
if(!exists('Hinvstruct')){
Hstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
  return(  H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction))
})  
}
if(!exists('ModelDiscrepancyMatrix')){
  ModelDiscrepancyMatrix <- apply(PriorNonImplausibleSet , 1 , function(X){sqrt(ModelDiscrepancy(X , Xstar , PsimulatorFunction))} )
}

# History Match    
  Implausability <- matrix(0 , dim(PriorNonImplausibleSet)[1] , dim(z)[1])
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    Betas = matrix(Hinvstruct[ , i] , 4  , 51)%*%t(z)
    
    Implausability[i , ] <- colMeans2( abs(t(z) - matrix(Hstruct[ , i] , 51  , 4)%*%Betas) / ModelDiscrepancyMatrix[,i])
    #apply(apply(abs(t(z) - H%*%Betas) , 2 ,  function(X){X / ModelDiscrepancyMatrix[,i]} ) , 2 , mean)
    #DP_WaitBar(i / dim(PriorNonImplausibleSet)[1])
  }
  
# Extract Useful Data

  XminStruct <- PriorNonImplausibleSet[apply(Implausability , 2 , which.min) , ]
  XminStruct <-XminStruct[apply(Implausability  , 2 , min) < 2, ]
  
  T_start <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff(PsimulatorFunction(X , Xstar)> 0.075))  )==1 )[1]]} )
  T_end <- apply(XminStruct , 1 , function(X){Xstar[which(abs(c(0, diff(PsimulatorFunction(X , Xstar)> 0.075))  )==1 )[2]]} )
  T_end[is.na(T_end)] <- 1
  T_start[is.na(T_start)] <- mean(T_start , na.rm=T)
  PWaveDurations <- as.numeric((T_end - T_start)*(median(diff(QS_Struct$t_start) , na.rm = T) - (20*0.005) ))*1000
  PAmplitudes <- apply(XminStruct , 1 ,function(X){which.max( PsimulatorFunction(X , Xstar) ) } )
  for(i in 1:length(PAmplitudes)){
    PAmplitudes[i] <- z[i , PAmplitudes[i]]
  }
  PLocations <- apply(XminStruct , 1 ,function(X){which.max( PsimulatorFunction(X , Xstar) ) } )
  for(i in 1:length(PAmplitudes)){
    PLocations[i] <- Xstar[PLocations[i]]
  }

  PWaveDurations <-PWaveDurations[!is.na(PWaveDurations)]
    
  output <- rep(0 , 8)
  output[1] <-  as.numeric(quantile(PWaveDurations , 0.98) - quantile(PWaveDurations , 0.02)) #dispersion
  output[2] <-  as.numeric(quantile(PWaveDurations , 0.98)) # max
  output[3] <-  as.numeric(quantile(PWaveDurations , 0.02)) # min
  output[4] <-  mean( PWaveDurations , na.rm = T) # mean duration
  output[5] <-  var( PWaveDurations , na.rm = T) # var duration
  output[6] <-  as.numeric((quantile(T_end , 0.99) - quantile(T_start , 0.01))*(median(diff(QS_Struct$t_start) , na.rm = T) - (20*0.005) )*1000) #total
  output[7] <-  mean(PAmplitudes, na.rm = T) # mean amp
  output[8] <-  var(PAmplitudes , na.rm = T) # var amp
  
  return(output)  
}
