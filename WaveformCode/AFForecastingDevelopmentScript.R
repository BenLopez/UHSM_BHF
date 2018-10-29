


AFF_SetInAFtozero <- function(Data, patientnames , PatIndex2017 ){
  for(i in 1:dim(Data)[1]){
    PatientRecord = DP_ExtractPatientRecordforIndex(PatIndex2017  , patientnames[[i]] )
    Data[i , Data[i,,12] > as.numeric(DP_StripTime(PatientRecord$ConfirmedFirstNewAF))  , 1:11] <- 0
  }
  return(Data)
}

AFF_ExtractMeanandCovarianceTimeSeries <- function(Data){
  mu = matrix(0 , dim(Data)[2] , dim(Data)[3] -1 )
  sigma =  array(0 , c(dim(Data)[2] , dim(Data)[3] -1, dim(Data)[3] -1) )
  
  for(i in 1:dim(Data)[2]){
   if(sum( Data[ , i , 1] != 0) > 3){
     mu[i , ] = apply(Data[ , i , 1:(dim(Data)[3] -1) ] , 2 , function(X){ mean(X[X != 0]) }) 
     sigma[i , , ] = cov( Data[Data[ , i , 1] != 0 , i , 1:(dim(Data)[3] -1)])
   }
  }
  return(setNames(list(mu , sigma) , c('mu' , 'Sigma') ) )
}


SecondOrderStructAF <- AFF_ExtractMeanandCovarianceTimeSeries( AFF_SetInAFtozero(DataBaseMaster$AFPatientsDatabase, DataBaseMaster$AFPatinetsNames , PatIndex2017))
SecondOrderStructNAF <- AFF_ExtractMeanandCovarianceTimeSeries(DataBaseMaster$NAFPatientsDatabase[sample(1:dim(DataBaseMaster$NAFPatientsDatabase)[1] ,  100 ) , , ])
SecondOrderStructNAF2 <- AFF_ExtractMeanandCovarianceTimeSeries(DataBaseMaster$NAFPatientsDatabase[sample(1:dim(DataBaseMaster$NAFPatientsDatabase)[1] ,  100 ) , , ])


i <- 3
plot(DP_NormaliseData(SecondOrderStructAF$mu[SecondOrderStructAF$mu[,i] !=0,i] -  SecondOrderStructNAF$mu[SecondOrderStructAF$mu[,i]!=0,i]) , type ='l' , col = rgb(0,0,1,alpha = 0.5) , ylab = 'Variable' , xlab='time')
abline(0,0)


DiscrepancyMatrix <- matrix( 0 , dim(SecondOrderStructAF$mu )[1] , 1)
DiscrepancyMatrix2 <- matrix( 0 , dim(SecondOrderStructAF$mu )[1] , 1)

for(i in 1: dim(SecondOrderStructAF$mu )[1]){
  DiscrepancyMatrix[i,] = mahalanobis(SecondOrderStructAF$mu[i,] - SecondOrderStructNAF$mu[i,] , 0 , SecondOrderStructAF$Sigma[i,,] + SecondOrderStructNAF$Sigma[i,,]) 
  DiscrepancyMatrix2[i,] = mahalanobis(SecondOrderStructNAF2$mu[i,] - SecondOrderStructNAF$mu[i,] , 0 , SecondOrderStructNAF2$Sigma[i,,] + SecondOrderStructNAF$Sigma[i,,]) 
  }

x11(20,14)
plot(rollmean(DiscrepancyMatrix[SecondOrderStructAF$mu[,1] !=0,] - DiscrepancyMatrix2[SecondOrderStructAF$mu[,1] !=0,] , k = 500) , type = 'l' , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'Time Realtive to FNAF' , ylab = 'Discrepancy Between AF and NAF populations' , ylim = c(0,100))
title('Time Series of Discrepancy between Population who go into AF and do not go into AF')
abline(0,0)

{patientnumeber = 21
DiscrepancyMatrix1 <- matrix( 0 , dim(SecondOrderStructAF$mu )[1] , 1)
DiscrepancyMatrix2 <- matrix( 0 , dim(SecondOrderStructAF$mu )[1] , 1)
for(i in 1: dim(SecondOrderStructAF$mu )[1]){
  DiscrepancyMatrix1[i,] = mahalanobis(DataBaseMaster$AFPatientsDatabase[patientnumeber , i , 1:11] , SecondOrderStructNAF$mu[i,]  , SecondOrderStructNAF$Sigma[patientnumeber,,]) 
  DiscrepancyMatrix2[i,] = mahalanobis(DataBaseMaster$NAFPatientsDatabase[patientnumeber , i , 1:11] , SecondOrderStructNAF$mu[i,]  , SecondOrderStructNAF$Sigma[patientnumeber,,]) 
}


DiscrepancyMatrix1[DataBaseMaster$AFPatientsDatabase[patientnumeber ,  , 1] == 0 , ] <- 0
DiscrepancyMatrix2[DataBaseMaster$NAFPatientsDatabase[patientnumeber ,  , 1] == 0 , ] <-0

RPeaks1 <- DP_LoadRpeaksfile(path , DataBaseMaster$AFPatinetsNames[[patientnumeber]])
RPeaks2 <- DP_LoadRpeaksfile(path , DataBaseMaster$NAFPatinetsNames[[patientnumeber]])


x11(20,14)
par(mfrow = c(2 , 1))
plot( RPeaks1$RRCombined$RR[1:min(which(SecondOrderStructAF$mu[,1] ==0))] , col = rgb(0,0,1,alpha = 0.05) , xlab= 'Time Realtive to FNAF' , ylab='RR times')
points( RPeaks2$RRCombined$RR[1:min(which(SecondOrderStructAF$mu[,1] ==0))] , col = rgb(1,0,0,alpha = 0.05) )

plot((DiscrepancyMatrix1[SecondOrderStructAF$mu[,1] !=0,]) , type = 'l' , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'Time Realtive to FNAF' , ylab = 'Discrepancy Between AF and NAF populations')
title('Time Series of Discrepancy between Population who go into AF and do not go into AF')
lines((DiscrepancyMatrix2[SecondOrderStructAF$mu[,1] !=0,]) , type = 'l' , col = rgb(1,0,0 , alpha = 0.5) , xlab = 'Time Realtive to FNAF' , ylab = 'Discrepancy Between AF and NAF populations')
}


