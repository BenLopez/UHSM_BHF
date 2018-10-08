{
  if(file.exists('CheckforDefaultsScript.R')){
    source('CheckforDefaultsScript.R')
  }else{
    pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
    source("LibrariesAndSettings.R" , print.eval  = TRUE )
    DP_LoadPatientIndex()
    DP_ChooseDataReps()
    FilestoProcess <- DP_ChooseECGstoProcess() 
    HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
  }
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

beta = c(0.3 , 0.7)

numbertrainingpoints = 100000


SamplePoints <- function(alpha){

MclustDistributionStruct <- BC_CreateDefaultmclustStruct()

MclustDistributionStruct$parameters$pro <- c(alpha , 1-2*alpha ,alpha)
MclustDistributionStruct$parameters$mean <- matrix(0, 1 , 3)
MclustDistributionStruct$parameters$mean[ , 1] <- 5
MclustDistributionStruct$parameters$mean[ , 2] <- 10
MclustDistributionStruct$parameters$mean[ , 3] <- 15

MclustDistributionStruct$parameters$variance$sigma <- array(0 , c(1 , 1 , 3))
MclustDistributionStruct$parameters$variance$sigma[ , , 1] <- c( 2.5) 
MclustDistributionStruct$parameters$variance$sigma[ , , 2] <- c( 5) 
MclustDistributionStruct$parameters$variance$sigma[ , , 3] <- c( 2.5)
SampleofPoints <<- BC_SampleGMM(MclustDistributionStruct = MclustDistributionStruct , numberofsamples =    beta[1]*numbertrainingpoints)

MclustDistributionStruct2 <- BC_CreateDefaultmclustStruct()

MclustDistributionStruct2$parameters$pro <- c(alpha , 1-2*alpha ,alpha)
MclustDistributionStruct2$parameters$mean <- matrix(0, 1 , 3)
MclustDistributionStruct2$parameters$mean[ , 1] <- 0
MclustDistributionStruct2$parameters$mean[ , 2] <- 5
MclustDistributionStruct2$parameters$mean[ , 3] <- 10

MclustDistributionStruct2$parameters$variance$sigma <- array(0 , c(1 , 1 , 3))
MclustDistributionStruct2$parameters$variance$sigma[ , , 1] <- c( 2.5) 
MclustDistributionStruct2$parameters$variance$sigma[ , , 2] <- c( 5) 
MclustDistributionStruct2$parameters$variance$sigma[ , , 3] <- c( 2.5)

SampleofPoints2 <<- BC_SampleGMM(MclustDistributionStruct = MclustDistributionStruct2 , numberofsamples =  beta[2]*numbertrainingpoints)

mu <<- mean(SampleofPoints)
Sigma <<- var(SampleofPoints)
mu1 <<- mean(SampleofPoints2)
Sigma1 <<- var(SampleofPoints2)

}


listofprops <- seq(0 , 0.25 , 0.01)


PPV1 <- matrix(0 , length(listofprops) , 1)
PPV2 <- matrix(0 , length(listofprops) , 1)


for( ii in 1:length(listofprops) ){

SamplePoints(listofprops[ii])

Z = rbind(as.matrix(SampleofPoints) , as.matrix(SampleofPoints2))
Im1 = mahalanobis(Z , mu , Sigma)
Im2 = mahalanobis(Z , mu1 , Sigma1)

ImThresh1 <- mean(Im1) + 3*sqrt(var(Im1))
ImThresh2 <- mean(Im2) + 3*sqrt(var(Im2))

LogicalVector = matrix(0 , dim(Z)[1] , 1)
LogicalVector[(Im1 < ImThresh1 )*(Im2 > ImThresh2) == 1] = 1
LogicalVector[(Im1 > ImThresh1 )*(Im2 < ImThresh2) == 1] = -1
LogicalVector[(Im1 < ImThresh1 )*(Im2 < ImThresh2) == 1] = 0

LogicalVector <- BC_CreateProbCalibrationStruct(LogicalVector , beta , numbertrainingpoints)


Sensivity1 <- sum((LogicalVector[ , 1] == 1)*(LogicalVector[ , 2] == 1)) / sum(LogicalVector[ , 2])
Specifictity1 <- sum((LogicalVector[ , 1] != 1)*(LogicalVector[ , 2] != 1)) / sum(LogicalVector[ , 2] !=1)

Sensivity2 <- sum((LogicalVector[ , 1] == -1)*(LogicalVector[ , 2] == 0)) / sum(LogicalVector[ , 2] !=1 )
Specifictity2 <- sum((LogicalVector[ , 1] != -1)*(LogicalVector[ , 2] == 1)) / sum(LogicalVector[ , 2] ==1)

PPV1[ii] <- beta[1]*Sensivity1/(beta[1]*Sensivity1 + beta[2]*(1-Specifictity1) )
PPV2[ii] <- beta[2]*Sensivity2/(beta[2]*Sensivity2 + beta[1]*(1-Specifictity2) )
}  


bayesPPV1 <- matrix(0 , length(listofprops) , 1)
bayesPPV2 <- matrix(0 , length(listofprops) , 1)


for( ii in 1:length(listofprops) ){

SamplePoints(listofprops[ii])
Z = rbind(as.matrix(SampleofPoints) , as.matrix(SampleofPoints2))
  
P = (beta[1]*dnorm(Z , mu , sqrt(Sigma))) / (beta[1]*dnorm(Z , mu , sqrt(Sigma)) + beta[2]*dnorm(Z , mu1 , sqrt(Sigma1)) )

LogicalVector = matrix(0 , dim(Z)[1] , 1)
LogicalVector[P > 0.95] = 1
LogicalVector[P < 0.05] = -1

LogicalVector <- BC_CreateProbCalibrationStruct(LogicalVector , beta , numbertrainingpoints)

Sensivity1 <- sum((LogicalVector[ , 1] == 1)*(LogicalVector[ , 2] == 1)) / sum(LogicalVector[ , 2])
Specifictity1 <- sum((LogicalVector[ , 1] != 1)*(LogicalVector[ , 2] != 1)) / sum(LogicalVector[ , 2] !=1)

Sensivity2 <- sum((LogicalVector[ , 1] == -1)*(LogicalVector[ , 2] == 0)) / sum(LogicalVector[ , 2] !=1 )
Specifictity2 <- sum((LogicalVector[ , 1] != -1)*(LogicalVector[ , 2] == 1)) / sum(LogicalVector[ , 2] ==1)

bayesPPV1[ii] <- beta[1]*Sensivity1/(beta[1]*Sensivity1 + beta[2]*(1-Specifictity1) )
bayesPPV2[ii] <- beta[2]*Sensivity2/(beta[2]*Sensivity2 + beta[1]*(1-Specifictity2) )

}

plot1 <- ggplot()+
         geom_point(data = data.frame(x = listofprops , y=PPV1)  , aes(x,y) , col = 'blue') +
         geom_point(data = data.frame(x = listofprops , y=PPV2) , aes(x,y) , col = 'red') +
         geom_point(data = data.frame(x = listofprops , y=bayesPPV1)  , aes(x,y) , col = 'black') +
         geom_point(data = data.frame(x = listofprops , y=bayesPPV2) , aes(x,y) , col = 'green')
