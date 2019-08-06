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
UseAnnotatedData = 0
source('FM_CreateRhythumPriors.R')

##### Iteration one ##### 

numberofsamples <- 1000000
PriorNonImplausibleSetBigeminy <- matrix(0 , numberofsamples , 10)

PriorNonImplausibleSetBigeminy[ , 1] <- 0.5
PriorNonImplausibleSetBigeminy[ , 2] <- 0.5
PriorNonImplausibleSetBigeminy[ , 3] <- 0

PriorNonImplausibleSetBigeminy[ , 4] <- runif(numberofsamples,0.2 , 2)
PriorNonImplausibleSetBigeminy[ , 5] <- runif(numberofsamples,0.2 , 2)
PriorNonImplausibleSetBigeminy[ , 6] <- runif(numberofsamples,0.2 , 2)

PriorNonImplausibleSetBigeminy[ , 7] <- runif(numberofsamples,0.000001 , 0.08)
PriorNonImplausibleSetBigeminy[ , 8] <- runif(numberofsamples,0.000001 , 0.08)
PriorNonImplausibleSetBigeminy[ , 9] <- runif(numberofsamples,0.000001 , 0.08)

PriorNonImplausibleSetBigeminy[ , 10] <- 1


Mu12_UB = function(x){  return( ((1 - 1.7)/(1 - 0.3))*x + 2 ) } 
Mu12_LB = function(x){  return( ((0.1 - 0.03)/(1 - 0.3))*x ) }
 
  
ValidVector <- (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 5]) & 
  ((PriorNonImplausibleSetBigeminy[ , 5]) < 2) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) > 0.2) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) > Mu12_LB(PriorNonImplausibleSetBigeminy[ , 4]) ) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) < Mu12_UB(PriorNonImplausibleSetBigeminy[ , 4]) )

x11()
plot(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000] , 4] , PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000] , 5] , col = rgb(0,0,1 , alpha = 0.1) , pch = 16 , xlim = c(0,2) , ylim = c(0,2))



Sigma12_UB = function(x){  return( ((0.08)/(0.9))*x + -0.008888888888888 ) } 

ValidVector <- (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 5]) & 
  ((PriorNonImplausibleSetBigeminy[ , 4] + PriorNonImplausibleSetBigeminy[ , 5]) < 2) &
  ((PriorNonImplausibleSetBigeminy[ , 4] + PriorNonImplausibleSetBigeminy[ , 5]) > 0.1) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) > Mu12_LB(PriorNonImplausibleSetBigeminy[ , 4]) ) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) < Mu12_UB(PriorNonImplausibleSetBigeminy[ , 4]) ) &  
  (abs(PriorNonImplausibleSetBigeminy[ , 7] + PriorNonImplausibleSetBigeminy[ , 8]) < 0.08)&
  (PriorNonImplausibleSetBigeminy[ , 7] < Sigma12_UB(abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) ))

BC_PlotPairs(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:10000] , c(4 , 5 ,7 , 8)] )

HREL_BigeminySampleECG <- function( X ){
  
RRTimes <-  FM_SampleGMMBigeminy(X , 100 )
t_observation = seq(0.25  , 10 , 0.005)
t = cumsum(RRTimes)
RRTimes[RRTimes<0.3] <- 0.3
RRTimes[RRTimes>2] <- 2
ECG <- PER_CreateECGReg( t , t_observation , RRTimes )
p3 <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
p3 <- p3 + theme(
  panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                  size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 1, linetype = 'solid',
                                  colour = rgb(1,0,0,0.25)), 
  panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                  colour = rgb(1,0,0,0.25))
) 

p3 <-p3 + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
return(p3)
}
HREL_SampleECG <- function( X ){
  
  RRTimes <-  FM_SampleGMM(X , 100 )
  t_observation = seq(0.25  , 10 , 0.005)
  t = cumsum(RRTimes)
  RRTimes[RRTimes<0.3] <- 0.3
  RRTimes[RRTimes>2] <- 2
  ECG <- PER_CreateECGReg( t , t_observation , RRTimes )
  p3 <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
  p3 <- p3 + theme(
    panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 1, linetype = 'solid',
                                    colour = rgb(1,0,0,0.25)), 
    panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                    colour = rgb(1,0,0,0.25))
  ) 
  
  p3 <-p3 + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
  return(p3)
}


TestSet <- list()
for(i in 1:5){
TestSet[[i]] <- HREL_BigeminySampleECG(X=PriorNonImplausibleSetBigeminy[which(ValidVector)[i],] )
}

FM_CalulateDistance <- function(X , ValidVector ){
  X <- apply( X , 2 , DP_NormaliseData)
  m <- apply(X[ValidVector,] , 2 , mean)
  return( apply(X , 1 , function(x){sum((x -m)^2)})) 
}

DistanceVector <- FM_CalulateDistance(X = PriorNonImplausibleSetBigeminy[ , c( 4, 5, 7 , 8)] , ValidVector = ValidVector )

Threshold1 <- quantile(DistanceVector[ValidVector] , 0.95)
Threshold2 <- quantile(DistanceVector[ValidVector] , 0.99)


EdgeValid <- which( (ValidVector)&(DistanceVector>Threshold1)&(DistanceVector<Threshold2))
EdgeNotValid <- which(  (PriorNonImplausibleSetBigeminy[ , 4] > 0.4) &  
                          (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 5]) & 
                          ((PriorNonImplausibleSetBigeminy[ , 4] + PriorNonImplausibleSetBigeminy[ , 5]) < 2)
                        &(!ValidVector)
                        &(DistanceVector<Threshold2))


BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000] , c(4 , 5 ,7 , 8)] ,
                               PriorNonImplausibleSetBigeminy[EdgeValid[1:1000] , c(4 , 5 ,7 , 8)] ,
                               PriorNonImplausibleSetBigeminy[EdgeNotValid[1:1000] , c(4 , 5 ,7 , 8)] )

 
for(i in c(6:10) ){
  TestSet[[i]] <- HREL_BigeminySampleECG(X = PriorNonImplausibleSetBigeminy[EdgeValid[i],] )
}

for(i in c(11:15)){
  TestSet[[i]] <- HREL_BigeminySampleECG(PriorNonImplausibleSetBigeminy[EdgeNotValid[i],] )
}

for(i in c(16:20)){
  TestSet[[i]] <- HREL_SampleECG(PriorNonImplausibleSetRegular[dim(PriorNonImplausibleSetRegular)[1] - i,] )
}



PermutedList <- sample(1:20 , 20)

#PriorNonImplausibleSetBigeminy[which(ValidVector)[1:5] , c(4 , 5 ,7 , 8)]
#PriorNonImplausibleSetBigeminy[EdgeValid[1:5] , c(4 , 5 ,7 , 8)]
#PriorNonImplausibleSetBigeminy[EdgeNotValid[1:5], c(4 , 5 ,7 , 8)]

x11(20,15)
grid.arrange(TestSet[[PermutedList[1]]] , TestSet[[PermutedList[2]]] , TestSet[[PermutedList[3]]] , TestSet[[PermutedList[4]]] , TestSet[[PermutedList[5]]] , nrow = 5)

x11(20,15)
grid.arrange(TestSet[[PermutedList[6]]] , TestSet[[PermutedList[7]]] , TestSet[[PermutedList[8]]] , TestSet[[PermutedList[9]]] , TestSet[[PermutedList[10]]] , nrow = 5)

x11(20,15)
grid.arrange(TestSet[[PermutedList[11]]] , TestSet[[PermutedList[12]]] , TestSet[[PermutedList[13]]] , TestSet[[PermutedList[14]]] , TestSet[[PermutedList[15]]] , nrow = 5)

x11(20,15)
grid.arrange(TestSet[[PermutedList[16]]] , TestSet[[PermutedList[17]]] , TestSet[[PermutedList[18]]] , TestSet[[PermutedList[19]]] , TestSet[[PermutedList[20]]] , nrow = 5)

##### End Elicitation One #####