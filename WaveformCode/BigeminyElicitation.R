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

# Some functions
{
FM_CalulateDistance <- function(X , ValidVector ){
    m <- apply(X[ValidVector,] , 2 , mean)
    CovarianceMatrix = cov(X[ValidVector,])
    return( return(mahalanobis(X , m , CovarianceMatrix)) ) 
  }
  
FM_FindWhyNotValid <- function(X){
    X <- t(as.matrix(X))
    output <- matrix(0  , 8 , 1)
    
    if(X[ , 4] < X[ , 6]){
      output[1,] <- 1
    } 
    if((X[ , 6]) < 2){
      output[2,] <- 1  
    }
    if((X[ , 4]) > 0.2){ 
      output[3,] <- 1
    }
    if(abs(X[ , 4] - X[ , 6]) > 0.2){
      output[4,] <- 1
    }
    if(abs(X[ , 4] - X[ , 6]) > Mu12_LB(X[ , 4]) ){
      output[5,] <- 1
      } 
    if(abs(X[ , 4] - X[ , 6]) < Mu12_UB(X[ , 4]) ){
      output[6,] <- 1
    }
    if(abs(X[ , 7] + X[ , 9]) < 0.16){
      output[7,] <- 1
    }
    if(X[ , 7] < Sigma12_UB(abs(X[ , 4] - X[ , 6]) )){
      output[8,] <- 1
    }
    return(output)
    }
  
}

##### Iteration one ##### 

numberofsamples <- 1000000
PriorNonImplausibleSetBigeminy <- matrix(0 , numberofsamples , 10)

PriorNonImplausibleSetBigeminy[ , 1] <- 0.5
PriorNonImplausibleSetBigeminy[ , 2] <- 0
PriorNonImplausibleSetBigeminy[ , 3] <- 0.5

PriorNonImplausibleSetBigeminy[ , 4] <- runif(numberofsamples,0.2 , 2)
PriorNonImplausibleSetBigeminy[ , 5] <- 0.00000001
PriorNonImplausibleSetBigeminy[ , 6] <- runif(numberofsamples,0.2 , 2)

PriorNonImplausibleSetBigeminy[ , 7] <- runif(numberofsamples,0.000001 , 0.16)
PriorNonImplausibleSetBigeminy[ , 8] <- 0.00000001
PriorNonImplausibleSetBigeminy[ , 9] <- runif(numberofsamples,0.000001 , 0.16)

PriorNonImplausibleSetBigeminy[ , 10] <- 1

Mu12_UB = function(x){  return( ((1 - 1.7)/(1 - 0.3))*x + 2 ) } 
Mu12_LB = function(x){  return( ((0.1 - 0.03)/(1 - 0.3))*x ) }
 
  
ValidVector <- (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 6]) & 
  ((PriorNonImplausibleSetBigeminy[ , 6]) < 2) &
  ((PriorNonImplausibleSetBigeminy[ , 4]) > 0.2) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 6]) > 0.2) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 6]) > Mu12_LB(PriorNonImplausibleSetBigeminy[ , 4]) ) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 6]) < Mu12_UB(PriorNonImplausibleSetBigeminy[ , 4]) )

Sigma12_UB = function(x){  return( ((-0.001 + 0.09)/(-0.03 + 1.7))*x +   -0.000598802395209581 ) } 

ValidVector <- (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 6]) & 
  ((PriorNonImplausibleSetBigeminy[ , 6]) < 2) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 6]) > 0.2) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 6]) > Mu12_LB(PriorNonImplausibleSetBigeminy[ , 4]) ) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 6]) < Mu12_UB(PriorNonImplausibleSetBigeminy[ , 4]) ) &  
  (abs(PriorNonImplausibleSetBigeminy[ , 7] + PriorNonImplausibleSetBigeminy[ , 9]) < 0.16)&
  ((PriorNonImplausibleSetBigeminy[ , 7] + PriorNonImplausibleSetBigeminy[ , 9]) < Sigma12_UB(abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 6]) ))

BC_PlotPairs(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000], c(4 , 7 ,6 , 9)] , alpha = 0.1)

BC_PlotPairs(cbind(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000], c(4)],
                   PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000], c(6 )]-PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000], c(4)],
             PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000], c( 7)],
             PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000], c( 7)] + PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000], c( 9)]), alpha = 0.1)



TestSet <- list()
for(i in 1:5){
TestSet[[i]] <- HREL_BigeminySampleECG(X=PriorNonImplausibleSetBigeminy[which(ValidVector)[i],] )
}

DistanceVector <- FM_CalulateDistance(X = PriorNonImplausibleSetBigeminy[ , c( 4, 6, 7 , 9)] , ValidVector = ValidVector )

Threshold1 <- quantile(DistanceVector[ValidVector] , 0.98)
Threshold2 <- quantile(DistanceVector[ValidVector] , 0.99)


EdgeValid <- which( (ValidVector)&(DistanceVector>Threshold1)&(DistanceVector<Threshold2))
EdgeNotValid <- which( (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 6]) &
                        (!ValidVector)&
                         (DistanceVector>Threshold1)&
                        (DistanceVector<Threshold2) )




BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:1000] , c(4 , 7 ,6 , 9)] ,
                               PriorNonImplausibleSetBigeminy[EdgeValid[1:1000] , c(4 , 7 ,6 , 9)] ,
                               PriorNonImplausibleSetBigeminy[EdgeNotValid[1:1000] , c(4 , 7 ,6 , 9)] )

 
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

GroundTruth <- PermutedList
GroundTruth[GroundTruth <= 10] <- T
GroundTruth[GroundTruth > 10] <- F


TestPoints <- rbind(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:5], ],
PriorNonImplausibleSetBigeminy[EdgeValid[1:5], ],
PriorNonImplausibleSetBigeminy[EdgeNotValid[1:5],],
PriorNonImplausibleSetRegular[(dim(PriorNonImplausibleSetRegular)[1] - 20):(dim(PriorNonImplausibleSetRegular)[1] - 16),])
NumericOutputs <- list(PermutedList , GroundTruth,TestPoints,TestSet)

#save(NumericOutputs , file = 'C:\\Users\\Ben\\Documents\\HeartRhythm Elicitation\\Bigeminy\\IterationTwo\\NumericOutputs.RData' )

x11(20,15)
grid.arrange( TestSet[[PermutedList[1]]] , TestSet[[PermutedList[2]]] , TestSet[[PermutedList[3]]] , TestSet[[PermutedList[4]]] , TestSet[[PermutedList[5]]] , nrow = 5)

x11(20,15)
grid.arrange( TestSet[[PermutedList[6]]] , TestSet[[PermutedList[7]]] , TestSet[[PermutedList[8]]] , TestSet[[PermutedList[9]]] , TestSet[[PermutedList[10]]] , nrow = 5)

x11(20,15)
grid.arrange( TestSet[[PermutedList[11]]] , TestSet[[PermutedList[12]]] , TestSet[[PermutedList[13]]] , TestSet[[PermutedList[14]]] , TestSet[[PermutedList[15]]] , nrow = 5)

x11(20,15)
grid.arrange( TestSet[[PermutedList[16]]] , TestSet[[PermutedList[17]]] , TestSet[[PermutedList[18]]] , TestSet[[PermutedList[19]]] , TestSet[[PermutedList[20]]] , nrow = 5)

##### End Elicitation One #####

xtable(read.csv(file = "C:\\Users\\Ben\\Documents\\HeartRhythm Elicitation\\Bigeminy\\IterationOne\\BigeminyItr1Form(withanswers).csv"))
