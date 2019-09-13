{
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
  precomputedfolderpath <- DP_SelectPrecomputedFolder()
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

UseAnnotatedData = 1
source('FM_CreateRhythumPriors.R')
}

HREL_NSRRuleFailure <- function(X){
  
  output <- matrix(0,3,1)
  
  if(!(X[1] > 0.6) ){output[1]<- 1}
  if(!(X[1] < 1) ){output[2]<- 1}
  if(!(X[2] < (0.04/3)) ){{output[3]<- 1}}
  return(output)
}

##### Iteration one ##### 

numberofsamples <- 1000000
PriorNonImplausibleSetNSR <- BE_SampleLHSinab(a = c(0.5,0.0000001,-10,1.8) , b = c(1.1,(0.05/3),10,40) , numberofsamples)
ValidVector <- (PriorNonImplausibleSetNSR[,3] > -sqrt((PriorNonImplausibleSetNSR[,4] -0.9) - 0.01)) & 
  (PriorNonImplausibleSetNSR[,3] < sqrt((PriorNonImplausibleSetNSR[,4] -0.9) - 0.01))&
  (PriorNonImplausibleSetNSR[,1] > 0.6 & PriorNonImplausibleSetNSR[,1] < 1  ) &
  (PriorNonImplausibleSetNSR[,2] < (0.04/3))
  
  
BC_PlotPairs(PriorNonImplausibleSetNSR[which(ValidVector)[1:1000],] , alpha = 0.1 , labels = c('Mean' , 'Variance' , 'Skewness' , 'Kurtosis'))

TestSet <- list()
for(i in 1:5){
  TestSet[[i]] <- HREL_RegularSampleECG(X=PriorNonImplausibleSetNSR[which(ValidVector)[i],] )
}

DistanceVector <- FM_CalulateDistance(X = PriorNonImplausibleSetNSR , ValidVector = ValidVector )

Threshold1 <- quantile(DistanceVector[ValidVector] , 0.98)
Threshold2 <- quantile(DistanceVector[ValidVector] , 0.99)


EdgeValid <- which( (ValidVector)&(DistanceVector>Threshold1)&(DistanceVector<Threshold2))
EdgeNotValid <- which(   (PriorNonImplausibleSetNSR[,3] > -sqrt((PriorNonImplausibleSetNSR[,4] -0.9) - 0.01)) & 
                          (PriorNonImplausibleSetNSR[,3] < sqrt((PriorNonImplausibleSetNSR[,4] -0.9) - 0.01))&
                         (!ValidVector)&
                         (DistanceVector>Threshold1)&
                         (DistanceVector<Threshold2) )




BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetNSR[which(ValidVector)[1:1000], ] ,
                               PriorNonImplausibleSetNSR[EdgeValid[1:1000], ] ,
                               PriorNonImplausibleSetNSR[EdgeNotValid[1:1000], ] )


for(i in c(6:10) ){
  TestSet[[i]] <- HREL_RegularSampleECG(X = PriorNonImplausibleSetNSR[EdgeValid[i],] )
}

for(i in c(11:15)){
  TestSet[[i]] <- HREL_RegularSampleECG(PriorNonImplausibleSetNSR[EdgeNotValid[i],] )
}

for(i in c(16:20)){
  TestSet[[i]] <- HREL_SampleECG(PriorNonImplausibleSetRegular[dim(PriorNonImplausibleSetRegular)[1] - i,] )
}
for(i in c(21:25)){
  TestSet[[i]] <- HREL_SampleECG(PriorNonImplausibleSetRegularyIreRegular[dim(PriorNonImplausibleSetRegularyIreRegular)[1] - i,] )
}

PermutedList <- sample(1:25 , 25)

GroundTruth <- PermutedList

TestPoints <- rbind(PriorNonImplausibleSetNSR[which(ValidVector)[1:5], ],
                    PriorNonImplausibleSetNSR[EdgeValid[1:5], ],
                    PriorNonImplausibleSetNSR[EdgeNotValid[1:5],])
                    NumericOutputs <- list(PermutedList , GroundTruth,TestPoints,TestSet)


#save(NumericOutputs , file = 'C:\\Users\\Ben\\Documents\\HeartRhythm Elicitation\\NSR\\IterationOne\\NumericOutputs.RData' )

x11(20,15)
grid.arrange( TestSet[[PermutedList[1]]] , TestSet[[PermutedList[2]]] , TestSet[[PermutedList[3]]] , TestSet[[PermutedList[4]]] , TestSet[[PermutedList[5]]] , nrow = 5)

x11(20,15)
grid.arrange( TestSet[[PermutedList[6]]] , TestSet[[PermutedList[7]]] , TestSet[[PermutedList[8]]] , TestSet[[PermutedList[9]]] , TestSet[[PermutedList[10]]] , nrow = 5)

x11(20,15)
grid.arrange( TestSet[[PermutedList[11]]] , TestSet[[PermutedList[12]]] , TestSet[[PermutedList[13]]] , TestSet[[PermutedList[14]]] , TestSet[[PermutedList[15]]] , nrow = 5)

x11(20,15)
grid.arrange( TestSet[[PermutedList[16]]] , TestSet[[PermutedList[17]]] , TestSet[[PermutedList[18]]] , TestSet[[PermutedList[19]]] , TestSet[[PermutedList[20]]] , nrow = 5)

x11(20,15)
grid.arrange( TestSet[[PermutedList[21]]] , TestSet[[PermutedList[22]]] , TestSet[[PermutedList[23]]] , TestSet[[PermutedList[24]]] , TestSet[[PermutedList[25]]] , nrow = 5)





##### Iteration Two #####