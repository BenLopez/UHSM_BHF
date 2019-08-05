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


# Iteration one. 

{numberofsamples <- 1000000
PriorNonImplausibleSetBigeminy <- matrix(0 , numberofsamples , 10)

PriorNonImplausibleSetBigeminy[ , 1] <- 0.5
PriorNonImplausibleSetBigeminy[ , 2] <- 0.5
PriorNonImplausibleSetBigeminy[ , 3] <- 0

PriorNonImplausibleSetBigeminy[ , 4] <- runif(numberofsamples,0.3 , 2)
PriorNonImplausibleSetBigeminy[ , 5] <- runif(numberofsamples,0.3 , 2)
PriorNonImplausibleSetBigeminy[ , 6] <- runif(numberofsamples,0.3 , 2)

PriorNonImplausibleSetBigeminy[ , 7] <- runif(numberofsamples,0.000001 , 0.08)
PriorNonImplausibleSetBigeminy[ , 8] <- runif(numberofsamples,0.000001 , 0.08)
PriorNonImplausibleSetBigeminy[ , 9] <- runif(numberofsamples,0.000001 , 0.08)

PriorNonImplausibleSetBigeminy[ , 10] <- 1


Mu12_UB = function(x){  return( ((1 - 1.7)/(1 - 0.3))*x + 2 ) } 
Mu12_LB = function(x){  return( ((0.1 - 0.03)/(1 - 0.3))*x ) }}
 
  
ValidVector <- (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 5]) & 
  ((PriorNonImplausibleSetBigeminy[ , 4] + PriorNonImplausibleSetBigeminy[ , 5]) < 2) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) > 0.1) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) > Mu1_LB(PriorNonImplausibleSetBigeminy[ , 4]) ) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) < Mu1_UB(PriorNonImplausibleSetBigeminy[ , 4]) )

x11()
plot(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:10000] , 4:5] , col = rgb(0,0,1 , alpha = 0.1) , pch = 16 , xlim = c(0,2) , ylim = c(0,2))


ValidVector <- (PriorNonImplausibleSetBigeminy[ , 4] < PriorNonImplausibleSetBigeminy[ , 5]) & 
  ((PriorNonImplausibleSetBigeminy[ , 4] + PriorNonImplausibleSetBigeminy[ , 5]) < 2) &
  ((PriorNonImplausibleSetBigeminy[ , 4] + PriorNonImplausibleSetBigeminy[ , 5]) > 0.1) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) > Mu1_LB(PriorNonImplausibleSetBigeminy[ , 4]) ) &
  (abs(PriorNonImplausibleSetBigeminy[ , 4] - PriorNonImplausibleSetBigeminy[ , 5]) < Mu1_UB(PriorNonImplausibleSetBigeminy[ , 4]) ) &  
  (abs(PriorNonImplausibleSetBigeminy[ , 7] + PriorNonImplausibleSetBigeminy[ , 8]) < 0.08)

Sigma12_UB = function(x){  return( ((1 - 1.7)/(1 - 0.3))*x + 2 ) } 
Sigma12_LB = function(x){  return( ((0.1 - 0.03)/(1 - 0.3))*x ) }
  


BC_PlotPairs(PriorNonImplausibleSetBigeminy[which(ValidVector)[1:10000] , c(4 , 5 ,7 , 8)] )
