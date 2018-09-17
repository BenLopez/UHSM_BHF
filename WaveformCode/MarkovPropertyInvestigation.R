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

{BCOptions <- setNames( list(1 , 1 , 1 , 1 , 1) , c('DataType' ,
                                                    'ClassifierType' , 
                                                    'DensityEstimationGlobal', 
                                                    'DensityEstimationlocal' , 
                                                    'GlobalUpdate') )
  BCOptions[[1]] <- 'DistributionSummaries'
  BCOptions[[2]] <- 'AFClassifier' 
  BCOptions[[3]] <- 'MVN' 
  BCOptions[[4]] <- 'GMM'
  BCOptions[[5]] <- 'No'
  
  BCParameters <- setNames( list(1,1,1,1) , c( 'TS_Likelihood_clique' , 'NumberComponentsforGMM' , 'ProbabilityThreshold' , 'Minute Threshold') )
  BCParameters[[1]] <- 200
  BCParameters[[2]] <- 25
  BCParameters[[3]] <- 0.95
  BCParameters[[4]] <- 8
  
  if(BCOptions[[2]] ==  'AFClassifier'){
    Priorprobabilities <- setNames( list(1 , 1 , 1 , 1 , 1 , 1) ,
                                    c('A' ,
                                      'A^c',
                                      'B|A',
                                      'B|A^c',
                                      'B',
                                      'B^c') )
    Priorprobabilities[[ 1 ]] <- 0.2
    Priorprobabilities[[ 2 ]] <- 1 - Priorprobabilities[[1]]
    Priorprobabilities[[ 3 ]] <- 0.0005
    Priorprobabilities[[ 4 ]] <- 0
    Priorprobabilities <- BC_CaculateBeatWiseProbabiltiesAFDetection(Priorprobabilities)
  }
  precomputedfolderpath <- DP_SelectPrecomputedFolder()
  load(paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]] , BCOptions[[4]] ,   'DistributionSummaries'  ,'.RData')))
}
Patinettotest <- 'z305'
DistributionSummaries <- DP_LoadDistributionSummaries(path , Patinettotest)
RPeaksStruct <- DP_LoadRpeaksfile(path , Patinettotest)

MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = Patinettotest)
AdjustedBeliefs <- BC_EstimateGlobalandLocalParametersDisSum(DistributionSummaries)
DiffStructure <- apply(AdjustedBeliefs$W , 2 , function(X){c(0,diff(X))} )

f_i <- BC_CalulateDenistiesGMM(DP_RemoveNaRows(AdjustedBeliefs$W) , LocalDistributionStruct = LocalDistributionStruct)

UpdateMatrix <- array(0 , c(dim(f_i ) , 250))

for(i in 1:dim(UpdateMatrix)[3]){
for(j in 1:dim(UpdateMatrix)[2]){
  if(i < 250){
   UpdateMatrix[,j,i] <- log(f_i[,j])*(1/250)
  }else{
    UpdateMatrix[,j,i] <- log(f_i[,j])
}
   }
}

UpdateMatrix <- exp(UpdateMatrix)

P_AF <- matrix(Priorprobabilities$B ,dim(f_i)[1] , 1)
plot((P_AF) , col ='red' , type = 'l' , ylim = c(0,1))

for(index in 1:250){
P_AF[index: (dim(P_AF)[1] ) ] <- P_AF[index: (dim(P_AF)[1]) ]*UpdateMatrix[1:(dim(P_AF)[1] -(index -1)),1 , index]/
                                 (P_AF[index: (dim(P_AF)[1] ) ]*UpdateMatrix[1:(dim(P_AF)[1] -(index -1)),1 , index] +
                                 (1-P_AF[index: (dim(P_AF)[1] ) ])*UpdateMatrix[1:(dim(P_AF)[1] -(index -1)),2 , index])
lines(P_AF , col =rgb(index/250,0,0 , alpha = 0.1))  
}

WeightFunction <- function(distance , n){
  output <- matrix(0 , length(distance) , 1)
  output[,1] <- distance/n 
  output <- apply(output , 1, function(X){min(X , 1)} )
}

BC_CaluluateCulmulativeLikelihoodKernel <- function( W , LocalDistributionStruct, n = 5 , weight = 250 , computationalWeight = 1 ){
  WeightKernel <- WeightFunction(1:n , weight)/computationalWeight
  num_classes <- length(LocalDistributionStruct)
  f_i <- matrix(0 , size(W)[1] , length(LocalDistributionStruct))
  for(i in 1:num_classes){
    f_i[ , i] <- predict( LocalDistributionStruct[[i]] , W , what = c('dens') , logarithm = TRUE)  
  } 
  
  for(i in 1:size(f_i)[2]){
    f_i[f_i[,i] <= 0,i] = 1e-10
    f_i[,i] <- (f_i[,i] + rollapply( f_i[,i]/computationalWeight , width = n , function(X){sum(X*WeightKernel)} , align = 'right' , na.pad =TRUE))
  }
  return( f_i )  
}

