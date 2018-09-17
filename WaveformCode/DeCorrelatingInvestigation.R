
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
    Priorprobabilities[[ 3 ]] <- 0.241
    Priorprobabilities[[ 4 ]] <- 0
    Priorprobabilities <- BC_CaculateBeatWiseProbabiltiesAFDetection(Priorprobabilities)
  }
  precomputedfolderpath <- DP_SelectPrecomputedFolder()
  load(paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]] , BCOptions[[4]] ,   'DistributionSummaries'  ,'.RData')))
}
Patinettotest <- 'z1026'
DistributionSummaries <- DP_LoadDistributionSummaries(path , Patinettotest)
RPeaksStruct <- DP_LoadRpeaksfile(path , Patinettotest)

MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = Patinettotest)
AdjustedBeliefs <- BC_EstimateGlobalandLocalParametersDisSum(DistributionSummaries)

index <- 3
window = 800
el = 400
Sigma_n = 0.1
X <- as.matrix(1:window)
Xstar <- window
CXXstar <- CF_ExponentialFamily( Xstar , X , l = el , p = 2)
CXX <- CF_ExponentialFamily( X , X , l = el , p = 1)
A = CXXstar%*%solve(CXX + Sigma_n*diag(dim(CXX)[1]))
A = A/sum(A)

ResolvedVariance = (CXX%*%solve(CXX + Sigma_n*diag(dim(CXX)[1])))%*%t(CXX)

L = solve(t(chol(ResolvedVariance)))

tmp <- AdjustedBeliefs$W[,index]
tmp[is.na(tmp)] <- mean(AdjustedBeliefs$W[,index][!is.na(AdjustedBeliefs$W[,index])])

blah <- rollapply(tmp , width = window , FUN = function(X){A%*%X} , na.pad = TRUE  , align = 'right')
plot(blah , type ='l' )
lines(AdjustedBeliefs$W[,index] , col ='red')
