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
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

# A structure for options of what type to process. 
# BCOptions[['DataType']] = 'DistributionSummaries' , 'CDFs'
# BCOptions[['Classifier Type']] = 'AFClassifier' , 'AFPreAFClassifier'
# BCOptions[[ 'Density Estimation Global']] = 'MVN' , 'GMM'
# BCOptions[['Density Estimation local']] = 'MVN' , 'GMM'

{BCOptions <- setNames( list(1 , 1 , 1 , 1 , 1) , c('DataType' ,
                                                    'ClassifierType' , 
                                                    'DensityEstimationGlobal', 
                                                    'DensityEstimationlocal' , 
                                                    'GlobalUpdate') )
BCOptions[[1]] <- 'DistributionSummaries'
BCOptions[[2]] <- 'AFClassifier' 
BCOptions[[3]] <- 'MVN' 
BCOptions[[4]] <- 'GMM'
BCOptions[[5]] <- 'Yes'

BCParameters <- setNames( list(1,1,1,1 , 1) , c( 'TS_Likelihood_clique' , 'NumberComponentsforGMM' , 'ProbabilityThreshold' , 'Minute Threshold' ,'GlobalNumberComponentsforGMM' ) )
BCParameters[[1]] <- 1
BCParameters[[2]] <- 25
BCParameters[[3]] <- 0.6
BCParameters[[4]] <- 6
BCParameters[[5]] <- 1


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
}
}

# Database is a list where each element is a matrix with data for the training population. 
# Rows are observations and columns are elements if an observation vector.
# Training
source( 'BC_TrainingScript.R' )
{
x11(20,14)
par(mfrow = c(2 , 5))
SampleFromModel <- BC_SampleGMM(LocalDistributionStruct[[1]] , 100000)
SampleFromModel2 <- BC_SampleGMM(LocalDistributionStruct[[2]] , 100000)
for(variabletoview in c(1:10)){
  tmp <- hist(DataBase[[1]][sample(1:size(DataBase[[1]])[1] , 100000) , variabletoview], col=rgb(0,0,1,alpha = 0.5) ,
              main = paste0(AFD_CreateDistributionSummaryNames()[variabletoview] , ' Histogram') ,
              xlab = AFD_CreateDistributionSummaryNames()[variabletoview] , freq = FALSE )
  tmp2 <-hist(DataBase[[2]][sample(1:size(DataBase[[2]])[1] , 100000) , variabletoview], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(min(DataBase[[2]][!is.na(DataBase[[2]][ , variabletoview]) , variabletoview] ) , tmp$breaks, max(DataBase[[2]][!is.na(DataBase[[2]][ , variabletoview]) , variabletoview]) ))
  tmp3 <-hist(SampleFromModel[ , variabletoview], col=rgb(0,1,0,alpha =0.5), add=T , freq = FALSE ,
       breaks = c(min(SampleFromModel[!is.na(SampleFromModel[ , variabletoview]) , variabletoview] ) , tmp$breaks, max(SampleFromModel[!is.na(SampleFromModel2[ , variabletoview]) , variabletoview]) ))
  tmp4 <-hist(SampleFromModel2[ , variabletoview], col=rgb(0.5,0.5,0.5,alpha =0.5), add=T , freq = FALSE ,
       breaks = c(min(SampleFromModel2[!is.na(SampleFromModel2[ , variabletoview]) , variabletoview] ) , tmp$breaks, max(SampleFromModel2[!is.na(SampleFromModel2[ , variabletoview]) , variabletoview]) ))
  
}
BC_plotValidateDensityEstimationMarginalHistograms(DataBase[[1]] , DataBase[[2]] , SampleFromModel , SampleFromModel2 )   
}


source( 'BC_GlobalProbabilisticCalibration.R' )
source( 'BC_LocalProbabilisticCalibration.R' )


