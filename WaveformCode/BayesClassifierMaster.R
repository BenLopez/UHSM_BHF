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
BCOptions[[2]] <- 'AFPreAFClassifier' 
BCOptions[[3]] <- 'GMM' 
BCOptions[[4]] <- 'GMM'
BCOptions[[5]] <- 'Yes'

BCParameters <- setNames( list(1,1,1,1 , 1) , c( 'TS_Likelihood_clique' , 'NumberComponentsforGMM' , 'ProbabilityThreshold' , 'Minute Threshold' ,'GlobalNumberComponentsforGMM' ) )
BCParameters[[1]] <- 100
BCParameters[[2]] <- 25
BCParameters[[3]] <- 0.6
BCParameters[[4]] <- 6
BCParameters[[5]] <- 1



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

precomputedfolderpath <- DP_SelectPrecomputedFolder()
}
}

# Database is a list where each element is a matrix with data for the training population. 
# Rows are observations and columns are elements if an observation vector.
# Training
source( 'BC_TrainingScript.R' )
source( 'BC_GlobalProbabilisticCalibration.R' )
source( 'BC_LocalProbabilisticCalibration.R' )
