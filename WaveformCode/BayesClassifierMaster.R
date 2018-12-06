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
  precomputedfolderpath <- DP_SelectPrecomputedFolder()

  source('BC_ChooseParametersandOptions.R')

}

# Database is a list where each element is a matrix with data for the training population. 
# Rows are observations and columns are elements if an observation vector.
# Training
source( 'BC_TrainingScript.R' )

if(BCOptions$ProbabilisticCalibration == 'Yes'){
source( 'BC_GlobalProbabilisticCalibration.R' )
source( 'BC_LocalProbabilisticCalibration.R' )
}

source('BC_TestingScript.R')