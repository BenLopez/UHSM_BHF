{if(file.exists('CheckforDefaultsScript.R')){
  source('CheckforDefaultsScript.R')
}else{
  pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
  source("LibrariesAndSettings.R" , print.eval  = TRUE )
  DP_LoadPatientIndex()
  DP_ChooseDataReps()
  FilestoProcess <- DP_ChooseECGstoProcess() 
  HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
}
}
options(warn=-1)
source( 'ASWF_ChooseLoadandProcessPatient.R' )
source( 'ASWF_InteractivePatientAnalysis.R'  )
