{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

DP_LoadPatientIndex()
DP_ChooseDataReps()

source( 'ASWF_ChooseLoadandProcessPatient.R' )
source( 'ASWF_InteractivePatientAnalysis.R'  )
