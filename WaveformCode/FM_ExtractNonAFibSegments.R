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
    set.seed( 1 )
  }
}

NSRNonImplausibleSets <- FM_GetNonImplausibleSetsFromlabelledDataset()

PreAFNonImplausibleSet <- matrix(0 , 0,10)
NPreAFNonImplausibleSet <- matrix(0 , 0,10)
AFLogical <- matrix(NA , length(listAllPatients))
for(PatientID in listAllPatients ){
  if( !file.exists(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData')) ){ next
  }
  if(DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017 , PatientID))){
  PreAFNonImplausibleSet <- unique(rbind(PreAFNonImplausibleSet ,FM_ExtractNonAFibNonImplausibleSets(PatientID  , PatIndex2017) ) ) 
  }else{
  NPreAFNonImplausibleSet <- unique(rbind(NPreAFNonImplausibleSet ,FM_ExtractNonAFibNonImplausibleSets(PatientID  , PatIndex2017) ) ) 
}
}

BC_PlotPairs( PreAFNonImplausibleSet[sample(1:dim(PreAFNonImplausibleSet)[1] , 1000) , ] , alpha = 0.05 )
BC_PlotPairs( NPreAFNonImplausibleSet[sample(1:dim(NPreAFNonImplausibleSet)[1] , 1000) , ] , alpha = 0.05 )


BC_PlotPairsFromTwoVariables(PreAFNonImplausibleSet[sample(1:dim(PreAFNonImplausibleSet)[1] , 1000) , ] , PreAFNonImplausibleSet[sample(1:dim(PreAFNonImplausibleSet)[1] , 1000) , ]  , alpha = 0.05 )

BC_PlotPairsFromTwoVariables(PreAFNonImplausibleSet[sample(1:dim(PreAFNonImplausibleSet)[1] , 1000) , ] , NPreAFNonImplausibleSet[sample(1:dim(NPreAFNonImplausibleSet)[1] , 1000) , ] , alpha = 0.05 )
BC_PlotPairsFromTwoVariables(PreAFNonImplausibleSet[sample(1:dim(PreAFNonImplausibleSet)[1] , 1000) , ],NSRNonImplausibleSets[sample(1:dim(NSRNonImplausibleSets)[1] , 1000),] , alpha = 0.05 )
