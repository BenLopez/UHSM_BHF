
{
if(exists('BCOptions') == FALSE && exists('LocalDistributionStruct') == FALSE){
  print('Training has not been completed, running Training Script.')
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
    
    # Run script to create parameters and options.
    source('BC_ChooseParametersandOptions.R')
    
  }
  
  # Database is a list where each element is a matrix with data for the training population. 
  # Rows are observations and columns are elements if an observation vector.
  # Training
  source( 'BC_TrainingScript.R' )
}
  
Patinettotest <- DP_choosepatient(listAllPatients = listAllPatients)
source('BC_LoadDataandTestSinglePatient.R')
if(nrow(AFLocations)>0){
  timetoview <- AFLocations$Start[1]
  }else{
  timetoview <- DP_SelectTimetoview(ECGs$ECGI$Date)
}

lengthtoview <- BC_SelectTimesofECGtoview()
source('BC_CreatePlots.R')

if( nrow(AFLocations)>0 ){
dev.off()  
  incidencedetected <- 1
  UserResponse <- winDialog(type = c('yesno') , message = paste0('Heart rate analysis has detected ',nrow(AFLocations) , ' an Atrial Fibrillation episode. Would you like to view P-wave graphics.') )
  while( (incidencedetected <= nrow(AFLocations)) && (UserResponse == 'YES') ){
      source('BC_PwaveGraphicalAnalysis.R')
      incidencedetected <- incidencedetected + 1
      if(incidencedetected <= nrow(AFLocations)){
        UserResponse <- winDialog(type = c('yesno') ,
                                  message = paste0('Would you like to view the next incidence?') )
        
      }
      
  }
if(UserResponse == 'YES'){source('BC_CreatePlots.R')}
  
}
  
if( nrow(AFLocations)>0 ){
UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to view Pairs analysis?')
if(UserResponse == 'YES'){source('BC_PairsAnalysis.R')}
}

UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to view another patient?')
if(UserResponse == 'YES'){source('BC_TestingScript.R')}
if(UserResponse == 'NO'){print('GUI has been exited. Run BC_TestingScript.R to continue viewing patient data.')}



}





