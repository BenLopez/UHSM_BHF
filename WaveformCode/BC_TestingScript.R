
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
}
  

UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to view another patient?')
if(UserResponse == 'YES'){source('BC_TestingScript.R')}
if(UserResponse == 'NO'){print('GUI has been exited. Run BC_TestingScript.R to continue viewing patient data.')}
}




#AFVev <- AdjustedBeliefs$W[ RPeaksStruct$RRCombined$t > DP_StripTime(MetaData$ConfirmedFirstNewAF) , ]
#AFVev <- AdjustedBeliefs$W[ (RPeaksStruct$RRCombined$t >  AFLocations$Start[1])*((RPeaksStruct$RRCombined$t <  AFLocations$End[1])) == 1 , ]

#x11(30,20)
#SAMPLENUM <- dim(AFVev)[1]
#variable = c(1:11)
#tmp <- rbind(DataBase[[1]][sample(1:size(DataBase[[1]])[1] , SAMPLENUM ),variable] , 
#             DataBase[[2]][sample(1:size(DataBase[[2]])[1] , SAMPLENUM ),variable],
#             AFVev[sample(1:size(AFVev)[1] , SAMPLENUM ),variable])
#al = 0.01
#colvector <- c(rep(rgb(1 ,0 , 0, alpha = al)   , SAMPLENUM) , 
     #          rep(rgb(0 , 0 , 1, alpha = al)   , SAMPLENUM),
      #         rep(rgb(0 , 1 , 0, alpha = al)   , SAMPLENUM))
#pairs(tmp , pch = 16 , col = colvector , labels = AFD_CreateDefaultSettings()$BinlimsScore[variable])

