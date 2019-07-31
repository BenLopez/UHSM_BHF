{if(file.exists('CheckforDefaultsScript.R')){
  source('CheckforDefaultsScript.R')
}else{
  pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
  source("LibrariesAndSettings.R" , print.eval  = TRUE )
  PatIndex2017 <- DP_LoadPatientIndex()
  DP_ChooseDataReps()
  FilestoProcess <- DP_ChooseECGstoProcess() 
  HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
  #HowtoFilterops <- read.csv(file.choose( caption = "Select files" ) , stringsAsFactors = FALSE  )
  set.seed(1)
}
}

# Options for processing.
{
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to check for processed files and only process waveforms which have not been processed?')
  if(UserResponse == 'NO'){
    checkforprocessedfiles <- 0
  }else{
    checkforprocessedfiles <- 1
  }
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to only process AFib patients?')
  if(UserResponse == 'NO'){
    OnlyAFPatientes <- 0
  }else{
    OnlyAFPatientes <- 1
  }
  
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to filter waveforms?')
  if(UserResponse == 'YES'){
    HowtoFilterops <- read.csv(file.choose(  ) , stringsAsFactors = FALSE  )
    listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  } 
}

listAllPatients <- select.list(listAllPatients , graphics = TRUE , multiple = TRUE , title = 'Select patients to process')

HoursBeforeandAfter$numberhoursbefore = 0
HoursBeforeandAfter$numberhoursafter = 6


listofpatientswithnotgoodtimes <- list()
counter2<-1

for( ii in  275:length(listAllPatients) ){
  # Some set up for the parameters used in the loop.
  { processdata <- 1
    outputdata  <- list()  
    rm(timestamp)
    sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[ii])}
  
  # Logic to skip patient if they should not (or can not) be processed.
  {if( DP_checkfilesprocessed(path , listAllPatients[ii] , FilestoProcess[1]) == 0  ){
    print(paste0('Skipping patient ', listAllPatients[ii] ,' as no ECGI processed.'))
    next
  }
    if(OnlyAFPatientes == 1){
      if(DP_CheckIfAFPatient(sub_pat) == FALSE){
        print(paste0('Skipping patient ' , listAllPatients[ii] , '.'))
        next
      }
    }}
  
  if(DP_CheckFileExists(path = path , listAllPatients[[ii]] , paste0(listAllPatients[[ii]] , '_ECGI' )) ){next}
  
  for( jj in 1:length(FilestoProcess) ){
    
   
    if(DP_CheckECGreducedfilesprocessed(path , listAllPatients[ii] , paste0(FilestoProcess[jj] , '_reduced') ) == FALSE || checkforprocessedfiles == 0){  
      print( paste0('Loading patient ' , listAllPatients[ii] , ' ' , FilestoProcess[jj] , '.'  ) )
      WaveData <- DP_LoadECG(path , listAllPatients[ii] , numberrep , FilestoProcess[jj] )
      print( paste0('Patient ' , listAllPatients[ii] , ' ' , FilestoProcess[jj]  , ' loaded.' ) ) 
      crop = 1
    }
    
    
      print('Cropping ECG.') 
      # Processing ECGI for a patient in AF.  
      if(jj == 1){
        timeindex <-  1000
        WaveData  <-  DP_CropWaveData(WaveData ,  timeindex , HoursBeforeandAfter)
        timestamp <-  c( WaveData$Date[1] , WaveData$Date[length(WaveData$Date)] )
      } 
      if( jj > 1 ){
          tmp <- DP_CropWaveData(WaveData ,  timestamp , data.frame(numberhoursbefore = 0 , numberhoursafter = (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter)))
          if( nrow(tmp) == 0 ){ 
            tmp <- DP_LoadECGReduced(path , listAllPatients[ii] , numberrep , 1 ) }else{
              WaveData <-  PE_ReturnWaveformwithPositiveOrientation(tmp)}
        }
      print('ECG cropped.')   
      Time <- WaveData$Date[1]  
    
    
    print('Extracting Rpeaks.')
    tmp <- PE_RPeakExtractionWavelet( PE_ReturnWaveformwithPositiveOrientation(WaveData) , Filter = wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5)
    if(nrow(tmp) < 1000){
      processdata = 0
      break}
    outputdata[[jj]] <- PE_CleanRpeaks( tmp , 2 )
    rm(tmp)
    print('Rpeaks extracted.')
    
    print('Saving Waveform.')
    save( WaveData , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_' , FilestoProcess[jj]  , '_reducedstart.RData' ) )
    print('Waveform saved.')
    rm(WaveData)  
    
  }
  
  
  # Combine peak informtion from all three ECGs.
  if(processdata == 1){ 
    print('Combining Rpeaks.')
    ECGs <- DP_LoadStartReducedECGs(path , listAllPatients[[ii]] , numberrep , FilestoProcess)
    outputdata[[4]] <- sub_pat 
    outputdata <- setNames( outputdata , c(FilestoProcess , 'Meta_Data') )
    outputdata[[5]] <- PE_MultipleECGRPeaks(outputdata = outputdata , ECGs = ECGs)
    outputdata <- setNames( outputdata , c(FilestoProcess , 'Meta_Data' , 'RRCombined') )
    print('Saving output.')
    save( outputdata , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_RPeaksStart.RData' ) )
    print('Output saved.')
    print('Rpeaks combined.')
    rm(outputdata)
  }
  
  DP_WaitBar(ii/length(listAllPatients))
}
