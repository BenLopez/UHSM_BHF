# Script to batch process wave data.

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
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you check for processed files and only process waveforms which have not been processed?')
  if(UserResponse == 'NO'){
    checkforprocessedfiles <- 0
  }else{
    checkforprocessedfiles <- 1
  }
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you to only process AFib patients?')
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


listofpatientswithnotgoodtimes <- list()
counter2<-1
for( ii in  1:length(listAllPatients) ){
  # Some set up for the parameters used in the loop.
  {processdata <- 1
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
  
  for( jj in 1:length(FilestoProcess) ){
    
    if(DP_CheckECGreducedfilesprocessed(path , listAllPatients[ii] , paste0(FilestoProcess[jj] , '_reduced') ) == FALSE || checkforprocessedfiles == 0){  
      print( paste0('Loading patient ' , listAllPatients[ii] , ' ' , FilestoProcess[jj] , '.'  ) )
      WaveData <-DP_LoadECG(path , listAllPatients[ii] , numberrep , FilestoProcess[jj] )
      print( paste0('Patient ' , listAllPatients[ii] , ' ' , FilestoProcess[jj]  , ' loaded.' ) ) 
      crop = 1
    }
    if(DP_CheckECGreducedfilesprocessed(path , listAllPatients[ii] , paste0(FilestoProcess[jj] , '_reduced')) & checkforprocessedfiles == 1){  
      print('Loading reduced file.')
      WaveData <- DP_LoadECGReduced(path , listAllPatients[ii] , numberrep , FilestoProcess[jj] )
      timestamp <- c(WaveData$Date[1] , WaveData$Date[length(WaveData$Date)])
      print('Reduced file loaded.')
      crop = 0
    }
    
    #sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[ii])
    
    if( crop == 1 ){
      print('Cropping ECG.') 
      # Processing ECGI for a patient in AF.  
      {if( DP_CheckIfAFPatient(sub_pat) & jj == 1 ){
        timeindex <-  DP_CalculateTimelimits( sub_pat ,  HoursBeforeandAfter)
        WaveData  <-  DP_CropWaveData(WaveData ,  timeindex , HoursBeforeandAfter)
        WaveData  <-  PE_ReturnWaveformwithPositiveOrientation(WaveData)
        timestamp <-  c( WaveData$Date[1] , WaveData$Date[length(WaveData$Date)] )
      }
      # If non AF Patient find n hour period with good data.  
      if( DP_CheckIfAFPatient(sub_pat) == FALSE & jj == 1){
        if( jj == 1 ){
          counter <- 1
          
          while(counter < 100){
            indent <- as.numeric(abs(median(diff(WaveData$Date[1:min(100000 , length(WaveData$Date))]) , na.rm = FALSE)))
            
            if(length(WaveData$Date) < ( (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter)*(60^2)/indent) ){
              counter <- 101
              break}    
            
            
            timeindex <- sample( (1 + (HoursBeforeandAfter$numberhoursbefore*(60^2)/indent)):(length(WaveData$Date) - (HoursBeforeandAfter$numberhoursafter*(60^2)/indent)) , 1  )
            timedifference <- abs(difftime( WaveData$Date[timeindex - (HoursBeforeandAfter$numberhoursbefore*(60^2)/indent)]  , WaveData$Date[timeindex + (HoursBeforeandAfter$numberhoursafter*(60^2)/indent)] , units ='hours'))
            
            if(timedifference > (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter + 1)){
              counter <- counter + 1  
              next}else{
                WaveData <-  PE_ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData ,  timeindex , HoursBeforeandAfter))
                timestamp <- c(WaveData$Date[1],WaveData$Date[length(WaveData$Date)]) 
                break}
          }
          
          if(counter > 99){
            print(paste0('No good time period found for patient ' , listAllPatients[ii] ,'.') ) 
            listofpatientswithnotgoodtimes[[counter2]] <- listAllPatients[[ii]]
            counter2 <- counter2
            processdata <- 0    
            break}  
        }
      } 
      if( jj > 1 ){
        tmp <- DP_CropWaveData(WaveData ,  timestamp , data.frame(numberhoursbefore = 0 , numberhoursafter = (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter)))
        if( nrow(tmp) == 0 ){ 
          tmp <- DP_LoadECGReduced(path , listAllPatients[ii] , numberrep , 1 ) }else{
            WaveData <-  PE_ReturnWaveformwithPositiveOrientation(tmp)}
      }}
      print('ECG cropped.')   
      Time <- WaveData$Date[1]  
    }
    
    if(processdata == 1){  
      print('Extracting Rpeaks.')
      tmp <- PE_RPeakExtractionWavelet( PE_ReturnWaveformwithPositiveOrientation(WaveData) , Filter = wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5)
      if(nrow(tmp) < 1000){
        processdata = 0
        break}
      outputdata[[jj]] <- PE_CleanRpeaks( tmp , 2 )
      rm(tmp)
      print('Rpeaks extracted.')
    }
    print('Saving Waveform.')
    save( WaveData , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_' , FilestoProcess[jj]  , '_reduced.RData' ) )
    print('Waveform saved.')
    rm(WaveData)  
    
  }
  
# Combine peak informtion from all three ECGs.
  if(processdata == 1){ 
    print('Combining Rpeaks.')
    ECGs <- DP_LoadReducedECGs(path , listAllPatients[[ii]] , numberrep , FilestoProcess)
    outputdata[[4]] <- sub_pat 
    outputdata <- setNames( outputdata , c(FilestoProcess , 'Meta_Data') )
    outputdata[[5]] <- PE_MultipleECGRPeaks(outputdata , ECGs = ECGs)
    outputdata <- setNames( outputdata , c(FilestoProcess , 'Meta_Data' , 'RRCombined') )
    print('Saving output.')
    save( outputdata , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_RPeaks.RData' ) )
    print('Output saved.')
    print('Rpeaks combined.')
    rm(outputdata)
  }
  
  DP_WaitBar(ii/length(listAllPatients))
}

