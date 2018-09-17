# Script to batch process wave data.

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

# Optionsforprocessing
checkforprocessedfiles <- 0
OnlyAFPatientes <- 0

listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
listAllPatients <- select.list(listAllPatients , graphics = TRUE , multiple = TRUE)
listofpatientswithnotgoodtimes <- list()
counter2<-1
#1:length(listAllPatients)
for( ii in  1:length(listAllPatients) ){
  processdata <- 1
  outputdata  <- list()  
  rm(timestamp)
  sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[ii])
  
  if(OnlyAFPatientes == 1){
    if(is.na(sub_pat$ConfirmedFirstNewAF) || sub_pat$ConfirmedFirstNewAF == 'CNAF'){
      print(paste0('skipping patient ' , listAllPatients[ii]))
      next
    }
  }
  
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
    sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[ii])
    
    if( crop == 1 ){
      print('Cropping ECG.') 
      # Patient in AF  
      if(!is.na(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1])) & sub_pat$ConfirmedFirstNewAF[1] != 'CNAF' & jj == 1){
        #timeindex <- which.min(abs(difftime(WaveData$Date , DP_StripTime(sub_pat$ConfirmedFirstNewAF[1]) , units = 'mins' )))[1]
        timeindex<- c(0,0)
        timeindex[1] <- DP_AddHour(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1] ) , -HoursBeforeandAfter$numberhoursbefore)
        timeindex[2] <- DP_AddHour(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1] ) , HoursBeforeandAfter$numberhoursafter)
        WaveData <-  PE_ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData ,  timeindex , HoursBeforeandAfter))
        timestamp <- c(WaveData$Date[1] , WaveData$Date[length(WaveData$Date)])
      }
      #if(is.na(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1])) & jj == 1){
      #  timeindex <- which.min(abs(difftime(WaveData$Date , DP_StripTime(sub_pat$FirstNewAF), units = 'mins' )))[1]
      #  WaveData <- DP_CropWaveData(WaveData ,  timeindex , HoursBeforeandAfter)
      #}
      if((is.na(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1])) || sub_pat$ConfirmedFirstNewAF[1] == 'CNAF') & jj == 1){
        if(jj == 1){
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
      if(jj > 1){
        tmp <- DP_CropWaveData(WaveData ,  timestamp , data.frame(numberhoursbefore = 0 , numberhoursafter = (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter)))
      if(nrow(tmp) == 0){ 
        tmp <- DP_LoadECGReduced(path , listAllPatients[ii] , numberrep , 1 ) }else{
        WaveData <-  PE_ReturnWaveformwithPositiveOrientation(tmp)}
      }
      print('ECG cropped.')   
      Time <- WaveData$Date[1]  
    }
    
    if(processdata == 1){  
      print('Extracting Rpeaks.')
      tmp <- RPeakExtractionWavelet( PE_ReturnWaveformwithPositiveOrientation(WaveData) , Filter = wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5)
      if(nrow(tmp) < 1000){
        processdata = 0
        break}
      outputdata[[jj]] <- CleanRpeaks(tmp , 2)
      print('Rpeaks extracted.')
    }
    print('Saving Waveform.')
    save( WaveData , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_' , FilestoProcess[jj]  , '_reduced.RData' ) )
    print('Waveform saved.')
    rm(WaveData)  
    
  }
  
#  if(counter > 99){
#   DP_WaitBar(ii/length(listAllPatients))
#    break} 
  
  if(processdata == 1){ 
    print('Combining Rpeaks.')
    ECGs <- DP_LoadReducedECGs(path , listAllPatients[[ii]] , numberrep , FilestoProcess)
    outputdata[[length(outputdata) + 1]] <- sub_pat 
    outputdata <- setNames( outputdata , c(FilestoProcess , 'Meta_Data') )
    outputdata[[length(outputdata) + 1]] <- PE_MultipleECGRPeaks(outputdata , ECGs = ECGs)
    outputdata <- setNames( outputdata , c(FilestoProcess , 'Meta_Data' , 'RRCombined') )
    print('Saving output.')
    save( outputdata , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_RPeaks.RData' ) )
    print('Output saved.')
    print('Rpeaks combined.')
    rm(outputdata)
  }
  DP_WaitBar(ii/length(listAllPatients))
}

