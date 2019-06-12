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

HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
OutputDirectory <- choose.dir()
Listoffiles <- choose.files()
FileNames <- lapply(Listoffiles , function(X){substr(basename(X) , start = 1 , stop = (nchar( basename(X) ) - 4) )})

for( i in 1:length(FileNames) ){

  dir.create( paste0(OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out'  ) ,  recursive = TRUE)
  tmp <- readMat(Listoffiles[[i]])
  timetmp <- as.POSIXct(Sys.time())
  WaveData <- data.frame(Date = as.POSIXct( tmp[[1]][[1]][,1] , origin = timetmp ) , Value = tmp[[1]][[1]][,2])
  AFStartTimes <-  as.POSIXct( tmp[[1]][[3]], origin = timetmp)
  AFEndTimes <-  as.POSIXct( tmp[[1]][[4]], origin = timetmp)
  
  if(sum(difftime(AFStartTimes , WaveData$Date[1] , units = 'hours') > 6) > 0){
    tmp1 <- AFStartTimes[difftime(AFStartTimes , WaveData$Date[1] , units = 'hours') > 5]
    tmp2 <- AFEndTimes[difftime(AFEndTimes , WaveData$Date[1] , units = 'hours') > 5]
    TimeIndexuse <- which.max(abs(AFStartTimes - AFEndTimes))
    indexOI <- which.min(abs(WaveData$Date - AFStartTimes[TimeIndexuse] ))
  }
  if(sum(difftime(AFStartTimes , WaveData$Date[1] , units = 'hours') > 6) == 0){
    indexOI <- which(difftime(WaveData$Date,WaveData$Date[1] , units = 'hours')>6)[1]
  }
  
  WaveData <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , indexOI , HoursBeforeandAfter))
  
  print('Saving ECGI')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , FileNames[[i]] , '_ECGI_Reduced.RData') )
  print('ECGI Saved')
  
  print('Saving ECGII')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , FileNames[[i]] , '_ECGII_Reduced.RData') )
  print('ECGII Saved')
  
  WaveData <- data.frame(Date = as.POSIXct( tmp[[1]][[1]][,1] , origin = timetmp )  , Value = tmp[[1]][[1]][,3])
  WaveData <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , indexOI , HoursBeforeandAfter))
  
  print('Saving ECGIII')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , FileNames[[i]] , '_ECGIII_Reduced.RData') )
  print('ECGIII Saved')
  
  print('Saving AFLocations')
  AFlocations = data.frame(AFStartTimes = AFStartTimes , AFEndTimes = AFEndTimes)
  save(AFlocations , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , FileNames[[i]] , '_AFTimes.RData') )
  print('AFLocations Saved')
  DP_WaitBar(i/length(FileNames))
  
}


