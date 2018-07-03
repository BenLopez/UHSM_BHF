{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
OutputDirectory <- choose.dir()
Listoffiles <- choose.files()
FileNames <- lapply(Listoffiles , function(X){substr(basename(X) , start = 1 , stop = (nchar( basename(X) ) - 4) )})

for( i in 1:length(FileNames) ){
  
  dir.create( paste0(OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out'  ) ,  recursive = TRUE)
  tmp <- readMat(Listoffiles[[i]])
  timetmp <- as.POSIXct(Sys.time())
  WaveData <- data.frame(Date =as.POSIXct(tmp[[1]][[1]][,1] , origin =  timetmp ) , Value = tmp[[1]][[1]][,2])
  AFlocations <- tmp[[1]][[3]]
  
  indexOI <- round(length(WaveData$Date)/2)
  WaveData <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , indexOI , HoursBeforeandAfter))
  
  print('Saving ECGI')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , FileNames[[i]] , '_ECGI_Reduced.RData') )
  print('ECGI Saved')
  
  print('Saving ECGIII')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , FileNames[[i]] , '_ECGIII_Reduced.RData') )
  print('ECGIII Saved')
  
  WaveData <- data.frame(Date = as.POSIXct(tmp[[1]][[1]][,1]  , origin =  timetmp ), Value = tmp[[1]][[1]][,3])
  WaveData <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , indexOI , HoursBeforeandAfter))
  
  print('Saving ECGII')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , FileNames[[i]] , '_ECGII_Reduced.RData') )
  print('ECGII Saved')
}
