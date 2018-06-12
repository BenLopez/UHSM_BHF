
{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}


HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
OutputDirectory <- choose.dir()
Listoffiles <- choose.files()
FileNames <- lapply(Listoffiles , function(X){substr(basename(X) , start = 1 , stop = (nchar( basename(X) ) - 4) )})

for( i in 1:length(FileNames) )
{

  dir.create( paste0(OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out'  ) ,  recursive = TRUE)
  tmp <- readMat(Listoffiles[[i]])
  WaveData <- data.frame(Date = tmp[[1]][[1]][,1] , Values = tmp[[1]][[1]][,2])
  AFlocations <- tmp[[1]][[3]]

  if(sum(AFlocations > 5*(60^2)) > 0){
    indexOI <- which.min(abs(WaveData$Date - AFlocations[AFlocations > 5*(60^2)][1]))
  }
  if(sum(AFlocations > 5*(60^2)) == 0){
    indexOI <- 5*(60^2)
  }
  WaveData <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , indexOI , HoursBeforeandAfter))
  
  print('Saving ECGI')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , 'ECGI_', FileNames[[i]] , '.RData') )
  print('ECGI Saved')
  
  WaveData <- data.frame(Date = tmp[[1]][[1]][,1] , Values = tmp[[1]][[1]][,3])
  WaveData <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , indexOI , HoursBeforeandAfter))
  
  print('Saving ECGII')
  save(WaveData , file = paste0( OutputDirectory , '\\' ,  FileNames[[i]] , '\\Zip_out\\' , 'ECGII_', FileNames[[i]] , '.RData') )
  print('ECGII Saved')
}


