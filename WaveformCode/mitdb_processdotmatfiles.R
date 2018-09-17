## Script to process.mat files into WaveData

##### Preamble #####
{if(file.exists('CheckforDefaultsScript.R')){
  source('CheckforDefaultsScript.R')
}else{
  pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
  source("LibrariesAndSettings.R" , print.eval  = TRUE )
  DP_LoadPatientIndex()
  DP_ChooseDataReps()
}
}
##### Preamble ##### 

files <- choose.files()
for(i in 1:length( files ) ){
  tmp <- readMat(files[[i]])
  
  WaveData <- setNames(  data.frame(tmp$Data[[2]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2)] ,
                                    tmp$Data[[1]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2) , 1] ) , c('Date' , 'Value'))
  WaveData$Date <- as.POSIXct(WaveData$Date , origin = Sys.time())
  WaveData[ , 2] <- WaveData[ , 2] - mean(WaveData[ , 2])
  WaveData <-ReturnWaveformwithPositiveOrientation( WaveData )
  Annotation <- setNames(data.frame( round(tmp$Data[[3]]/2) ,  tmp$Data[[4]])  , c('Index' , 'Code'))
  
  folderpath <- paste0('D:\\mtbd_Waveform' ,substr( files[[i]] ,  nchar(files[[i]]) - 13 , nchar(files[[i]]) - 4 ) , '\\Zip_out')
  dir.create(path = folderpath   , recursive = TRUE)
  
  save(WaveData , file = paste0( folderpath , '\\' , 'ECGI_' ,  file = substr( files[[i]] ,  nchar(files[[i]]) - 12 , nchar(files[[i]]) - 4 ) , '.RData' ))
  save(Annotation , file = paste0( folderpath ,file =  substr( files[[i]] ,  nchar(files[[i]]) - 13 , nchar(files[[i]]) - 4 ) , '_Annotation.RData' ))
  save(WaveData , file = paste0( folderpath , '\\' , 'ECGIII_' ,  file = substr( files[[i]] ,  nchar(files[[i]]) - 12 , nchar(files[[i]]) - 4 ) , '.RData' ))
  
  
  WaveData <- setNames(  data.frame( WaveData$Date ,
                                    tmp$Data[[1]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2) , 2] ) , c('Date' , 'Value'))
  WaveData[ , 2] <- WaveData[ , 2] - mean(WaveData[ , 2])
  WaveData <-ReturnWaveformwithPositiveOrientation( WaveData )
  Annotation <- setNames(data.frame( round(tmp$Data[[3]]/2) ,  tmp$Data[[4]])  , c('Index' , 'Code'))

  save(WaveData , file = paste0( folderpath , '\\' , 'ECGII_' ,  file = substr( files[[i]] ,  nchar(files[[i]]) - 12 , nchar(files[[i]]) - 4 ) , '.RData' ))
  save(Annotation , file = paste0( folderpath , '\\' , 'Annotation_' ,  file = substr( files[[i]] ,  nchar(files[[i]]) - 12 , nchar(files[[i]]) - 4 ) , '.RData' ))
  print(i/length( files ))
}




