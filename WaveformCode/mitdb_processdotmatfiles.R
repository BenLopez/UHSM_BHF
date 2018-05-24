## Script to process.mat files into wavedata
pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))

source("LibrariesAndSettings.R" , print.eval  = TRUE )


files <- choose.files()

for(i in 1:length( files ) )
{
  
  tmp <- readMat(files[[i]])
  
  waveData <- setNames(  data.frame(tmp$Data[[2]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2)] ,
                                    tmp$Data[[1]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2) , 1] ) , c('Date' , 'Value'))
  waveData[ , 2] <- waveData[ , 2] - mean(waveData[ , 2])
  waveData <-ReturnWaveformwithPositiveOrientation( waveData )
  Annotation <- setNames(data.frame( round(tmp$Data[[3]]/2) ,  tmp$Data[[4]])  , c('Index' , 'Code'))
  
  # save(waveData , file = paste0( 'D:\\mtbd_Waveform' , file = substr( files[[i]] ,  nchar(files[[i]]) - 13 , nchar(files[[i]]) - 4 ) , '_ECGI.RData' ))
  # save(Annotation , file = paste0( 'D:\\mtbd_Waveform' ,file =  substr( files[[i]] ,  nchar(files[[i]]) - 13 , nchar(files[[i]]) - 4 ) , '_Annotation.RData' ))
  # RData <- RPeakExtractionWavelet( waveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8)
  
  Results <- mitdb_ComputePeakresults( WaveData , Annotation , stdthresh = 2.5)
  
  save(Results , file = paste0( 'D:\\mtdb_peakdetectionresults' ,file =  substr( files[[i]] ,  nchar(files[[i]]) - 13 , nchar(files[[i]]) - 4 ) , '_PeakDetectionResults.RData' ))
  print(i/length( files ))
}




