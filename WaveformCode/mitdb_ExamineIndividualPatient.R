pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))

source("LibrariesAndSettings.R" , print.eval  = TRUE )


files <- choose.files(multi = FALSE)


  tmp <- readMat(files)
  waveData <- setNames(  data.frame(tmp$Data[[2]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2)] ,
                                    tmp$Data[[1]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2) , 1] ) , c('Date' , 'Value'))
  waveData <-ReturnWaveformwithPositiveOrientation( waveData )
  Annotation <- setNames(data.frame( round(tmp$Data[[3]]/2) ,  tmp$Data[[4]])  , c('Index' , 'Code'))
  stdthresh <- 2.5
  RData <- RPeakExtractionWavelet( waveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = stdthresh)
  Results <- mitdb_ComputePeakresults( WaveData , Annotation , stdthresh = stdthresh )
  

interestingpoint <- 3
centre <- Results$Missedpeaks[interestingpoint]
regionofinterest <- max(1 , centre - 1000): min(length(waveData[ , 1]) , centre + 1000)
plot(waveData[regionofinterest , 1] ,waveData[regionofinterest , 2] , type = 'l'  )
points(RData$t , RData$RA , col = 'red' ,pch = 20)
points(waveData[Annotation$Index , 1] , waveData[Annotation$Index , 2] , col ='blue' )
title(paste0('Patient ' , substr(files , 21 , nchar(files) - 4)))
abline(v=waveData[centre , 1])
