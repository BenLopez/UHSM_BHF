pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# load data % Tested on z1139.Rdata
load(choose.files(caption = "Select ECG.Rdata file"))

# take a number of hours before the end of recording
HoursBeforeEnd = 6;
WaveData <- WaveData[ ( WaveData$Date )> (WaveData$Date[length(WaveData$Date)] - HoursBeforeEnd*(60^2))  , 1:2]

# create wavelet filter object
Filter = wt.filter(filter = "d6" , modwt=TRUE, level=1)

incriment <- 1200
RWaveExtractedData<-RPeakExtractionWavelet( WaveData[1:incriment,] , Filter)

par(mfrow = c( 3 , 1))
plot(WaveData[1:incriment , 1] , WaveData[1:incriment , 2], type = 'l', xlab = 't' , ylab = 'Hz')
points(RWaveExtractedData$t , RWaveExtractedData$RA  , col = 'blue', xlab = 't' , ylab = 'Hz')
title('ECG Wave Form')
plot(RWaveExtractedData$t , RWaveExtractedData$RA, xlab = 't' , ylab = 'Hz')
title('R-peak Amplitudes')
plot(RWaveExtractedData$t  , RWaveExtractedData$RR , xlab = 't' , ylab = 'dt' , ylim = c(0.6 , 0.85))
title('R-peak Times')
