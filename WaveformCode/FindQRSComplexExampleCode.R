pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# load data % Tested on z1126.Rdata
load(choose.files(caption = "Select ECG.Rdata file"))

# take a number of hours before the end of recording
HoursBeforeEnd = 6;
WaveData <- WaveData[ ( WaveData$Date )> (WaveData$Date[length(WaveData$Date)] - HoursBeforeEnd*(60^2))  , 1:2]

# create wavelet filter object
Filter = wt.filter(filter = "d6" , modwt=TRUE, level=1)

incriment <- 1200
QRSoutput <- FindQRSComplex(WaveData[1:incriment,] , Filter)

par(mfrow = c(1 , 1))
plot(WaveData[1:incriment,1] , WaveData[1:incriment,2] , type = 'l' , xlab = 't' , ylab = 'Hz')
points(QRSoutput$Qt , QRSoutput$QA  , col = 'blue')
points(QRSoutput$St , QRSoutput$SA  , col = 'green')
title('ECG WaveForm')
abline(0,0)
points(QRSoutput$Rt , QRSoutput$RA  , col = 'red')
