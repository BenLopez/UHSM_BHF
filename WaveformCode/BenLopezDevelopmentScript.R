# Load in libraries and settings

pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

source("LibrariesAndSettings.R" , print.eval  = TRUE )

#
load("C:/Users/Ben/Desktop/WaveformSample/z1139/Zip_out/ECGI_z1139.RData")

# Settings
HoursBeforeEnd = 1;
WaveData <- WaveData[ ( WaveData$Date )> (WaveData$Date[length(WaveData$Date)] - HoursBeforeEnd*(60^2))  , 1:2]

incriment <- 1000
x <- WaveData$Date[1:incriment]
WaveForm1 = WaveData$Value[1:incriment]
RPeakThreshold <-100

# plot wave form
plot(x , (WaveForm1) , type = 'l')

# plot waveform first order differences
lines(x[2:length(x)],abs(diff(WaveForm1)) , type = 'l' , col = 'red')

Spectrum = fft(WaveForm1)
plot(abs(fftshift(fft(WaveForm1))) , type = 'l')
