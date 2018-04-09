# Load in libraries and settings

pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# load data % Tested on z1139.Rdata
FileLocation <- choose.files()
load(FileLocation)

# take a number of hours before the end of recording to process
HoursBeforeEnd = 1;
WaveData <- WaveData[ ( WaveData$Date )> (WaveData$Date[length(WaveData$Date)] - HoursBeforeEnd*(60^2))  , 1:2]

# Check first chunk of data.
regiontocheck <- c(1:1000);
t <- WaveData$Date[ regiontocheck ]
f_t = WaveData$Value[ regiontocheck ]

RWaveExtractedDataTest <- RPeakExtraction(t, f_t)
par(mfrow = c(3 , 1))
plot(t , f_t , type = 'l', ylab="H-z")
points( RWaveExtractedDataTest[,1] ,   RWaveExtractedDataTest[,2] , col = 'blue' )
plot(   RWaveExtractedDataTest[,1] ,   RWaveExtractedDataTest[,2] , xlab="t" , ylab="R-Amplitude" )
plot(   RWaveExtractedDataTest[,1] ,   RWaveExtractedDataTest[,3] , xlab="t" , ylab="R-R Times" )

# If the plot above look good, run this section.

t <-  WaveData$Date
f_t <- WaveData$Value
RWaveExtractedData <- RPeakExtraction(t, f_t)

# Save file
SaveLocation <- choose.dir()
SaveLocation <- paste0(SaveLocation , '\\ECG1_RwaveExtraction_z1139.RData' )
save('RWaveExtractedData' , file = SaveLocation)
