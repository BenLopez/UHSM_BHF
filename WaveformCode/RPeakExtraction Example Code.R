# Load in libraries and settings

pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# load data % Tested on z1139.Rdata
load(choose.files())

# take a number of hours before the end of recording
HoursBeforeEnd = 1;
WaveData <- WaveData[ ( WaveData$Date )> (WaveData$Date[length(WaveData$Date)] - HoursBeforeEnd*(60^2))  , 1:2]

incriment <- 12000
t <- WaveData$Date[1:incriment]
f_t = WaveData$Value[1:incriment]

# function to perfrom peak extraction 
functionoutputtest <- RPeakExtraction(t, f_t)

# Diagnostic plots
par(mfrow = c(3 , 1))
plot(t , f_t , type = 'l')
points( functionoutputtest[,1] , functionoutputtest[,2] , col = 'blue' )
plot( functionoutputtest[,1] ,   functionoutputtest[,2] , xlab="t", ylab="R-Amplitude" )
plot( functionoutputtest[,1] ,   functionoutputtest[,3] , xlab="t", ylab="R-R Times" )

