# Load in libraries and settings

pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# load data % Tested on z1139.Rdata
FileLocation <- choose.files(caption="Select .Rdata file of cleaned ECG data")
load(FileLocation)

# take a number of hours before the end of recording to process
HoursBeforeEnd = 6;
WaveData <- WaveData[ ( WaveData$Date )> (WaveData$Date[length(WaveData$Date)] - HoursBeforeEnd*(60^2))  , 1:2]

# Check first chunk of data.
regiontocheck <- c(1:2000);
t <- WaveData$Date[ regiontocheck ]
f_t = WaveData$Value[ regiontocheck ]

RWaveExtractedDataTest <- RPeakExtraction(t, f_t)
par(mfrow = c(3 , 1))
plot(t , f_t , type = 'l', ylab="H-z")
points( RWaveExtractedDataTest[,1] ,   RWaveExtractedDataTest[,2] , col = 'blue' )
plot(   RWaveExtractedDataTest[,1] ,   RWaveExtractedDataTest[,2] , xlab="t" , ylab="R-Amplitude" )
plot(   RWaveExtractedDataTest[,1] ,   RWaveExtractedDataTest[,3] , xlab="t" , ylab="R-R Times" )

# If the plot above looks good, run this section.

t <-  WaveData$Date
f_t <- WaveData$Value
RWaveExtractedData <- RPeakExtraction(t, f_t)

# Save file
Temp <-gregexpr('Zip_out' , FileLocation)
SaveLocation <- substr(FileLocation , 1 , Temp[[1]][1] + nchar('Zip_out'))
Temp <- gregexpr('Zip_out' , FileLocation)
Temp2 <- gregexpr('.RData' , FileLocation)
SaveLocation <- paste0(SaveLocation , substr(FileLocation , Temp[[1]][1] + nchar('Zip_out') +1 , Temp2[[1]][1] -1 ) ,  '_RPeaks.RData' )
rm(Temp,Temp2)
save('RWaveExtractedData' , file = SaveLocation)

## Wavelets


Filter = wt.filter(filter = "d6" , modwt=TRUE, level=1)
features <- attributes(Filter)

lengthoftime <- 0.001
startoftime <- 1.9

#rangetotest = ( ((t > (t[length(t)] - ((startoftime+lengthoftime)*(60^2))))*((t < (t[length(t)] - ((startoftime)*(60^2)))))) == 1 );
rangetotest = c(1:1000);

modoutput <-  modwt( f_t[rangetotest] , Filter , 12)
modoutputattributes <- attributes(modoutput)
W <- slot(modoutput , 'W')
V <- slot(modoutput , 'V')
W <- SetElementsOfListoToZero(W , c(1:2 , 5:12) )
V <- SetElementsOfListoToZero(V , c(1:2 , 5:12) )
slot(modoutput , 'W')<- W
slot(modoutput , 'V')<- V
imodoutput <- imodwt(modoutput, fast=TRUE)

stdresid = imodoutput/sqrt(var(imodoutput))
stdresid[ stdresid < 0] = 0 
Peakslogical = stdresid>2.8
Temp <- f_t[rangetotest]
for (i in 1:length(Peakslogical))
{
  if(Peakslogical[i] == FALSE){ next }
  if(Peakslogical[i] == TRUE && Peakslogical[i + 1] == TRUE )
  {# count logicals
    j <- 1
    while(Peakslogical[i + j] ==  TRUE ){ j <- (j+1) }
  }  
  # Find location of max stdresid in set
  maxindex <- which.max(Temp[i:(i+j)])
  Peakslogical[i:(i+j)] <- FALSE
  Peakslogical[i + (maxindex -1)] <- TRUE
}
Temp <- t[rangetotest]
RPeakLocations <- Temp[Peakslogical == TRUE]
Temp <- f_t[rangetotest]
RAmplitudes <- Temp[Peakslogical == TRUE ]

par(mfrow = c( 2 , 1))
plot(RPeakLocations , RAmplitudes, xlab = 't' , ylab = 'Hz')
title('R-peak Amplitudes')
plot(RPeakLocations  , c(diff(RPeakLocations) , mean(diff(RPeakLocations))) , xlab = 'Index' , ylab = 't' , ylim = c(0.6 , 0.85))
title('R-peak Times')

RWaveExtractedData<-RPeakExtractionWavelet( WaveData[1:10000,] , Filter)

par(mfrow = c( 3 , 1))
plot(WaveData[1:10000 , 1] , WaveData[1:10000 , 2], type = 'l')
points(RWaveExtractedData$t , RWaveExtractedData$RA  , col = 'blue', xlab = 't' , ylab = 'Hz')
title('ECG Wave Form')
plot(RWaveExtractedData$t , RWaveExtractedData$RA, xlab = 't' , ylab = 'Hz')
title('R-peak Amplitudes')
plot(RWaveExtractedData$t  , RWaveExtractedData$RR , xlab = 'Index' , ylab = 't' , ylim = c(0.6 , 0.85))
title('R-peak Times')


par(mfrow = c( 4 , 1))
plot(t[rangetotest] , f_t[rangetotest] , type = 'l' , xlab = 't' , ylab = 'Hz')
title('ECG WaveForm')
abline(0,0)
points(RPeakLocations , RAmplitudes  , col = 'blue')
plot(t[rangetotest] , (imodoutput) , type = 'l', xlab = 't' , ylab = 'Hz')
title('Wavelet Reconstruction (Noise and Drift Removed)')
abline(0,0)
abline(mean(imodoutput + 3*sqrt(var(imodoutput))),0 , col = 'red')
abline(mean(imodoutput - 3*sqrt(var(imodoutput))),0 , col = 'red')
plot(t[rangetotest] , stdresid , type = 'l', xlab = 't' , ylab = 'Std Error')
abline(3,0)
title('Standardised Residuals for Peak Detection')
plot(t[rangetotest] ,-imodoutput + f_t[rangetotest] , type = 'l', xlab = 't' , ylab = 'Hz')
abline(0,0)
title('Difference Between Signal and Wavelet Reconstruction')


multiresolutionanalysis <- attributes(mra(f_t[rangetotest] , Filter ,  12  , boundary = 'periodic' , method = "modwt"))
par(mfrow = c(6 , 1))
plot(t[rangetotest] , f_t[rangetotest] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[1]])] , multiresolutionanalysis$S[[1]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[2]])] , multiresolutionanalysis$S[[2]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[3]])] , multiresolutionanalysis$S[[3]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[4]])] , multiresolutionanalysis$S[[4]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$S[[5]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[1]])] , multiresolutionanalysis$S[[6]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[2]])] , multiresolutionanalysis$S[[7]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[3]])] , multiresolutionanalysis$S[[8]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[4]])] , multiresolutionanalysis$S[[9]] , type = 'l')
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$S[[10]] , type = 'l')

par(mfrow = c(3 , 1))
plot(t[rangetotest] , f_t[rangetotest] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[1]])] , multiresolutionanalysis$D[[1]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[2]])] , multiresolutionanalysis$D[[2]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[3]])] , multiresolutionanalysis$D[[3]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[4]])] , multiresolutionanalysis$D[[4]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[6]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[6]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[7]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[8]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[9]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[10]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[11]] , type = 'l',ann = FALSE)
plot(t[1:length( multiresolutionanalysis$S[[5]])] , multiresolutionanalysis$D[[12]] , type = 'l',ann = FALSE)
