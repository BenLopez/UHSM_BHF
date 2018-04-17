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

Filter <- wt.filter(filter = "d6" , modwt=TRUE, level=1)
RWaveExtractedData <- RPeakExtractionWavelet( WaveData , Filter)

par(mfrow = c( 3 , 1))
plot(WaveData[14000:18000  , 1] , WaveData[14000:18000 , 2], type = 'l')
points(RWaveExtractedData$t , RWaveExtractedData$RA  , col = 'blue', xlab = 't' , ylab = 'Hz')
title('ECG Wave Form')
plot(RWaveExtractedData$t , RWaveExtractedData$RA, xlab = 't' , ylab = 'Hz')
title('R-peak Amplitudes')
plot(RWaveExtractedData$t  , RWaveExtractedData$RR , xlab = 't' , ylab = 't' , ylim = c(0.6 , 1.2))
title('RR-peak Times')

# Save file
Temp <-gregexpr('Zip_out' , FileLocation)
SaveLocation <- substr(FileLocation , 1 , Temp[[1]][1] + nchar('Zip_out'))
Temp <- gregexpr('Zip_out' , FileLocation)
Temp2 <- gregexpr('.RData' , FileLocation)
SaveLocation <- paste0(SaveLocation , substr(FileLocation , Temp[[1]][1] + nchar('Zip_out') +1 , Temp2[[1]][1] -1 ) ,  '_RPeaks.RData' )
rm(Temp,Temp2)
save('RWaveExtractedData' , file = SaveLocation)


#rangetotest = ( ((t > (t[length(t)] - ((startoftime+lengthoftime)*(60^2))))*((t < (t[length(t)] - ((startoftime)*(60^2)))))) == 1 );
rangetotest = c(1000:2000);
tt <- WaveData$Date[rangetotest]
f_tt <- WaveData$Value[rangetotest]

modoutput <-  modwt( f_tt  , Filter , 12)
modoutputattributes <- attributes(modoutput)
W <- slot(modoutput , 'W')
V <- slot(modoutput , 'V')
W <- SetElementsOfListoToZero(W , c(1:2 , 5:12) )
V <- SetElementsOfListoToZero(V , c(1:2 , 5:12) )
slot(modoutput , 'W')<- W
slot(modoutput , 'V')<- V
imodoutput <- imodwt(modoutput, fast=TRUE)

stdresid <- imodoutput/sqrt(var(imodoutput))
stdresid2 <- stdresid
stdresid[ stdresid < 0] <- 0
Rlogical <- FindLocalTurningPoints( stdresid>2.8 , f_tt )
RPeakLocations <- tt[Rlogical == TRUE]
RAmplitudes <- f_tt[Rlogical == TRUE ]

par(mfrow = c(1 , 1))
plot(tt , f_tt , type = 'l' , xlab = 't' , ylab = 'Hz')

title('ECG WaveForm')
abline(0,0)
points(RPeakLocations , RAmplitudes  , col = 'red')

stdresid2[stdresid2 > 0] = 0
stdresid2 <- abs(stdresid2)
stdresid2 <- stdresid2 > ( mean(stdresid2) + 2.5*sqrt( mean(stdresid2)) )

# Use a filter to create a local clique for every R peak.
Filter2 =  rep(1 , 41)
Rregion <- (abs(convolve(Rlogical , Filter2 , type = "open")) < 0.5)
Rregion <- Rregion[(length(Filter2)/2):(length(Rregion)-(length(Filter2)/2))]
stdresid2[Rregion] = 0

QSlogical <- FindLocalTurningPoints(stdresid2 , f_tt , 0)
QSlocations <- tt[QSlogical == TRUE]
QSValues <- f_tt[QSlogical == TRUE]

# Find closest Neighbours
Temp <- get.knnx(QSlocations , RPeakLocations , k = 2)
Temp$nn.index <- t(apply(Temp$nn.index, 1 , sort))
QLocations <- QSlocations[Temp$nn.index[ , 1]]
SLocations <- QSlocations[Temp$nn.index[ , 2]]

QAmplitudes <- QSValues[Temp$nn.index[ , 1]]
SAmplitudes <- QSValues[Temp$nn.index[ , 2]]
QStime <- t(apply(Temp$nn.dist , 1 , sum))

par(mfrow = c(1 , 1))
plot(tt , f_tt , type = 'l' , xlab = 't' , ylab = 'Hz')
points(QLocations , QAmplitudes  , col = 'blue')
points(SLocations , SAmplitudes  , col = 'green')
title('ECG WaveForm')
abline(0,0)
points(RPeakLocations , RAmplitudes  , col = 'red')


par(mfrow = c( 3 , 1))
plot(RPeakLocations , RAmplitudes, xlab = 't' , ylab = 'Hz')
title('R-peak Amplitudes')
plot(RPeakLocations  , c(diff(RPeakLocations) , mean(diff(RPeakLocations))) , xlab = 'Index' , ylab = 't' , ylim = c(0.6 , 1))
title('R-peak Times')
plot(RPeakLocations , QStime)
abline(mean(QStime) , 0 )
abline(mean(QStime) + 2*sqrt(var(t(QStime))) , 0  , col = 'red')
abline(mean(QStime) - 2*sqrt(var(t(QStime))) , 0 )
title('QSTime')

QRSLogical <- matrix(0 , length(tt) , 1)
bandincrement <- 10*0.005

for(i in 1:length(QLocations))
{
  QRSLogical[ ((tt >= ( QLocations[i] - bandincrement ) )*((tt <= (SLocations[i] + bandincrement)  ))) == 1 ] <- 1 +  QRSLogical[ ((tt >= (QLocations[i] - bandincrement)  )*((tt <= (SLocations[i] + bandincrement)  ))) == 1 ]
}

QRSWaveForm <- matrix(0 , length(tt) , 1)
QRSWaveForm[QRSLogical == 1] <- f_tt[QRSLogical == 1]
QRSWaveForm[QRSLogical == 0] <- NA
QRSWaveForm[1] <- 0
QRSWaveForm[length(QRSWaveForm)] <- 0
QRSWaveForm <- na.approx(QRSWaveForm)

ResidualWaveFrom <- matrix(0 , length(tt) , 1)
ResidualWaveFrom[QRSLogical == 0] <- f_tt[QRSLogical == 0]
ResidualWaveFrom[QRSLogical == 1] <- NA
ResidualWaveFrom <- na.approx(ResidualWaveFrom)



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
