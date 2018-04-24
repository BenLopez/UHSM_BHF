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
RWaveExtractedData <- RPeakExtractionWavelet( WaveData , Filter , stdthresh = 1.5)

par(mfrow = c( 3 , 1))
regionofinterest = regionofinterest + 5000
plot(WaveData[regionofinterest  , 1] , WaveData[regionofinterest , 2], type = 'l', xlab = 't' , ylab = 'Hz')
points(RWaveExtractedData$t , RWaveExtractedData$RA  , col = 'blue', xlab = 't' , ylab = 'Hz')
title('z1127 ECG Wave Form')
plot(RWaveExtractedData$t , RWaveExtractedData$RA, xlab = 't' , ylab = 'Hz')
title('R-peak Amplitudes')
plot(RWaveExtractedData$t  , RWaveExtractedData$RR , xlab = 't' , ylab = 't' , ylim = c(0.65 , 0.8))
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
rangetotest = c(3000:5000)
QRSoutput <- FindQRSComplex(WaveData[rangetotest,] , Filter)
tt<- WaveData[rangetotest , 1]
f_tt <- WaveData[rangetotest  , 2]

par(mfrow = c(2 , 1))
plot(tt , f_tt , type = 'l' , xlab = 't' , ylab = 'Hz')
title('z 1126 ECG WaveForm')

points(QRSoutput$Qt , QRSoutput$QA  , col = 'blue')
points(QRSoutput$St , QRSoutput$SA  , col = 'green')
abline(0,0)
points(QRSoutput$Rt , QRSoutput$RA  , col = 'red')

SeparatedWaveForms <- SeparateWaveQRSandPTWaveforms( WaveData[rangetotest,] , QRSoutput ,  bandincrement = (0*0.005))

output <- waveletremovecompenentsandreconstruct(SeparatedWaveForms$PTWaveFrom , Filter  , 12 , 7)
stdresid <- output/sqrt( var(output) )
stdresid[stdresid<0] = 0
TLocations <- FindLocalTurningPoints(stdresid>mean(stdresid) , SeparatedWaveForms$PTWaveFrom , 1)
TAmplitudes <- f_tt[TLocations]
TLocations <- tt[TLocations]

output <- waveletremovecompenentsandreconstruct(SeparatedWaveForms$PTWaveFrom , Filter  , 12 , 6)
stdresid <- output/sqrt(var(output))
stdresid[stdresid<0] = 0
PTLocations <- FindLocalTurningPoints(stdresid>mean(stdresid) , SeparatedWaveForms$PTWaveFrom , 1)
PTAmplitudes <- f_tt[PTLocations]
PTLocations <- tt[PTLocations]

par(mfrow = c(1,1))
plot(tt , f_tt , type = 'l')
points(QRSoutput$Qt , QRSoutput$QA  , col = 'blue')
points(QRSoutput$St , QRSoutput$SA  , col = 'green')
title('z 1126 ECG WaveForm')
abline(0,0)
points(QRSoutput$Rt , QRSoutput$RA  , col = 'red')
points(PTLocations , PTAmplitudes , col = 'black')
points(TLocations , TAmplitudes , col = 'yellow')

par(mfrow = c(3,1))
plot(tt , output , type = 'l')
firstdiff = c(diff(output) , 0)
seconddiff <- c(diff(diff(output)) ,0, 0)
logicalforpeaks = ((abs(firstdiff) < (mean(firstdiff) + 2*sqrt(var(firstdiff)) ))*(seconddiff<0) == 1)
points(tt[logicalforpeaks] , output[logicalforpeaks] , col = 'red')
plot(tt , firstdiff , type = 'l')
points(tt[logicalforpeaks ] , firstdiff[logicalforpeaks] , col = 'red')
abline(0,0)
plot(tt , seconddiff , type = 'l')
points(tt[logicalforpeaks ] , seconddiff[logicalforpeaks] , col = 'red')

PTLocations <- FindLocalTurningPoints(logicalforpeaks , SeparatedWaveForms$PTWaveFrom , 1)
PTAmplitudes <- f_tt[PTLocations]
PTLocations <- tt[PTLocations]

par(mfrow = c(1,1))
plot(tt , f_tt , type = 'l')
points(QRSoutput$Qt , QRSoutput$QA  , col = 'blue')
points(QRSoutput$St , QRSoutput$SA  , col = 'green')
title('z 1126 ECG WaveForm')
abline(0,0)
points(QRSoutput$Rt , QRSoutput$RA  , col = 'red')
points(PTLocations , PTAmplitudes , col = 'black')
points(TLocations , TAmplitudes , col = 'yellow')

PQRSToutput <- ExtractPQRST( WaveData[1:1000000,] , Filter )
par(mfrow = c(5 , 1))
plot(PQRSToutput[['Pt']] , PQRSToutput[['PA']] , col = 'black' , xlab = 't', ylab = 'Hz' )
title('P-Amplitudes')
plot(PQRSToutput[['Qt']] , PQRSToutput[['QA']]  , col = 'blue' , xlab = 't', ylab = 'Hz')
title('Q-Amplitudes')
plot(PQRSToutput[['Rt']] , PQRSToutput[['RA']]  , col = 'red', xlab = 't', ylab = 'Hz')
title('R-Amplitudes')
plot(PQRSToutput[['St']] , PQRSToutput[['SA']]  , col = 'green', xlab = 't', ylab = 'Hz')
title('S-Amplitudes')
plot(PQRSToutput[['Tt']] , PQRSToutput[['TA']] , col = 'yellow', xlab = 't', ylab = 'Hz')
title('T-Amplitudes')
