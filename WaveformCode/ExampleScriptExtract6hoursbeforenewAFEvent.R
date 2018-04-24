pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# Load in patientindex data
##### Select file PatientIndex.csv if it exists
path_PatIndex = choose.files(caption="Select 2017 PatientIndex.csv file")

if(length(path_PatIndex)>0){
  PatIndex2017 = read.csv(file=path_PatIndex, stringsAsFactors = FALSE)
  # PatIndex2017$FirstITUEntry=as.POSIXct(PatIndex2017$FirstITUEntry, format="%d/%m/%Y %H:%M")
  # PatIndex2017$LastITUEntry=as.POSIXct(PatIndex2017$LastITUEntry, format="%d/%m/%Y %H:%M")
} else {
  warning("No Patient Info provided")
  sub_pat = list()
}

# Load in ECG data.
FileLocation <- choose.files(caption="Select .Rdata file of cleaned ECG data")
load(FileLocation)

# get patient pusedo ID
Temp <- gregexpr('Zip_out' , FileLocation)
Temp2 <- gregexpr('.RData' , FileLocation)
Filename <- substr(FileLocation , Temp[[1]][1] + nchar('Zip_out') +1 , Temp2[[1]][1] -1 )
Temp2 <- gregexpr('_' , Filename)
PatientCode = substr(Filename , Temp2[[1]][1] +1 , nchar(Filename)  )
rm(Temp,Temp2)

sub_pat = subset(PatIndex2017, PseudoId %in% PatientCode)
HoursBeforeEnd <- 5

# Extract last 6 hours before newAF event.
if(is.na(sub_pat$FirstNewAF)  )
{
  TimeofNAFEvent<-0
  index = round(length(WaveData$Date)/2)
  WaveData <- WaveData[  max(( index - ( HoursBeforeEnd*( (60^2) /0.005 ) ) ) , 1) : min(( index + (( (60^2) /0.005 ) ) )  , length(WaveData$Date)),  ]
}

if(is.na(sub_pat$FirstNewAF) == FALSE  )
{
  TimeofNAFEvent <- as.POSIXct(strptime(sub_pat$FirstNewAF , "%d/%m/%Y %H:%M" ) )
  index <- which.min( abs(WaveData$Date - TimeofNAFEvent) )
  WaveData <- WaveData[  max(( index - ( HoursBeforeEnd*( (60^2) /0.005 ) ) ) , 1) : min(( index + (( (60^2) /0.005 ) ) )  , length(WaveData$Date)),  ]
}  
  

# Extact peak infomation.

Filter <-  wt.filter(filter = "d6" , modwt=TRUE, level=1)
WaveData <- ReturnWaveformwithPositiveOrientation(WaveData)
RWaveExtractedData <- RPeakExtractionWavelet( WaveData , Filter , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2)

par(mfrow = c( 3 , 1))
regionofinterest = regionofinterest #c( max((length(WaveData[,1]) -  ((60^2)/0.005)  - 5000) , 1) :( length(WaveData[,1])  -  ((60^2)/0.005) ) )
if(regionofinterest[1]<0){error('Negative region of interst indicies')}
plot(WaveData[regionofinterest  , 1] , WaveData[regionofinterest , 2], type = 'l', xlab = 't' , ylab = 'Hz')
points(RWaveExtractedData$t , RWaveExtractedData$RA  , col = 'blue', xlab = 't' , ylab = 'Hz')
title(paste0(PatientCode , ' ECG Wave Form'))
plot(RWaveExtractedData$t , RWaveExtractedData$RA, xlab = 't' , ylab = 'Hz')
title('R-peak Amplitudes')
abline(v = TimeofNAFEvent , col = 'red')
abline(v = WaveData[regionofinterest[1]  , 1] , col = 'blue')
plot(RWaveExtractedData$t  , RWaveExtractedData$RR , xlab = 't' , ylab = 't' , ylim = c(0.4 , 1.5))
title('RR-peak Times')
abline(v = TimeofNAFEvent , col = 'red')
abline(v = WaveData[regionofinterest[1]  , 1] , col = 'blue')
lines( RWaveExtractedData$t , (medfilt1(RWaveExtractedData$RR , 6/0.005)) , col ='red' )


# Save file
Temp <-gregexpr('Zip_out' , FileLocation)
SaveLocation <- substr(FileLocation , 1 , Temp[[1]][1] + nchar('Zip_out'))
Temp <- gregexpr('Zip_out' , FileLocation)
Temp2 <- gregexpr('.RData' , FileLocation)
SaveLocation <- paste0(SaveLocation , substr(FileLocation , Temp[[1]][1] + nchar('Zip_out') +1 , Temp2[[1]][1] -1 ) ,  '_RPeaks.RData' )
rm(Temp,Temp2)
save('RWaveExtractedData' , file = SaveLocation)


PQRSToutput <- ExtractPQRST( WaveData[1:100000,] , Filter )
par(mfrow = c(5 , 1))
plot(PQRSToutput[['Pt']] , PQRSToutput[['PA']] , col = 'black' , xlab = 't', ylab = 'Hz' )
title('P-Amplitudes')
abline(v = TimeofNAFEvent , col = 'red')
plot(PQRSToutput[['Qt']] , PQRSToutput[['QA']]  , col = 'blue' , xlab = 't', ylab = 'Hz')
title('Q-Amplitudes')
abline(v = TimeofNAFEvent , col = 'red')
plot(PQRSToutput[['Rt']] , PQRSToutput[['RA']]  , col = 'red', xlab = 't', ylab = 'Hz')
title('R-Amplitudes')
abline(v = TimeofNAFEvent , col = 'red')
plot(PQRSToutput[['St']] , PQRSToutput[['SA']]  , col = 'green', xlab = 't', ylab = 'Hz')
title('S-Amplitudes')
abline(v = TimeofNAFEvent , col = 'red')
plot(PQRSToutput[['Tt']] , PQRSToutput[['TA']] , col = 'yellow', xlab = 't', ylab = 'Hz')
title('T-Amplitudes')
abline(v = TimeofNAFEvent , col = 'red')


