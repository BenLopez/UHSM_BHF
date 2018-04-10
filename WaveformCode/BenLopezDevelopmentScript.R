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


# Camila's method
# priors 
PriorTimePeriod <- 2
n <- 30
SecondOrderSpecifiction <- BayesLinearDynamicUpdateAutomatedPrior(RWaveExtractedData , PriorTimePeriod , n)

E_D_z <- matrix(0 , length(RWaveExtractedData$RA) , 2)
stdresid <-matrix(0 , length(RWaveExtractedData$RA) , 2)
discrepancyadjustedversion <-matrix(0 , length(RWaveExtractedData$RA) , 1)


for(i in 1:length(RWaveExtractedData$RA))
{
if(i < (length(RWaveExtractedData$RA) - (n-1)))
{  
diff <- as.vector(cbind( RWaveExtractedData$RA[(1 + (i-1)):((n-1) + (i-1))] , RWaveExtractedData$RR[(1 + (i-1)):((n-1) + (i-1))]) - E_z )

E_D_z[i , 1:2] <- E_x +  W%*%diff
diff <- (E_D_z[i , 1:2] - cbind(RWaveExtractedData$RA[(n-1) + i] , RWaveExtractedData$RR[(n-1) + i]))
stdresid[i,1:2] <- diff/(sqrt(cbind(V_D_z[1 , 1] , V_D_z[2 , 2])))
discrepancyadjustedversion[i,1] <-diff%*%inv_V_D_z%*%t(diff)
}
}

indexupperbound = 2000
par(mfrow = c(3 , 1))
plot(RWaveExtractedData$t[1:indexupperbound], stdresid[1:indexupperbound , 1] , ylim=c(-10, 10)  , xlab = 't' , ylab = 'Std')
title('Standardised Error R-Wave Amplitude')
abline( 3 ,0)
abline(-3 ,0)
plot(RWaveExtractedData$t[1:indexupperbound], stdresid[1:indexupperbound , 2] , ylim=c(-10, 10), xlab = 't', ylab = 'Std')
abline( 3 ,0)
abline(-3 ,0)
title('Standardised Error RR-Wave Time')
plot(RWaveExtractedData$t[1:indexupperbound], discrepancyadjustedversion[1:indexupperbound , 1] , ylim=c(0, 20) , xlab = 't' , ylab = 'Dis')
abline( 8 ,0)
title('Joint Discepancy')

SecondOrderSpecifiction <- BayesLinearDynamicUpdateAutomatedPrior(RWaveExtractedData , PriorTimePeriod , n)
AdjustedBeliefs <- BayesLinearDynamicUpdateCalulateAdjustedBeliefs(RWaveExtractedData , SecondOrderSpecifiction , n)

par(mfrow = c(3 , 1))
plot(RWaveExtractedData$t[1:indexupperbound], AdjustedBeliefs[["std_D_z"]][1:indexupperbound , 1] , ylim=c(-10, 10)  , xlab = 't' , ylab = 'Std')
title('Standardised Error R-Wave Amplitude')
abline( 3 ,0)
abline(-3 ,0)
plot(RWaveExtractedData$t[1:indexupperbound], AdjustedBeliefs[["std_D_z"]][1:indexupperbound , 2] , ylim=c(-10, 10), xlab = 't', ylab = 'Std')
abline( 3 ,0)
abline(-3 ,0)
title('Standardised Error RR-Wave Time')
plot(RWaveExtractedData$t[1:indexupperbound], AdjustedBeliefs[["dis_D_z"]][1:indexupperbound , 1] , ylim=c(0, 20) , xlab = 't' , ylab = 'Dis')
abline( 8 ,0)
title('Joint Discepancy')

