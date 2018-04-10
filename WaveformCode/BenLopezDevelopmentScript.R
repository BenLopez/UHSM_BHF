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

## Wavelets


Filter = wt.filter(filter = "d20" , modwt=TRUE, level=1)
features <- attributes(Filter)

par(mfrow = c( 3 , 1))
modoutput <-  modwt( f_t[1000:2000] , Filter , 12)
modoutputattributes <- attributes(modoutput)
W <- slot(modoutput , 'W')
V <- slot(modoutput , 'V')
W <- SetElementsOfListoToZero(W , c(5:12) )
V <- SetElementsOfListoToZero(V , c(5:12) )
slot(modoutput , 'W')<- W
slot(modoutput , 'V')<- V
imodoutput <- imodwt(modoutput, fast=TRUE)
plot(t[1:1001] , f_t[1000:2000] , type = 'l')
abline(0,0)
plot(t[1:1001] , (imodoutput) , type = 'l')
abline(0,0)
plot(t[1:1001] , imodoutput - f_t[1000:2000] , type = 'l')
abline(0,0)

multiresolutionanalysis <- attributes(mra(f_t[1000:2000] , Filter ,  12  , boundary = 'periodic' , method = "modwt"))
par(mfrow = c(6 , 1))
plot(t[1:1001] , f_t[1000:2000] , type = 'l')
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
plot(t[1:1001] , f_t[1000:2000] , type = 'l',ann = FALSE)
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
