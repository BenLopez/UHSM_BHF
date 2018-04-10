pathFiles <- choose.dir(caption = "Select folder with source code")
pathFiles <- paste0( pathFiles , "\\")
setwd( pathFiles )

source("LibrariesAndSettings.R" , print.eval  = TRUE )

# Load in R-R time and R-amplitude extracted data. (The script was tested on 1139 z-file)
DataPath <- choose.files(caption = 'Choose R-R time and R-amplitube .Rdata files')
load(DataPath)

PriorTimePeriod <- 4 # Choose time period in minutes to specify priors
n<- 30

SecondOrderSpecifiction <- BayesLinearDynamicUpdateAutomatedPrior(RWaveExtractedData , PriorTimePeriod , n)
AdjustedBeliefs <- BayesLinearDynamicUpdateCalulateAdjustedBeliefs(RWaveExtractedData , SecondOrderSpecifiction , n)

indexupperbound <- 1000
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

