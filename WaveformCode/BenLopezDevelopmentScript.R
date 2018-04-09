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
HoursBeforeEnd = 6;
WaveData <- WaveData[ ( WaveData$Date )> (WaveData$Date[length(WaveData$Date)] - HoursBeforeEnd*(60^2))  , 1:2]

# Check first chunk of data.
regiontocheck <- c(1:500);
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
PriorRA <- as.matrix(RWaveExtractedData$RA[RWaveExtractedData$t < (RWaveExtractedData$t[1] + 60*PriorTimePeriod ) ])
PriorRR <- as.matrix(RWaveExtractedData$RR[RWaveExtractedData$t < (RWaveExtractedData$t[1] + 60*PriorTimePeriod ) ])
PriorRRt <- as.matrix(RWaveExtractedData$t[RWaveExtractedData$t < (RWaveExtractedData$t[1] + 60*PriorTimePeriod ) ])
PriorData <- cbind(PriorRRt , PriorRA , PriorRR)

N <- dim(PriorRA)
n <- 30
# producing the priors

nminusonematrix <- matrix(0 , N[1] - (n) ,  (n-1))
nplus1matrix <- matrix(0 , N[1] - (n) ,  1)

for(i in 1:(N[1] - n))
{
  
  nminusonematrix[i,] <- PriorRA[(1 + (i-1)) :((n-1)+ (i-1))]
  nplus1matrix[i , 1] <- PriorRA[(n-1) + i]  

}

E_D <- as.matrix(apply(cbind(nminusonematrix , nplus1matrix) , 1 , mean))
V_D <- cov(cbind(nminusonematrix , nplus1matrix))


E_z <- E_D[1:(n-1),]
V_Z <- V_D[1:(n-1) , 1:(n-1)]
E_x <- E_D[n,]
V_x <- V_D[n , n]
C_xz <- V_D[1:(n-1) , n]

L = t(chol(V_Z))
W = C_xz%*%solve(t(L) , solve(L))

E_D_z <- 0*matrix(0 , length(RWaveExtractedData$RA) , 1)
V_D_z <- 0*matrix(0 , length(RWaveExtractedData$RA) , 1)
stdresid <-0*matrix(0 , length(RWaveExtractedData$RA) , 1)

for(i in 1:length(RWaveExtractedData$RA))
{
E_D_z[i , 1] <- E_x +  W%*%(RWaveExtractedData$RA[(1 + (i-1)):((n-1) + (i-1))] - E_z )
V_D_z[i , 1] <- V_x - W%*%(C_xz)
stdresid[i,1] <- (E_D_z[i , 1] - RWaveExtractedData$RA[(n-1) + i])/(sqrt(V_D_z[i , 1] + 0.1))
}

plot(RWaveExtractedData$t[1:2000], stdresid[1:2000] , ylim=c(-10, 10))
abline( 3 ,0)
abline(-3 ,0)
