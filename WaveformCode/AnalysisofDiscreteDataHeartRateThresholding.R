pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

source("LibrariesAndSettings.R" , print.eval  = TRUE )

# Load patient index
path_PatIndex <- choose.files(caption="Select 2017 PatientIndex.csv file")
PatIndex2017  <- read.csv(file=path_PatIndex, stringsAsFactors = FALSE)
path_PatIndex <- choose.files(caption="Select 2017 PatientIndex.csv file")
PatIndex2017  <- rbind(PatIndex2017 , read.csv(file=path_PatIndex, stringsAsFactors = FALSE))

# PatIndex2017$FirstITUEntry=as.POSIXct(PatIndex2017$FirstITUEntry, format="%d/%m/%Y %H:%M")
# PatIndex2017$LastITUEntry=as.POSIXct(PatIndex2017$LastITUEntry, format="%d/%m/%Y %H:%M")

# Choose data to process
path <- choose.files()
load(path)

AFlogical <- matrix(0 , length(DataSet[["MetaData"]]) , 1)
for( i in 1:length(DataSet[["MetaData"]]) )
{
  if(!is.na(DataSet[["MetaData"]][[i]]$FirstNewAF[1])){  AFlogical[i] <- 1 } 
  if( DataSet[["MetaData"]][[i]]$Pre_OperativeHeartRhythm[1] != "Sinus Rhythm" ){ AFlogical[i] <- 1 }
}

HeartRatethreshold = c(100 , 110 , 120 , 130) 
Window = 100 
IntervalThreshold = c(1 , 25 , 50 , 75 , 99)

Sensitivity <- matrix(0 , length(HeartRatethreshold) , length(IntervalThreshold))
Specificity <-  matrix(0 , length(HeartRatethreshold) , length(IntervalThreshold))
Accuracy <-  matrix(0 , length(HeartRatethreshold) , length(IntervalThreshold))
Errorpatients <- list()
counter <- 1

for(ii in 1:length(HeartRatethreshold))
{
  
output <- DetectAFDiscreteHeartRateApproach(DataSet[["DiscreteDataSet"]] , HeartRatethreshold[ii] , Window )

for(jj in 1:length(IntervalThreshold))  
{
  
AFDetection <- matrix(0 , length(DataSet[["MetaData"]]) , 1)
TimeofDetection <- matrix(NA , length(DataSet[["MetaData"]]) , 1)
for( i in 1:length(DataSet[["MetaData"]]) )
{
    if( sum(output[[i]] > IntervalThreshold[jj]) > 0 ){ 
    AFDetection[i] <- 1
    Index <- which((output[[i]] > IntervalThreshold[jj]))
    TimeofDetection[i] <-Index #DataSet[["DiscreteDataSet"]][[1]][[1]][Index[1]]
  }
}
# Calulate sensitivity, specificity and accuracy.  
  NumberAF <- sum(AFlogical)
  N <- length(AFlogical)
  NumberCorrectAFDetections <- sum(AFDetection*AFlogical)
  NumberCorrectNonAFDetections <- sum((AFDetection == 0)*(AFlogical ==0))
  NumberInCorrectAFDetections <- NumberAF - NumberCorrectAFDetections
  Sensitivity[ii,jj]  <-  NumberCorrectAFDetections/NumberAF
  Specificity[ii,jj]  <-  NumberCorrectNonAFDetections/sum(AFlogical ==0)
  Accuracy[ii,jj]     <-  (NumberCorrectAFDetections + NumberCorrectNonAFDetections)/N
    
# Extract false positive and negative patients
  FalsePositivePatients <- list()
  FalseNegativePatients <- list()
  counter1 <- 1
  counter2 <- 1
  
  for( i in 1:length(DataSet[["MetaData"]]) )
  {
    if( ((AFDetection[i] == 1)*(AFlogical[i] ==0)) == 1){
      FalsePositivePatients[[counter1]] <- DataSet[["MetaData"]][[i]]$PseudoId
      counter1 <- counter1 + 1
    }
    if( ((AFDetection[i] == 0)*(AFlogical[i] == 1)) == 1){
      FalseNegativePatients[[counter2]] <- DataSet[["MetaData"]][[i]]$PseudoId
      counter2 <- counter2 + 1
    }
  }
  
Errorpatients[[counter]] <- list( FalsePositivePatients 
                                  , FalseNegativePatients 
                                  , c(HeartRatethreshold[ii] 
                                  , IntervalThreshold[jj]) 
                                  , Sensitivity[ii,jj] 
                                  , Specificity[ii,jj]
                                  , Accuracy[ii,jj]
                                  ,TimeofDetection)
counter <- counter + 1  
}
}

## Diagnostic plots

par(mfrow = c(1 , 3))
image.plot(Sensitivity*100 , col = heat.colors(100) , axes = F , xlab = 'Heart rate' , ylab = '5 Min Threshold' , zlim = c(0,100))
axis(1 , at = seq(0 , 1 , 1/3) , labels = HeartRatethreshold)
axis(2 , at = seq(0 , 1 , 1/4) , labels = IntervalThreshold)
title('Sensitivity')

image(Specificity*100, col = heat.colors(100), axes = F, xlab = 'Heart rate' , ylab = '5 Min Threshold', zlim = c(0,100))
axis(1 , at = seq(0 , 1 , 1/3) , labels = HeartRatethreshold)
axis(2 , at = seq(0 , 1 , 1/4) , labels = IntervalThreshold)
title('Specificity')

image(Accuracy*100, col = heat.colors(100), axes = F, xlab = 'Heart rate' , ylab = '5 Min Threshold', zlim = c(0,100))
title('Accuracy')
axis(1 , at = seq(0 , 1 , 1/3) , labels = HeartRatethreshold)
axis(2 , at = seq(0 , 1 , 1/4) , labels = IntervalThreshold)


AlwaysFalsePositive <- Errorpatients[[1]][[1]]
for( i in 2:length(Errorpatients))
{
  AlwaysFalsePositive <- intersect( as.matrix(AlwaysFalsePositive) , as.matrix(Errorpatients[[i]][[1]]) )
}

IndexAlwaysFalsePositive <- matrix(0 , length(AlwaysFalsePositive) , 1)
counter <-1
for( i in 1:length(DataSet[[2]]) )
{
if(length(intersect(AlwaysFalsePositive , DataSet[[2]][[i]]$PseudoId)) > 0){
  IndexAlwaysFalsePositive[counter]<- i 
  counter <- counter + 1}

}

waveformindex <- IndexAlwaysFalsePositive[11]
par(mfrow = c(1 , 1))
plot( DataSet[[1]][[waveformindex]][[2]] , ylab = 'Heart Rate' , xlab = 'index')
abline(130,0 , col = 'red')
abline(100,0 , col = 'pink')
abline(v = TimeofDetection[waveformindex])
title(DataSet[[2]][[waveformindex]]$PseudoId)

TotalErrorRates <- ExtractNumberofErrors(Errorpatients)

par( mfrow = c(1 , 2) )
plot( TotalErrorRates$FalsePositives$n , ylab = 'Number of Errors' , xaxt='n' )
title('False Positive')
axis(1 , TotalErrorRates$FalsePositives$values , at = seq(1 , length(TotalErrorRates$FalsePositives$values) , 1 ) , las =2 ,
     cex.axis=0.5, cex.lab=0.5, cex.main=0.5, cex.sub=0.5)
plot( TotalErrorRates$FalseNegatives$n , ylab = 'Number of Errors' , xaxt='n')
axis(1 , TotalErrorRates$FalseNegatives$values , at = seq(1 , length(TotalErrorRates$FalseNegatives$values) , 1 ) , las =2 )
title('False Negative')
