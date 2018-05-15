pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

source("LibrariesAndSettings.R" , print.eval  = TRUE )

# Load patient index
path_PatIndex <- choose.files(caption="Select 2017 PatientIndex.csv file")
PatIndex2017  <- read.csv(file=path_PatIndex, stringsAsFactors = FALSE)
path_PatIndex <- choose.files(caption="Select 2017 (2) PatientIndex.csv file")
PatIndex2017  <- rbind(PatIndex2017 , read.csv(file=path_PatIndex, stringsAsFactors = FALSE))

# Load in discrete data processed data.
path <- choose.files(caption="Select ReducedExtractedDiscreteData.RData")
load(path)

# Remove all patients with under 5 hours of data.
DataSet <- DataSet[(as.matrix(lapply(DataSet , function(X){length(X[[1]]$tt)})) > 5000)]
#remove all patients with pre_op AF
DataSet <- DataSet[(as.matrix(lapply(DataSet , function(X){X$MetaData$Pre_OperativeHeartRhythm[[1]] == "Sinus Rhythm"})) == TRUE )]


AFlogical <- matrix(0 , length(DataSet) , 1)
PreOpSRLogical <- matrix(0 , length(DataSet) , 1)
for( i in 1:length(DataSet) )
{
  if(!is.na(DataSet[[i]][["MetaData"]]$FirstNewAF[1])){  AFlogical[i] <- 1 } 
  if( DataSet[[i]][["MetaData"]]$Pre_OperativeHeartRhythm[1] == "Sinus Rhythm" ){ PreOpSRLogical[i] <- 1 }
}

HeartRatethreshold = c(100 , 110 , 120 , 130) 
Window = 1000 
IntervalThreshold = seq(25, 999 , 50)
Sensitivity <- matrix(0 , length(HeartRatethreshold) , length(IntervalThreshold))
Specificity <-  matrix(0 , length(HeartRatethreshold) , length(IntervalThreshold))
Accuracy <-  matrix(0 , length(HeartRatethreshold) , length(IntervalThreshold))
Errorpatients <- list()
counter <- 1

for(ii in 1:length(HeartRatethreshold))
{
  
output <- DetectAFDiscreteHeartRateApproach(DataSet , HeartRatethreshold[ii] , Window )

for(jj in 1:length(IntervalThreshold))  
{
  
AFDetection <- matrix(0 , length(DataSet) , 1)
TimeofDetection <- list()
DetectionDifference <- matrix(NA , length(DataSet) , 1)

for( i in 1:length(DataSet) )
{
    if( sum(output[[i]] > IntervalThreshold[jj]) > 0 ){ 
    AFDetection[i] <- 1
    if( DataSet[[i]][["MetaData"]]$Pre_OperativeHeartRhythm[1] != "Sinus Rhythm" ){ AFDetection[i] <- 0 }
    Index <- which((output[[i]] > IntervalThreshold[jj]))
    TimeofDetection[[i]] <- DataSet[[i]]$HeartRate$tt[Index[[1]]]
    if(AFlogical[i] == 1){
      DetectionDifference[i] <- difftime(DataSet[[i]]$HeartRate$tt[Index[[1]]] , as.POSIXct(DataSet[[i]]$MetaData$FirstNewAF[[1]]) , units = 'hours')}
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
  
  for( i in 1:length(DataSet) )
  {
    if( ((AFDetection[i] == 1)*(AFlogical[i] ==0)) == 1){
      FalsePositivePatients[[counter1]] <- names(DataSet[i])
      counter1 <- counter1 + 1
    }
    if( ((AFDetection[i] == 0)*(AFlogical[i] == 1)) == 1){
      FalseNegativePatients[[counter2]] <- names(DataSet[i])
      counter2 <- counter2 + 1
    }
  }
  
Errorpatients[[counter]] <- setNames(list( FalsePositivePatients 
                                  , FalseNegativePatients 
                                  , c(HeartRatethreshold[ii] , IntervalThreshold[jj]) 
                                  , Sensitivity[ii,jj] 
                                  , Specificity[ii,jj]
                                  , Accuracy[ii,jj]
                                  , TimeofDetection
                                  , DetectionDifference) , 
                                  c('FalsePositivePatients' 
                                    ,'FalseNegativePatients'
                                    , 'Parameters' 
                                    , 'Sensitivity' 
                                    , 'Specificity' 
                                    , 'Accuracy' 
                                    , 'TimeofDetection '
                                    , 'DetectionTimeDifference')
                                  )
counter <- counter + 1  
}
}

## Diagnostic plots

par(mfrow = c(1 , 3))
image.plot(Sensitivity*100 , col = heat.colors(100) , axes = F , xlab = 'Heart rate' , ylab = '5 Min Threshold' , zlim = c(0,100))
axis(1 , at = seq(0 , 1 , 1/3) , labels = HeartRatethreshold)
axis(2 , at = seq(0 , 1 , 1/(length(IntervalThreshold) - 1) ) , labels = IntervalThreshold)
title('Sensitivity')

image(Specificity*100, col = heat.colors(100), axes = F, xlab = 'Heart rate' , ylab = '5 Min Threshold', zlim = c(0,100))
axis(1 , at = seq(0 , 1 , 1/3) , labels = HeartRatethreshold)
axis(2 , at = seq(0 , 1 , 1/(length(IntervalThreshold) - 1) ) , labels = IntervalThreshold)
title('Specificity')

image(Accuracy*100, col = heat.colors(100), axes = F, xlab = 'Heart rate' , ylab = '5 Min Threshold', zlim = c(0,100))
title('Accuracy')
axis(1 , at = seq(0 , 1 , 1/3) , labels = HeartRatethreshold)
axis(2 , at = seq(0 , 1 , 1/(length(IntervalThreshold) - 1) ) , labels = IntervalThreshold)


AlwaysFalsePositive <- Errorpatients[[1]][[1]]
for( i in 2:length(Errorpatients))
{
  AlwaysFalsePositive <- intersect( as.matrix(AlwaysFalsePositive) , as.matrix(Errorpatients[[i]][[1]]) )
}

AlwaysFalseNegative <- Errorpatients[[1]][[2]]
for( i in 2:length(Errorpatients))
{
  AlwaysFalseNegative <- intersect( as.matrix(AlwaysFalseNegative) , as.matrix(Errorpatients[[i]][[2]]) )
}



waveformindex <- which(names(DataSet)==AlwaysFalsePositive[4])
par(mfrow = c(1 , 1))
df<- data.frame(x<- DataSet[[waveformindex]]$HeartRate$tt , y <- DataSet[[waveformindex]]$HeartRate$HeartRate )
ggplot( df , aes(x,y)) + 
  geom_point(colour="blue", alpha=0.009) + 
  ggtitle(names(DataSet[waveformindex])) +
  xlab("Time") + ylab("Heart Rate") +
  geom_hline( yintercept = 130 , linetype="dashed" , color = "red" ) + 
  geom_hline( yintercept = 100 , linetype="dashed" , color = "black" ) +
  geom_hline( yintercept = 60  , linetype="dashed" , color = "blue" )  +
  geom_vline( xintercept = as.numeric(TimeofDetection[[waveformindex]]) , linetype="dashed" , color = "black" ) 


waveformindex <-  which(names(DataSet)==AlwaysFalseNegative[2])
par(mfrow = c(1 , 1))
df <- data.frame(x <- as.POSIXct(DataSet[[waveformindex]]$HeartRate$tt) , y <- DataSet[[waveformindex]]$HeartRate$HeartRate )
ggplot(df , aes(x,y)) + 
  geom_point(colour="blue", alpha=0.009) + 
  ggtitle(names(DataSet[waveformindex])) +
  xlab("Time") + ylab("Heart Rate") +
  geom_hline( yintercept = 130 , linetype="dashed" , color = "red" ) + 
  geom_hline( yintercept = 100 , linetype="dashed" , color = "black" ) +
  geom_hline( yintercept = 60  , linetype="dashed" , color = "blue" ) + 
  geom_vline( xintercept = as.numeric(as.POSIXct(DataSet[[waveformindex]]$MetaData$FirstNewAF , tz = 'GMT')) , linetype="dashed" , color = "black" ) 

  
TotalErrorRates <- ExtractNumberofErrors(Errorpatients)

par( mfrow = c(1 , 2) )
plot( TotalErrorRates$FalsePositives$n , ylab = 'Number of Errors' , xaxt='n' )
title('False Positive')
plot( TotalErrorRates$FalseNegatives$n , ylab = 'Number of Errors' , xaxt='n')
title('False Negative')

par(mfrow = c(1 , 1))
plot(Errorpatients[[1]][[4]] + 0*Errorpatients[[1]][[8]] , Errorpatients[[1]][[8]] , xlab = 'Sensitivity' , ylab = 'Time Diff' , xlim = c(0,1) , ylim = c(-100 , 100))
for( i in 2:length(Errorpatients) )
{
points(Errorpatients[[i]][[4]] + 0*Errorpatients[[i]][[8]] ,   Errorpatients[[i]][[8]] )
}
title('Forecast Potential')  

#save(Errorpatients , file = 'C:\\Users\\Ben\\Documents\\BHF Results\\DiscreteDataHeartRateThresholdResults.RData')
