pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

source("LibrariesAndSettings.R" , print.eval  = TRUE )

# Load patient index
path_PatIndex <- choose.files(caption="Select 2017 PatientIndex.csv file")
PatIndex2017  <- read.csv(file=path_PatIndex, stringsAsFactors = FALSE)
path_PatIndex <- choose.files(caption="Select 2017 (2) PatientIndex.csv file")
PatIndex2017  <- rbind(PatIndex2017 , read.csv(file=path_PatIndex, stringsAsFactors = FALSE))

path <- choose.files(caption="choose ReducedExtractedDiscreteData.RData")
load( path )
path <- choose.files(caption="choose GroupofNonAFPatients.RData")
load( path )

# Remove data with under 100
DataSet <- DataSet[(as.matrix(lapply(DataSet , function(X){length(X[[1]]$tt)})) > 5000)]
DataSet <- DataSet[(as.matrix(lapply(DataSet , function(X){X$MetaData$Pre_OperativeHeartRhythm[[1]] == "Sinus Rhythm"})) == TRUE )]

AFlogical <- matrix(0 , length(DataSet) , 1)
PreOpSRLogical <- matrix(0 , length(DataSet) , 1)
for( i in 1:length(DataSet) )
{
  if(!is.na(DataSet[[i]][["MetaData"]]$FirstNewAF[1])){  AFlogical[i] <- 1 } 
  if( DataSet[[i]][["MetaData"]]$Pre_OperativeHeartRhythm[1] == "Sinus Rhythm" ){ PreOpSRLogical[i] <- 1 }
}

# Process Data
LocalAndGlobalComponents <- setNames(lapply( DataSet , ComputeLocalandGlobalSecondOrderStatistics ) , names(DataSet)) 


tslengths <- matrix(0 , length(GoodGroup) , 1)
for(i in 1:length(GoodGroup))
{
  tslengths[i] <- length(DataSet[[GoodGroup[i]]]$HeartRate$tt)
}

SecondOrderStatisticsMatrix <- matrix( 0 , max(tslengths) , length(GoodGroup) )

for( i in 1:length(GoodGroup) )
{
  SecondOrderStatisticsMatrix[1:length( LocalAndGlobalComponents[[GoodGroup[i]]]$Glo),i] <-  LocalAndGlobalComponents[[GoodGroup[i]]]$Glo
}

E_Glo <- apply( SecondOrderStatisticsMatrix , 1 , sum ) / apply( SecondOrderStatisticsMatrix != 0 , 1 , sum )
V_GLo <- ((apply( SecondOrderStatisticsMatrix^2 , 1 , sum ) / apply( SecondOrderStatisticsMatrix != 0 , 1 , sum ) - E_Glo^2))

SecondOrderStatisticsMatrix <- matrix( 0 , max(tslengths) , length(GoodGroup) )

for( i in 1:length(GoodGroup) )
{
  SecondOrderStatisticsMatrix[1:length( LocalAndGlobalComponents[[GoodGroup[i]]]$Lo),i] <-  LocalAndGlobalComponents[[GoodGroup[i]]]$Lo
}

E_Lo <- apply( SecondOrderStatisticsMatrix , 1 , sum ) / apply( SecondOrderStatisticsMatrix != 0 , 1 , sum )
V_Lo <- (apply( SecondOrderStatisticsMatrix^2 , 1 , sum ) / apply( SecondOrderStatisticsMatrix != 0 , 1 , sum ) - E_Lo^2)

SecondOrderStatisticsMatrix <- matrix( 0 , max(tslengths) , length(GoodGroup) )

for( i in 1:length(GoodGroup) )
{
  SecondOrderStatisticsMatrix[1:length( LocalAndGlobalComponents[[GoodGroup[i]]]$V_L),i] <-  LocalAndGlobalComponents[[GoodGroup[i]]]$V_L
}

E_V_Lo <- apply( SecondOrderStatisticsMatrix , 1 , sum ) / apply( SecondOrderStatisticsMatrix != 0 , 1 , sum )
V_V_Lo <- (apply( SecondOrderStatisticsMatrix^2 , 1 , sum ) / apply( SecondOrderStatisticsMatrix != 0 , 1 , sum ) - E_V_Lo^2)


stdresid <- list()
for(i in 1:length(DataSet))
{
  stdresid2 <- (E_Lo[1:length(LocalAndGlobalComponents[[i]]$Lo)] - LocalAndGlobalComponents[[i]]$Lo) / sqrt((V_Lo[1:length(LocalAndGlobalComponents[[i]]$Lo)]))
  stdresid3 <- (E_V_Lo[1:length(LocalAndGlobalComponents[[i]]$V_L)] - LocalAndGlobalComponents[[i]]$V_L) / sqrt((V_V_Lo[1:length(LocalAndGlobalComponents[[i]]$V_L)]))
  #stdresid[[i]] <- abs(stdresid3)
  stdresid[[i]] <- apply( abs(rbind(stdresid2 , stdresid3)) , 2 , mean ) 
}

stdresid <- setNames(stdresid , names(DataSet))

stdthresholdthreshold = c(2 , 3 , 4 , 5) 
Window = 1000 
IntervalThreshold = seq(25, 999 , 50)
Sensitivity <- matrix(0 , length(stdthresholdthreshold) , length(IntervalThreshold))
Specificity <-  matrix(0 , length(stdthresholdthreshold) , length(IntervalThreshold))
Accuracy <-  matrix(0 , length(stdthresholdthreshold) , length(IntervalThreshold))
Errorpatients <- list()
counter <- 1

for(ii in 1:length(stdthresholdthreshold))
{
  
  output <- lapply( stdresid , function(X){ cumsum(X > stdthresholdthreshold[ii]) } )
  output <- lapply(output , function(X){X[(Window+1):length(X)] - X[1:( length(X) - Window )] } )
  
  for(jj in 1:length(IntervalThreshold))  
  {
    
    AFDetection <- matrix(0 , length(DataSet) , 1)
    TimeofDetection <- list()
    DetectionDifference <- matrix(NA , length(DataSet) , 1)
    
    
    for( i in 1:length(DataSet) )
    {
      if( sum( output[[i]] > IntervalThreshold[jj],  na.rm = TRUE) > 0 ){ 
        AFDetection[i] <- 1
        if( DataSet[[i]][["MetaData"]]$Pre_OperativeHeartRhythm[1] != "Sinus Rhythm" ){ AFDetection[i] <- 0 }
        Index <- which((output[[i]] > IntervalThreshold[jj]))
        TimeofDetection[[i]] <- LocalAndGlobalComponents[[i]]$t[Index[[1]]]
        if(AFlogical[i] == 1){
          DetectionDifference[i] <- difftime(TimeofDetection[[i]] , as.POSIXct(DataSet[[i]]$MetaData$FirstNewAF[[1]]) , units = 'hours')}
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
                                               , c(stdthresholdthreshold[ii] , IntervalThreshold[jj]) 
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


par(mfrow = c(1 , 3))
image.plot(Sensitivity*100 , col = heat.colors(100) , axes = F , xlab = 'Heart rate' , ylab = '5 Min Threshold' , zlim = c(0,100))
axis(1 , at = seq(0 , 1 , 1/3) , labels = stdthresholdthreshold)
axis(2 , at = seq(0 , 1 , 1/(length(IntervalThreshold) - 1) ) , labels = IntervalThreshold)
title('Sensitivity')

image(Specificity*100, col = heat.colors(100), axes = F , xlab = 'Heart rate' , ylab = '5 Min Threshold', zlim = c(0,100))
axis(1 , at = seq(0 , 1 , 1/3) , labels = stdthresholdthreshold)
axis(2 , at = seq(0 , 1 , 1/(length(IntervalThreshold) - 1) ) , labels = IntervalThreshold)
title('Specificity')

image(Accuracy*100, col = heat.colors(100), axes = F , xlab = 'Heart rate' , ylab = '5 Min Threshold', zlim = c(0,100))
title('Accuracy')
axis(1 , at = seq(0 , 1 , 1/3) , labels = stdthresholdthreshold)
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

waveformindex <- which(names(DataSet)==AlwaysFalsePositive[10])
par(mfrow = c(1 , 1))
df<- data.frame(x<- DataSet[[waveformindex]]$HeartRate$tt , y <- DataSet[[waveformindex]]$HeartRate$HeartRate )
ggplot( df , aes(x,y)) + 
  geom_point(colour="blue", alpha=0.009) + 
  ggtitle(names(DataSet[waveformindex])) +
  xlab("Time") + ylab("Heart Rate") +
  geom_hline( yintercept = 130 , linetype="dashed" , color = "red" ) + 
  geom_hline( yintercept = 100 , linetype="dashed" , color = "black" ) +
  geom_hline( yintercept = 60  , linetype="dashed" , color = "blue" )  +
  geom_vline( xintercept = as.numeric(Errorpatients[[7]]$TimeofDetection[[waveformindex]]) , linetype="dashed" , color = "black" ) 

waveformindex <-  which(names(DataSet)== AlwaysFalseNegative[1])
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
plot( Errorpatients[[7]]$DetectionTimeDifference[!is.na( Errorpatients[[7]]$DetectionTimeDifference)] , xlab = 'Sensitivity' , ylab = 'Time Diff Hours' , ylim = c(-100 , 100))
abline(mean(Errorpatients[[7]]$DetectionTimeDifference[!is.na( Errorpatients[[7]]$DetectionTimeDifference)]) , 0)
title('Forecast Potential')  

save(Errorpatients , file = 'C:\\Users\\Ben\\Documents\\BHF Results\\DiscreteDataHeartRateStatisticalResults.RData')
