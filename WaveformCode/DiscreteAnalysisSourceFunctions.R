ExtractHeartratefromDiscreteData <- function(DiscDataTrim)
{
  # function to extract heartrate from discretedata and sort into appropriate order.
  HeartRate <- DiscDataTrim$Value[DiscDataTrim$VarName == DiscDataTrim$VarName[3]]
  tt <- DiscDataTrim$Date[DiscDataTrim$VarName == DiscDataTrim$VarName[3]]
  HeartRate <- data.frame(tt , HeartRate)
  HeartRate <- HeartRate[order(HeartRate$tt) , ]
  output <- data.frame(HeartRate)
  return(output)
}

LoadSetofDiscreteDataHeartRates <- function(subList , path)
{
# Function to load and process a set of discrete data files.
  
counter <- 1
DiscreteDataSet <- list()
MetaData <- list()
for( PatientCode in subList)
# Check ECG1Rdata has been prcossed
{
    if(!file.exists( paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'Discrete_', PatientCode , '.RData') )){
      print(paste0(PatientCode , ' has no discrete.rData skipping to next file'))
      next
    }
    
    load(paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'Discrete_', PatientCode , '.RData'))
    sub_pat = subset(PatIndex2017, PseudoId %in% PatientCode)
    
    HeartRate <- ExtractHeartratefromDiscreteData(DiscDataTrim)
    
    TimeofNAFEvent <- as.POSIXct(strptime(sub_pat$FirstNewAF , "%Y-%m-%d %H:%M" ) )

    DiscreteDataSet[[counter]] <- HeartRate
    MetaData[[counter]] <- sub_pat
    counter <- counter + 1
}
output <- list(DiscreteDataSet , MetaData)
output <- setNames(output , c('DiscreteDataSet' , 'MetaData'))
return(output)

}

LoadSingleDiscreteDataHeartRates <- function(subList , path)
{
  # Function to load and process a set of discrete data files.
  
  PatientCode <- subList
  load(paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'Discrete_', PatientCode , '.RData'))
  sub_pat = subset(PatIndex2017, PseudoId %in% PatientCode)
    
  HeartRate <- ExtractHeartratefromDiscreteData(DiscDataTrim)
  DiscreteDataSet <- HeartRate
  MetaData <- sub_pat

output <- list(DiscreteDataSet , MetaData)
output <- setNames(output , c('Data' , 'MetaData'))
return(output)
}  


DetectAFDiscreteHeartRateApproach <- function( DiscreteDataSet , HeartRatethreshold = 130 , Window = 100 )
{
  
#Filter <- c(rep(1 , Window) , rep(0 , Window) )
AbnormalHeartRateLogical <- list()

for( i in c(1:length(DiscreteDataSet)) )
{
  if(length(DiscreteDataSet[[i]]$HeartRate$tt) < 1000)
  {
    AbnormalHeartRateLogical[[i]]<-0*DiscreteDataSet[[i]]$HeartRate$HeartRate
    next
  }
    
  AbnormalHeartRateLogical[[i]] <- cumsum(DiscreteDataSet[[i]]$HeartRate$HeartRate > HeartRatethreshold)
  AbnormalHeartRateLogical[[i]] <- AbnormalHeartRateLogical[[i]][(Window+1):length(AbnormalHeartRateLogical[[i]])] - AbnormalHeartRateLogical[[i]][1:(length(AbnormalHeartRateLogical[[i]]) - Window)] 
  #AbnormalHeartRateLogical[[i]] <- imfilter1D(AbnormalHeartRateLogical[[i]] , Filter)
}
  return(AbnormalHeartRateLogical)
}

FindNumberUniques <- function(X)
{
  # Function to find number of unique values in a column vector.
  values <- unique(X)
  n <- matrix(0 , length(values) , 1)
  nn <- matrix(0 , length(values) , 1)
  for(i in 1:length(values))
  {
    n[i]<-sum( X == values[[i]])  
    nn[i] <- values[[i]]
  }
 output <- data.frame( nn , n)
 output <- setNames(output , c('values' , 'n'))
 return(output) 
}

ExtractNumberofErrors <- function(Errorpatients)
{

Srinklist <- function(X){
  X <- X[[1]]
  return(X)
}

Totalfalsepositives <- as.matrix(Errorpatients[[1]][[1]])
Totalfalsenegatives <- as.matrix(Errorpatients[[1]][[2]])

# Count number of false positive and flase negatives 
for( i in 2:length( Errorpatients ) )
{
  Totalfalsepositives <- rbind( Totalfalsepositives , as.matrix(Errorpatients[[i]][[1]]) )
  Totalfalsenegatives <- rbind( Totalfalsenegatives , as.matrix(Errorpatients[[i]][[2]]) )
}

Totalfalsepositives <- FindNumberUniques(as.vector(lapply(Totalfalsepositives , Srinklist )))
Totalfalsenegatives <- FindNumberUniques(as.vector(lapply(Totalfalsenegatives , Srinklist )))
output <- list(Totalfalsepositives , Totalfalsenegatives)
output <- setNames(output , c('FalsePositives' , 'FalseNegatives'))
}

ComputeLocalandGlobalSecondOrderStatistics <- function(DiscreteDataSet , Window = 500  , Window2 = 100)
{

  DiscreteDataSet$HeartRate$HeartRate[is.na(DiscreteDataSet$HeartRate$HeartRate)] <- mean(DiscreteDataSet$HeartRate$HeartRate , rm.na = TRUE)
  GlobalComponent <- smth( DiscreteDataSet$HeartRate$HeartRate , Window , method = 'sma')
  LocalComponent <- smth( DiscreteDataSet$HeartRate$HeartRate - GlobalComponent , Window2, method = 'sma')
  VarLocalCompnent <- smth( (DiscreteDataSet$HeartRate$HeartRate - GlobalComponent)^2 , Window2  , method = 'sma' ) - LocalComponent^2
  output <- setNames( 
    data.frame( DiscreteDataSet$HeartRate$tt[(Window + Window2) : length(GlobalComponent)]
                , GlobalComponent[(Window + Window2) : length(GlobalComponent)] 
                , LocalComponent[(Window + Window2) : length(GlobalComponent)] 
                , VarLocalCompnent[(Window + Window2) : length(GlobalComponent)]) 
                , c( 't', 'Glo' , 'Lo' , 'V_L'))
  
  return(output)
}
