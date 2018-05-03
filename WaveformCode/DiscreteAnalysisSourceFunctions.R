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

DetectAFDiscreteHeartRateApproach <- function( DiscreteDataSet , HeartRatethreshold = 130 , Window = 100 )
{
  
Filter <- c(rep(1 , Window) , rep(0 , Window) )
AbnormalHeartRateLogical <- list()

for( i in c(1:length(DiscreteDataSet)) )
{
  AbnormalHeartRateLogical[[i]] <- DiscreteDataSet[[i]][[2]] > HeartRatethreshold
  AbnormalHeartRateLogical[[i]] <- imfilter1D(AbnormalHeartRateLogical[[i]] , Filter)
}

  return(AbnormalHeartRateLogical)
}
