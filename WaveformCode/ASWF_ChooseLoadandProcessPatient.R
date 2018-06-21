source('ASWF_ChoosePatient.R')
if(PatientRecord$TotalITUTimeHRS > 100 ){print(paste0('Total hours over 100'))
  source('ASWF_ChoosePatient.R')}
# Plot discrete data
if( DP_checkfilesprocessed(path , subList , 'Discrete') == 1 ){
plot(1)
dev.off()
myplot <- ggplot( data.frame(x<- DataSet$Data$tt[seq(1,length(DataSet$Data$tt),3)] , y <- DataSet$Data$HeartRate[seq(1,length(DataSet$Data$tt),3)] ) , aes(x,y))  + geom_point(colour="blue", alpha=0.03) + 
  ggtitle(PatientRecord$PseudoId) +
  xlab("Time") + ylab("Heart Rate") +
  geom_hline( yintercept = 130 , linetype="dashed" , color = "red" ) + 
  geom_hline( yintercept = 100 , linetype="dashed" , color = "black" ) +
  geom_hline( yintercept = 60  , linetype="dashed" , color = "blue" )  +
  geom_vline( xintercept = as.numeric(as.POSIXct(DP_StripTime(PatientRecord$FirstNewAF[1]))) , linetype="dashed" , color = "black" ) +
  geom_vline( xintercept = as.numeric(as.POSIXct(DP_StripTime(PatientRecord$ConfirmedFirstNewAF[1]))) , linetype="dashed" , color = "grey" ) +
  geom_vline( xintercept = as.numeric(as.POSIXct(DP_StripTime(PatientRecord$LastITUEntry[1]))) , linetype="dashed" , color = "red" ) 
x11(12 , 7)
print(myplot)
}

# Ask user whether to used reduced file
if( DP_CheckECGreducedfilesprocessed( path , subList  , "ECGI_reduced") ){
  UseReduced <- winDialog(type = c('yesno') , message = 'A reduced ECG has been processed. Would you like to use this?')
}  
if( !DP_CheckECGreducedfilesprocessed( path , subList  , "ECGI_reduced") ){
  UseReduced <- "NO"
}

# Load ECGI
if(UseReduced == "YES"){
 print('Loading ECGI reduced.')
 ECGI <- DP_LoadECGReduced(path , subList , numberrep , 1 )
 print('ECGI reduced loaded.')
}  
if(UseReduced == "NO"){
  # Load wave form data 
  print('Loading ECGI.')
  if(DP_checkfilesprocessed(path , subList , 'ECGI') == 0){
    if(PatientRecord$TotalITUTimeHRS > 100 )
    {print(paste0('Total hours over 100'))}else{  
    print('No ECGI data processed.')}
    source('ASWF_ChooseLoadandProcessPatient.R')
  }else{  
  WaveData <- DP_LoadECG(path , subList , numberrep , ECGNum = 1 )}
  print('ECGI Loaded.')
  interestingtimepoint <- which.min( abs( difftime(
    as.POSIXct( as.vector(as.character(round.POSIXt(DataSet$Data$tt , units = c('hours'))))),
    as.POSIXct(select.list(as.vector(as.character(unique(round.POSIXt(
      WaveData$Date[seq(from = 1, to = length(WaveData$Date) , by =1000)] , units = c('hours')))))
      , preselect = PatientRecord$FirstNewAF[1]
      , multiple = FALSE
      , title = 'Select Interest Time Point'
      , graphics = TRUE ) ), units ='hours') ))
  
  timeindex <- which.min( abs(difftime( WaveData$Date ,  DataSet$Data$tt[interestingtimepoint[1]] , units = 'secs')) )
  HoursBeforeAndAFter <- DP_SelectHoursBeforeandAfter()
  
  ECGI <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(WaveData , timeindex , HoursBeforeAndAFter))
  
}  


# Load ECGII and ECGIII 
if(UseReduced == "NO"  || !DP_CheckECGreducedfilesprocessed( path , subList  , "ECGII_reduced")){  
print('Loading ECGII.')
  if(DP_checkfilesprocessed(path , subList , 'ECGII') == 1){
    ECGII <- DP_LoadECG( path , subList , numberrep , ECGNum = 3 )
  }
  if(DP_checkfilesprocessed(path , subList , 'ECGII') == 0){
    ECGII <- ECGI}  
  
if(DP_checkfilesprocessed(path , subList , 'ECGII') == 1 ){
  
{ 
    numberhoursbefore = 0
    numberhoursafter = abs( difftime( ECGI$Date[length(ECGI$Date)] , ECGI$Date[1]  , units ='hours' ))
    HoursBeforeAndAFter = data.frame(numberhoursbefore , numberhoursafter)
}
  
timeindex <- which.min( abs(difftime( ECGII$Date , ECGI$Date[1] , units = 'secs')) )
ECGII <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(ECGII , timeindex , HoursBeforeAndAFter))
print('ECGII loaded.')
}

if(DP_checkfilesprocessed(path , subList , 'ECGII') == 0 ){
  ECGII = ECGI
  warning('No ECGII processed.')
  }
}
if(UseReduced == "NO"  || !DP_CheckECGreducedfilesprocessed( path , subList  , "ECGIII_reduced")){  
  print('Loading ECGIII.')
  if(DP_checkfilesprocessed(path , subList , 'ECGIII') == 1){
  ECGIII <- DP_LoadECG( path , subList , numberrep , ECGNum = 3 )
  }
  if(DP_checkfilesprocessed(path , subList , 'ECGIII') == 0){
    ECGIII <- ECGI}  
  if(DP_checkfilesprocessed(path , subList , 'ECGIII') == 1 ){
    { 
      numberhoursbefore = 0
      numberhoursafter = abs( difftime( ECGI$Date[length(ECGI$Date)] , ECGI$Date[1]  , units ='hours' ))
      HoursBeforeAndAFter = data.frame(numberhoursbefore , numberhoursafter)
    }
    
    timeindex <- which.min( abs(difftime( ECGIII$Date , ECGI$Date[1] , units = 'secs')) )
    ECGIII <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(ECGIII , timeindex , HoursBeforeAndAFter))
    print('ECGIII loaded.')
  }
  
  if(DP_checkfilesprocessed(path , subList , 'ECGIII') == 0 ){
    ECGIII = ECGI
    warning('No ECGIII processed.')
  }
}
if(UseReduced == "YES" &   DP_CheckECGreducedfilesprocessed( path , subList  , "ECGII_reduced")){
  print('Loading ECGII reduced.')
  ECGII <- DP_LoadECGReduced(path , subList , numberrep , 2 )
  print('ECGII reduced loaded.')
}
if(UseReduced == "YES" &   DP_CheckECGreducedfilesprocessed( path , subList  , "ECGIII_reduced")){
  print('Loading ECGIII reduced.')
  ECGIII <- DP_LoadECGReduced(path , subList , numberrep , 3 )
  print('ECGIII reduced loaded.')
}

if( DP_checkfilesprocessed(path , subList , 'Discrete') == 1 ){  
DataSet$Data <-  DataSet$Data[ ((DataSet$Data$tt > ECGI$Date[1])*(DataSet$Data$tt < ECGI$Date[length( ECGI$Date)]) == 1), ]
}

if( UseReduced == "NO"  || !DP_checkRpeaksfilesprocessed(path , subList) ){
  print('Extracting RA and R-R times.')
  outputdata <- list()
  outputdata[[1]] <- CleanRpeaks(RPeakExtractionWavelet( ECGI   , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  outputdata[[2]] <- CleanRpeaks(RPeakExtractionWavelet( ECGII  , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  outputdata[[3]] <- CleanRpeaks(RPeakExtractionWavelet( ECGIII , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  outputdata[[4]] <- PatientRecord
  outputdata[[5]] <- 1
  outputdata <- setNames( outputdata , c('ECGI' ,'ECGII' ,'ECGIII' , 'Meta_Data' , 'RRCombined') )
  outputdata[[5]] <- PE_MultipleECGRPeaks(outputdata)
  print('RA and R-R times Extracted.')  
}
if( UseReduced == "YES" &  DP_checkRpeaksfilesprocessed(path , subList) ){
print( 'Loading Rpeaks.' )  
  outputdata <- DP_LoadRpeaksfile(path , subList) 

if(DP_ValidateRPeaks(outputdata) == FALSE){
  print('Extracting RA and R-R times.')
  outputdata <- list()
  outputdata[[1]] <- CleanRpeaks(RPeakExtractionWavelet( ECGI   , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  outputdata[[2]] <- CleanRpeaks(RPeakExtractionWavelet( ECGII  , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  outputdata[[3]] <- CleanRpeaks(RPeakExtractionWavelet( ECGIII , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  outputdata[[4]] <- PatientRecord
  outputdata[[5]] <- 1
  outputdata <- setNames( outputdata , c('ECGI' ,'ECGII' ,'ECGIII' , 'Meta_Data' , 'RRCombined') )
  outputdata[[5]] <- PE_MultipleECGRPeaks(outputdata)
  print('RA and R-R times Extracted.')
}
print( 'Rpeaks loaded.' )  
}

InferenceOutput <- AFD_DetectionWrapper( outputdata$RRCombined )
AFScore <- InferenceOutput$AFScore
if(length(InferenceOutput$StartEndTimesAF$Start) >0 ){
  StartEndTimesAF <- AFD_Checkformissingdata(InferenceOutput$StartEndTimesAF , AFScore , ECGI , ECGII , ECGIII)}else{
  StartEndTimesAF <- InferenceOutput$StartEndTimesAF
  }
if(length(InferenceOutput$StartEndTimesMM$Start) >0 ){
  StartEndTimesMM <- AFD_Checkformissingdata(InferenceOutput$StartEndTimesMM , AFScore , ECGI , ECGII , ECGIII)}else{
  StartEndTimesMM <- InferenceOutput$StartEndTimesMM
  }


timelist <- as.vector(as.character(round.POSIXt(ECGI[seq(from = 1 , to = length(ECGI[ , 1]) , by = 1000), 1] , units = 'mins')))

startindex = which.min( abs( as.POSIXct( ECGI$Date ) - as.POSIXct(select.list(unique(timelist)
                                            , preselect = NULL
                                            , multiple = FALSE
                                            , title = 'Select time to view ECGI'
                                            , graphics = TRUE ) ) ) )[1]

timeinterval <- as.numeric(select.list(unique(as.character(c(1:100)))
                                       , preselect = '10'
                                       , multiple = FALSE
                                       , title = 'Select number of seconds of full ECG to be viewed'
                                       , graphics = TRUE ))

endindex <- startindex + (round(timeinterval / as.numeric(abs(ECGI$Date[1]- ECGI$Date[2]))))
regionofinterest <- startindex:endindex

regionofinterest2  <- ASWF_AlignRegionofInterests(ECGI , ECGII , regionofinterest)
regionofinterest3  <- ASWF_AlignRegionofInterests(ECGI , ECGIII , regionofinterest)



p1 <- ggplot(outputdata$ECGI[seq(1 , length(outputdata$ECGI$t) , 3) , ] , aes(t , RA)) +
  geom_point( colour = "blue" ,  alpha=0.03 ) +
  geom_point(data=outputdata$ECGII[seq(1 , length(outputdata$ECGI$t) , 3) , ] , colour = "red"  ,  alpha=0.03 ) +
  geom_point(data= outputdata$ECGIII[seq(1 , length(outputdata$ECGI$t) , 3) , ] , colour = "black",  alpha=0.03 ) +
  ggtitle('R-amplitudes') +
  xlab("t") +
  ylab("RA") + coord_cartesian(ylim = c(0, 200)) 


p2 <- ggplot() + 
      geom_line( data = AFScore , aes(x = t , y = IHAVFScore/150) , colour ='red' , alpha = 0.25 )  + 
      geom_point(data = outputdata$RRCombined  , aes(x = t , y = RR) , colour="blue", alpha=0.01 ) +
      scale_y_continuous(sec.axis = sec_axis(~.*150, name = "AF Score")) +
      xlab("t") +
      ylab("RR") + coord_cartesian(ylim = c(0, 1.5))


if(length(StartEndTimesAF$Start) > 0)
{
for( i in ( 1:length(StartEndTimesAF$Start) ) )
{
  p2 <- p2 + annotate("rect" , xmin = StartEndTimesAF$Start[i], xmax = StartEndTimesAF$End[i], ymin = -1000, ymax= 1000 , fill = 'pink' , alpha = 0.25)
}
}

if(length(StartEndTimesMM$Start) > 0)
{
  for( i in ( 1:length(StartEndTimesMM$Start) ) )
  {
    p2 <- p2 + annotate("rect" , xmin = StartEndTimesMM$Start[i], xmax = StartEndTimesMM$End[i], ymin = -1000, ymax= 1000 , fill = 'black' , alpha = 0.25)
  }
}


#p4 <- ggplot(DataSet$Data , aes(tt , HeartRate)) +
# geom_point(colour="blue", alpha=0.1) +
#  ggtitle('Discrete Heart Rate') +
#  xlab("t") +
#  ylab("HR") 

rm(WaveData)
