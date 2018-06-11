# Script to batch process wave data.

pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )


DP_LoadPatientIndex()
HoursBeforeAndAFter <- DP_SelectHoursBeforeandAfter()
FilestoProcess <- DP_ChooseECGstoProcess() 
HowtoFilterops <- read.csv(choose.files(caption = "Select listofopsSH") , stringsAsFactors = FALSE)

DP_ChooseDataReps()
listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
listAllPatients <- select.list(listAllPatients , graphics = TRUE , multiple = TRUE)

for( ii in 1:length(listAllPatients) )
{
processdata <- 1
outputdata  <- list()  
    
for( jj in 1:length(FilestoProcess) )
{
  
  print( paste0('Loading patient ' , listAllPatients[ii] , ' ' , FilestoProcess[jj] , '.'  ) )
  WaveData <- DP_LoadECG(path , listAllPatients[ii] , numberrep , FilestoProcess[jj] )
  print( paste0('Patient ' , listAllPatients[ii] , ' ' , FilestoProcess[jj]  , ' loaded.' ) )
  
  sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[ii])

print('Cropping ECG.')  
  if(!is.na(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1])))
  {
    timeindex <- which.min(abs(difftime(WaveData$Date , DP_StripTime(sub_pat$ConfirmedFirstNewAF[1]) , units = 'mins' )))[1]
    WaveData <- DP_CropWaveData(WaveData ,  timeindex , HoursBeforeAndAFter)
  }
  
  if(is.na(sub_pat$ConfirmedFirstNewAF[1]) & !is.na(DP_StripTime(sub_pat$FirstNewAF[1])))
  {
    timeindex <- which.min(abs(difftime(WaveData$Date , DP_StripTime(sub_pat$FirstNewAF), units = 'mins' )))[1]
    WaveData <- DP_CropWaveData(WaveData ,  timeindex , HoursBeforeAndAFter)
  }
  
  if(is.na(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1])) & is.na(DP_StripTime(sub_pat$FirstNewAF[1]) ))
  {
  
  if(jj == 1)
  {
    counter <- 1
    while(counter < 10)
  {
   
  if(length(WaveData$Date) < ( (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter)*(60^2)/0.005) ){
    break
    counter <- 11}    
      
   timeindex <- sample( (1 + (HoursBeforeandAfter$numberhoursbefore*(60^2)/0.005)):(length(WaveData$Date) - (HoursBeforeandAfter$numberhoursafter*(60^2)/0.005)) , 1  )
   timedifference <- abs(difftime( WaveData$Date[timeindex - (HoursBeforeandAfter$numberhoursbefore*(60^2)/0.005)]  , WaveData$Date[timeindex + (HoursBeforeandAfter$numberhoursafter*(60^2)/0.005)] , units ='hours'))
   
   if(timedifference > (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter + 1)){
   counter <- counter + 1  
   next}
   
   if(timedifference < (HoursBeforeandAfter$numberhoursbefore + HoursBeforeandAfter$numberhoursafter + 1)){
     WaveData <- DP_CropWaveData(WaveData ,  timeindex , HoursBeforeAndAFter)
     timestamp <- WaveData$Date[timeindex]
     break}
  }
  
if(counter > 10){
      warning(paste0('No good time period found for patient ') , listAllPatients[ii] ,'.' ) 
  processdata <- 0    
  break}  
  }
    
  if(jj > 1)
  {
    timeindex <- which.min(abs(WaveData$Date -  timestamp))
    WaveData <- DP_CropWaveData(WaveData ,  timeindex , HoursBeforeAndAFter)
  }
    
}  
print('ECG cropped.')    

if(processdata == 1)
{  
print('Extracting Rpeaks.')
outputdata[[jj]] <- CleanRpeaks(RPeakExtractionWavelet( WaveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
print('Rpeaks extracted.')
}

save( WaveData , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_' , FilestoProcess[jj]  , '_reduced.RData' ) )
rm(WaveData)  

}
 
if(processdata == 1)
{ 
  outputdata[[length(outputdata) + 1]] <- sub_pat 
  outputdata <- setNames( outputdata , c(FilestoProcess , 'Meta_Data') )
  print('Saving output.')
  save( outputdata , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_RPeaks.RData' ) )
  print('Output saved.')
  rm(outputdata)
}
}