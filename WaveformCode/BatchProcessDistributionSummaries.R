{if(file.exists('CheckforDefaultsScript.R')){
  source('CheckforDefaultsScript.R')
}else{
  pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
  source("LibrariesAndSettings.R" , print.eval  = TRUE )
  DP_LoadPatientIndex()
  DP_ChooseDataReps()
  FilestoProcess <- DP_ChooseECGstoProcess() 
  HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
}
}

SettingsAFDetection <- AFD_CreateDefaultSettings()
AFPatientsOnly <- 0

for(ii in 1:length(listAllPatients)){
  
sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[[ii]])

if(nrow(sub_pat) == 0 ){
  DP_WaitBar(ii/length(listAllPatients))
  next} 
if(AFPatientsOnly == 1){
if( is.na(sub_pat$ConfirmedFirstNewAF) || sub_pat$ConfirmedFirstNewAF == 'CNAF'){
  print(paste0('Skipping Patient ' , listAllPatients[[ii]]))
  next}
}
if(DP_checkRpeaksfilesprocessed(path , listAllPatients[[ii]])){
  outputdata <- DP_LoadRpeaksfile(path , listAllPatients[ii])
  }else{
    next
}
if( length(outputdata$RRCombined$t) < 1000){
  next}
  
if(DP_CheckFieldExists(outputdata , 'RRCombined')){
  print('Calulating Distribution Summaries')
  DistributionSummaries <- AFD_ExtractDistributionSummaries(RRStruct = outputdata$RRCombined , MetaData = outputdata$Meta_Data)  
  print('Distribution Summaries calculated')
  
  print('Saving Distribution Summari}es')
  DP_SaveFile(DistributionSummaries , path , listAllPatients[[ii]] , Name = paste0(listAllPatients[[ii]] , '_DistributionSummaries' ) )
  print('Distribution Summaries Saved')
  
  DP_WaitBar(ii/length(listAllPatients))
  next
}
}
