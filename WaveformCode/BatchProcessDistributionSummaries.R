{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )
DP_GetDirectories()}

for(ii in 1:length(listAllPatients)){
if(DP_checkRpeaksfilesprocessed(path , listAllPatients[[ii]])){
  outputdata <- DP_LoadRpeaksfile(path , listAllPatients[ii])}else{
    next
  }

if(DP_CheckFieldExists(outputdata , 'RRCombined')){
  print('Calulating Distribution Summaries')
  DistributionSummaries <- AFD_ExtractDistributionSummaries(outputdata$RRCombined)  
  print('Distribution Summaries calculated')
  
  print('Saving Distribution Summaries')
  DP_SaveFile(DistributionSummaries , path , listAllPatients[[ii]] , Name = paste0(listAllPatients[[ii]] , '_DistributionSummaries' ) )
  print('Distribution Summaries Saved')
  
  DP_WaitBar(ii/length(listAllPatients))
  next
}
}

