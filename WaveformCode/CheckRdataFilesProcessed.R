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

listoffiles <- list.files(choose.dir())
listofpatients <- unique(lapply(listoffiles , function(X){substr(X , start = ( regexpr('_',X  )[[1]] +1) , stop = (nchar(X) - 6) )}))
  
FilesProcess <- matrix( 0,  length(listofpatients) , 3 )
rownames(FilesProcess) <- listofpatients
colnames(FilesProcess) <- c('ECGI' , 'ECGII' , 'ECGIII')
for(i in 1:length(listofpatients)){
  if(sum( paste0('ECGI_' , listofpatients[i] , '.RData') %in% listoffiles) ==1){
    FilesProcess[i,1] <- 1  
  }
  if(sum( paste0('ECGII_' , listofpatients[i] , '.RData') %in% listoffiles)==1){
    FilesProcess[i,2] <- 1  
  }
  if(sum( paste0('ECGIII_' , listofpatients[i] , '.RData') %in% listoffiles)==1){
    FilesProcess[i,3] <- 1 
  }
}

write.csv(FilesProcess , file = 'DownloadedFiles.csv')
