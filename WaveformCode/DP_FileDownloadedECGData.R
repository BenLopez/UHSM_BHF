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

startpath <- choose.dir('Choose ECGs to file.')
listoffiles <- list.files(startpath)
listofpatients <- lapply(listoffiles , function(X){substr(X , start = ( regexpr('_',X  )[[1]] +1) , stop = (nchar(X) - 6) )})
Targetpath <- choose.dir(caption = 'Choose location of waveform files.')

for(i in 1:length(listoffiles)){
dir.create( paste0(Targetpath , '\\' , listofpatients[[i]] , '\\Zip_out' ) )
file.copy( from = paste0(startpath , '\\',listoffiles[i] ) ,
           to = paste0(Targetpath , '\\' , listofpatients[[i]] , '\\Zip_out' )  , recursive = TRUE )
DP_WaitBar(i/length(listoffiles))
}
