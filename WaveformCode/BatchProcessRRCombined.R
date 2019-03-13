{
  {
    if(file.exists('CheckforDefaultsScript.R')){
      source('CheckforDefaultsScript.R')
    }else{
      pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
      source("LibrariesAndSettings.R" , print.eval  = TRUE )
      DP_LoadPatientIndex()
      DP_ChooseDataReps()
      FilestoProcess <- DP_ChooseECGstoProcess() 
      HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
    }
    listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
    set.seed( 1 )
  }
}  

for( ii in 1:length(listAllPatients)){
outputdata <- list()
  for( jj in 1:3){  
  print(paste0('Extracting Rpeaks ' , jj , '.'))
  WaveData <- DP_LoadECGReduced(path , subList = listAllPatients[ii] ,numberrep , jj )
  outputdata[[jj]] <- CleanRpeaks(RPeakExtractionWavelet( WaveData , wt.filter( filter = "d6" , modwt=TRUE , level=1 ) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  print(paste0('Rpeaks ' , jj , ' extracted.'))
  }
  outputdata[[4]] <- DP_CreateDummyMetaData(PatIndex2017  , listAllPatients[ii] )
  print(paste0('Combining Rpeaks '))
  outputdata <- setNames( outputdata , c(FilestoProcess , 'MetaData') )
  outputdata[[5]] <- PE_MultipleECGRPeaks(outputdata)
  outputdata <- setNames( outputdata , c(FilestoProcess , 'MetaData' , 'RRCombined') )
  save( outputdata , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_RPeaks.RData' ) )
  print(paste0('Rpeaks Combined'))
DP_WaitBar(ii/length(listAllPatients))
}

