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

for( ii in 1:length(listAllPatients) ){
outputdata <- list()
if( !DP_CheckECGreducedfilesprocessed(path , listAllPatients[ii] , Filestoprocess = 'ECGI_reduced') ){
  next
}
  for( jj in 1:3){  
  print(paste0('Extracting Rpeaks ' , jj , '.'))
  WaveData <- DP_LoadECGReduced(path , subList = listAllPatients[ii] ,numberrep , jj )
  outputdata[[jj]] <- CleanRpeaks(RPeakExtractionWavelet( WaveData , wt.filter( filter = "d6" , modwt=TRUE , level=1 ) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
  outputdata[[jj]] <- RPeakExtractionWavelet( WaveData , wt.filter( filter = "d6" , modwt=TRUE , level=1 ) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) 
  print(paste0('Rpeaks ' , jj , ' extracted.'))
  }
  outputdata[[4]] <- DP_CreateDummyMetaData(PatIndex2017  , listAllPatients[ii] )
  print(paste0('Combining Rpeaks '))
  outputdata <- setNames( outputdata , c(FilestoProcess , 'MetaData') )
  ECGs <- DP_LoadReducedECGs(path , subList = listAllPatients[ii] ,numberrep , FilestoProcess = c('ECGI' , 'ECGII' , 'ECGIII') )
  outputdata[[5]] <- PE_MultipleECGRPeaks(outputdata = outputdata , ECGs = ECGs)
  outputdata <- setNames( outputdata , c(FilestoProcess , 'MetaData' , 'RRCombined') )
  save( outputdata , file = paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' ,  listAllPatients[ii]  , '_RPeaks.RData' ) )
  print(paste0('Rpeaks Combined'))
DP_WaitBar(ii/length(listAllPatients))
}



#x11()
#par(mfrow = c(2 , 2))
#plot(outputdata[[1]]$t  , outputdata[[1]]$RR  , col = rgb(0,0,1,alpha = 0.05) , pch = 16)
#plot(outputdata[[2]]$t  , outputdata[[2]]$RR  , col = rgb(0,0,1,alpha = 0.05) , pch = 16)
#plot(outputdata[[3]]$t  , outputdata[[3]]$RR  , col = rgb(0,0,1,alpha = 0.05) , pch = 16)
#plot(outputdata[[5]]$t  , outputdata[[5]]$RR  , col = rgb(0,0,1,alpha = 0.05) , pch = 16)


