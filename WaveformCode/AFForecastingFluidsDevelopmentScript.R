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
  #listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}



{DataStructure <- matrix(0 , length(FluidsIndex2017) , 1)
AFLogical <- matrix(0 , length(FluidsIndex2017) , 1)
PreOPAF <-matrix(0 , length(FluidsIndex2017) , 1)
setofops <- unique(PatIndex2017$ProcDetails)[c(1 , 2 ,3 , 6 , 9 ,10, 19,20, 32 , 45)]
for(i in 1:length(FluidsIndex2017)){
  
  PatientID <- DP_ExtractPatientIDFromNewPatinetID(names(FluidsIndex2017)[i])
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientCode = PatientID)

  if(is.na(MetaData$Pre_OperativeHeartRhythm[1])){MetaData$Pre_OperativeHeartRhythm[1] = "Sinus Rhythm"}
  if(MetaData$Pre_OperativeHeartRhythm[1] == "Atrial fibrillation/flutter"){
    PreOPAF[i] = T  
  }else{
    PreOPAF[i] =F
  }
  if(sum(MetaData$ProcDetails[1] %in% setofops) == 0  ){PreOPAF[i] =T}
  
  if(is.na(MetaData$ConfirmedFirstNewAF[1])){MetaData$ConfirmedFirstNewAF[1] = 'NAF'}
  if(!is.na(MetaData$FirstNewAF[1]) & (MetaData$ConfirmedFirstNewAF[1] != 'CNAF')){
    diff <- FluidsIndex2017[[i]]$TimeSeriesData$Hourly.In - FluidsIndex2017[[i]]$TimeSeriesData$Hourly.Out
    #diff <- FluidsIndex2017[[i]]$TimeSeriesData$Hourly.Out
    timevec <- FluidsIndex2017[[i]]$TimeSeriesData$time
    timevec <- timevec[!is.na(diff)]
    diff <- diff[!is.na(diff)]
    diff <- diff[!is.na(timevec)]
    timevec <- timevec[!is.na(timevec)] 
    
    if(sum(timevec < DP_StripTime(MetaData$FirstNewAF[1])) ==0){
      DataStructure[i ,1] <- diff[which(!is.na(diff))[1] ]
      PreOPAF[i] = T
      AFLogical[i] <- TRUE
      }else{
    #DataStructure[i ,1] <- diff[which(timevec < DP_StripTime(MetaData$FirstNewAF))[length(which(timevec < DP_StripTime(MetaData$FirstNewAF)))]]
    DataStructure[i ,1] <- sqrt(var(diff[which(timevec < DP_StripTime(MetaData$FirstNewAF))[(length(which(timevec < DP_StripTime(MetaData$FirstNewAF))) - min(12 , length(which(timevec < DP_StripTime(MetaData$FirstNewAF))))):length(which(timevec < DP_StripTime(MetaData$FirstNewAF)))]] , na.rm = T))
    #DataStructure[i ,1] <- mean(diff[24:min(48 , length(diff))] , na.rm = T)
    AFLogical[i] <- TRUE}
      }else{
        
    diff <- FluidsIndex2017[[i]]$TimeSeriesData$Hourly.In - FluidsIndex2017[[i]]$TimeSeriesData$Hourly.Out
    #diff <- FluidsIndex2017[[i]]$TimeSeriesData$Hourly.Out
    if(sum(!is.na(diff)) < 32){
      DataStructure[i ,1] = 0
      PreOPAF[i] = T
    }
    else{
    #DataStructure[i ,1] <- sample(diff[!is.na(diff)] , 1)
    diff <- diff[!is.na(diff)]
    sampleindex <- sample(24:(length(diff) -7) , 1)
    DataStructure[i ,1] <- sqrt(var(diff[(sampleindex -3):(sampleindex +3)] , na.rm = TRUE))
    #DataStructure[i ,1] <- mean(diff[24:min(48 , length(diff))] , na.rm = T)
    
    AFLogical[i] <- FALSE}
  }
}
}
AFLogical <- AFLogical[PreOPAF == 0]
DataStructure <- DataStructure[PreOPAF == 0 , ]


x11()
BC_PlotCompareSingleHists(DataStructure[ ((AFLogical ==0)) ==1  ],DataStructure[((AFLogical ==1)) ==1  ] , breaks = 15 , main = 'Fluid Balance 12 hours Before AFib' , xlab = 'Variance Fluid Balance')



PatientID <- DP_ExtractPatientIDFromNewPatinetID(names(FluidsIndex2017)[1])
MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientCode = PatientID)

plot(FluidsIndex2017[[1]]$TimeSeriesData$time - DP_StripTime(MetaData$FirstITUEntry) , FluidsIndex2017[[1]]$TimeSeriesData$Hourly.In - FluidsIndex2017[[1]]$TimeSeriesData$Hourly.Out , type ='l'  , xlim = c(0 , 3600*5*24) , ylim = c(-1000,1000), col =rgb(0,0,1, alpha = 0.05))

for(i in 2:length(FluidsIndex2017)){
  PatientID <- DP_ExtractPatientIDFromNewPatinetID(names(FluidsIndex2017[i]))
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientCode = PatientID)
  
  if( sum(MetaData$ProcDetails[1] %in% setofops) == 0  ){next}
  if( is.na(MetaData$ConfirmedFirstNewAF[1])){MetaData$ConfirmedFirstNewAF[1] = 'NAF'}
  if( !is.na(MetaData$FirstNewAF[1]) & (MetaData$ConfirmedFirstNewAF[1] != 'CNAF')){
  lines(FluidsIndex2017[[i]]$TimeSeriesData$time - DP_StripTime(MetaData$FirstITUEntry[1]) , FluidsIndex2017[[i]]$TimeSeriesData$Hourly.In - FluidsIndex2017[[i]]$TimeSeriesData$Hourly.Out , type ='l' , col =rgb(1,0,0, alpha = 0.05))
  }else{
  #  lines(FluidsIndex2017[[i]]$TimeSeriesData$time - DP_StripTime(MetaData$FirstITUEntry[1]) , FluidsIndex2017[[i]]$TimeSeriesData$Hourly.In - FluidsIndex2017[[i]]$TimeSeriesData$Hourly.Out , type ='l' , col =rgb(0,0,1, alpha = 0.05))
    
    }
}
