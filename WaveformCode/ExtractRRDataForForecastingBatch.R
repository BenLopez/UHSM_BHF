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
  set.seed(1)
}

FilestoProcess = DP_ChooseECGstoProcess()



for(PatientID in listAllPatients[1:length(listAllPatients)]){
  
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = PatientID)
  if(!DP_CheckECGreducedfilesprocessed(path = path , PatientsId = PatientID , Filestoprocess = 'ECGI_reduced')){next}
  RPeakData <- DP_LoadRpeaksfile(path , PatientID)
  
  if(DP_CheckIfAFPatient(MetaData)){
    RPeakData$RRCombined <- RPeakData$RRCombined[RPeakData$RRCombined$t < DP_StripTime(MetaData$ConfirmedFirstNewAF[1]) , ]
  }
  if(dim(RPeakData$RRCombined)[1] < 2000){next}
  
  
  NumberofBeats <- 500
  StartBeat = 1000
RRDistributionSummariesOutput <- matrix(0 , 13 , floor((length(RPeakData$RRCombined$RR) - 1000)/NumberofBeats) )
#for(kk in 1:dim(RRDistributionSummariesOutput)[2]){
for(kk in 1:50){
    if(length(RPeakData$RRCombined$RR[ max(0,StartBeat):min((StartBeat + NumberofBeats) , length(RPeakData$RRCombined$RR) - 1000) ]) < NumberofBeats){next}

    RRTimes <- RPeakData$RRCombined$RR[ max(0,StartBeat):min((StartBeat + NumberofBeats) , length(RPeakData$RRCombined$RR) - 1000) ]
    tmpKde <- kde(RRTimes + sqrt(0.00001)*rnorm(length(RRTimes)))
    
    RRDistributionSummariesOutput[1 , kk] <- mean(RRTimes)
    RRDistributionSummariesOutput[2 , kk] <- median(RRTimes)
    RRDistributionSummariesOutput[3 , kk] <- var(RRTimes)
    RRDistributionSummariesOutput[4 , kk] <- IQR(RRTimes)
    RRDistributionSummariesOutput[5 , kk] <- skewness(RRTimes)
    RRDistributionSummariesOutput[6 , kk] <- kurtosis(RRTimes)
    RRDistributionSummariesOutput[7 , kk] <- max(tmpKde$estimate)
    RRDistributionSummariesOutput[8 , kk] <- var(tmpKde$estimate/max(tmpKde$estimate))
    PeaksLogical <- PE_FindLocalTurningPoints(tmpKde$estimate/max(tmpKde$estimate) > 0.1 , tmpKde$estimate) 
    RRDistributionSummariesOutput[9 , kk] <- sum(PeaksLogical) # number of modes
    if(sum(PeaksLogical) > 1 ){
      RRDistributionSummariesOutput[10 , kk] <- var((tmpKde$estimate/max(tmpKde$estimate))[PeaksLogical])
      RRDistributionSummariesOutput[11 , kk] <- var(tmpKde$eval.points[PeaksLogical])
    }else{
      RRDistributionSummariesOutput[10 , kk] <- 0
      RRDistributionSummariesOutput[11 , kk] <- 0
    }
    RRDistributionSummariesOutput[12 , kk] <- approx_entropy(RRTimes)
    RRDistributionSummariesOutput[13 , kk] <- mean(RPeakData$RRCombined$t[ max(0,StartBeat):min((StartBeat + NumberofBeats) , length(RPeakData$RRCombined$RR) - 1000) ])
    # Extract P-wave Data
    StartBeat <- StartBeat + NumberofBeats + 1
    
  }
  # Save
  output <- RRDistributionSummariesOutput
  DP_SaveFile(output , path = path , PatientID = PatientID , paste0(PatientID , '_RRDistributionSummariesForeCasting'))
DP_WaitBar(which(listAllPatients == PatientID)/length(listAllPatients))
  }
