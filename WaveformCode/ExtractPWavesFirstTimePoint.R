EmulatorParameters <- PWaveHM_CreateDefaultEmulationclass()
FilestoProcess = DP_ChooseECGstoProcess()

{Xstar = seq(0.5 ,1 , 0.01)
  PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.001 , 0.95, 0.001 ) , b = c(0.6  , 0.04 , 0.6 , 0.04 ) , numbersamples = 100000 )
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.1 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.02 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]
  QSwidth = 10 }


for(PatientID in listAllPatients[1:length(listAllPatients)]){
  
{  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = PatientID)
  #if(!DP_CheckIfAFPatient(MetaData)){next}
  if(!DP_CheckECGreducedfilesprocessed(path = path , PatientsId = PatientID , Filestoprocess = 'ECGI_reduced')){next}
  ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
  RPeakData <- DP_LoadRpeaksfile(path , PatientID)
  if(DP_CheckIfAFPatient(MetaData)){
    RPeakData$RRCombined <- RPeakData$RRCombined[RPeakData$RRCombined$t < DP_StripTime(MetaData$ConfirmedFirstNewAF[1]) , ]
  }
  if(dim(RPeakData$RRCombined)[1] < 2000){next}
  
  QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[1000:1500,] , QSwidth = QSwidth)
  
  if(dim(as.matrix(QS_Struct$Date))[1] < 250){next}
  if(length(QS_Struct$numvalues) == 0){next}
  if(length(QS_Struct$Date) == 0){next}
  
  NumberofBeats <- 500
  StartBeat = 1000
  
  PWaveSummariesOutputs <- matrix(0 , 12 ,  1)
  
  for(kk in 1){
    if(length(RPeakData$RRCombined$RR[ max(0,StartBeat):min((StartBeat + NumberofBeats) , length(RPeakData$RRCombined$RR) - 450) ]) < NumberofBeats){next}
    
    QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[ max(0,StartBeat):min((StartBeat + NumberofBeats) , length(RPeakData$RRCombined$RR) - 1000), ], QSwidth = QSwidth)
    if(length(QS_Struct$t_start) < 200){next}
    if(sum(!is.na(QS_Struct$Value[1,])) < 50 ){next}
    EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct = QS_Struct , EmulatorParameters = EmulatorParameters ,Xstar =  Xstar )
    mQS <- apply(DP_RemoveNaRows(EmulatedQS) , 2 , median) 
    mQS <- mQS - mean(mQS)
    vQS <- apply(DP_RemoveNaRows(EmulatedQS) , 2 , IQR)
    
    Implausability <- PWaveHM_CalulateImplausabilityTQSegment( Xstar = Xstar , z = mQS ,  PriorNonImplausibleSet =  PriorNonImplausibleSet , ImplausabilityFunction = CalculateImplausability )
    Xmin <- PriorNonImplausibleSet[which.min(Implausability) , ] 
    
    PWaveSummariesOutputs[ 1 , kk] <- mQS[which.max( PsimulatorFunction(Xmin , Xstar))]  # P-amplitude
    PWaveSummariesOutputs[ 2 , kk] <- Xstar[which.max( PsimulatorFunction(Xmin , Xstar))]  # P-amplitude location
    PWaveSummariesOutputs[ 3 , kk] <- Xstar[which(abs(c(0, diff(PsimulatorFunction(Xmin , Xstar)> 0.1))  )==1 )[1]]# P-start
    if(length(which(abs(c(0, diff(PsimulatorFunction(Xmin , Xstar)> 0.1))  )==1) ) == 1){
      PWaveSummariesOutputs[ 4 , kk] <- 1}else{
        PWaveSummariesOutputs[ 4 , kk] <- Xstar[which(abs(c(0, diff(PsimulatorFunction(Xmin , Xstar)> 0.1))  )==1 )[2]]}   # P - end
    PWaveSummariesOutputs[ 5 , kk] <- PWaveSummariesOutputs[4,kk] - PWaveSummariesOutputs[3,kk] # width
    PWaveSummariesOutputs[ 6 , kk] <- Xmin[3] - Xmin[1] # Difference between beats  
    PWaveSummariesOutputs[ 7 , kk] <- mean(mQS[Xstar > PWaveSummariesOutputs[3,kk]]) # PR level
    tmpdiff <- c(0, diff(PsimulatorFunction(Xmin , Xstar)))
    PWaveSummariesOutputs[ 8 , kk] <- max(tmpdiff[((Xstar > PWaveSummariesOutputs[3,kk])*(Xstar < PWaveSummariesOutputs[4,kk]) )==1])
    PWaveSummariesOutputs[ 9 , kk] <- min(tmpdiff[((Xstar > PWaveSummariesOutputs[3,kk])*(Xstar < PWaveSummariesOutputs[4,kk]) )==1])
    PWaveSummariesOutputs[ 10 , kk] <- mean(vQS[((Xstar > PWaveSummariesOutputs[3,kk])*(Xstar < PWaveSummariesOutputs[4,kk]) )==1])
    PWaveSummariesOutputs[ 11 , kk] <- max(vQS[((Xstar > PWaveSummariesOutputs[3,kk])*(Xstar < PWaveSummariesOutputs[4,kk]) )==1])
    PWaveSummariesOutputs[ 12 , kk] <- min(Implausability)
    
    # Extract P-wave Data
    #DP_WaitBar(kk/dim(RRDistributionSummariesOutput)[2])
    StartBeat <- StartBeat + NumberofBeats + 1
  }
}  
  # Save
  output <- PWaveSummariesOutputs
  DP_SaveFile(output , path = path , PatientID = PatientID , paste0(PatientID , '_PwaveSummariesfirst'))
  
  DP_WaitBar(which(listAllPatients == PatientID)/length(listAllPatients))
  
}
