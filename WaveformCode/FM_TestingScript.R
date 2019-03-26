
PatientID <- DP_choosepatient( listAllPatients )

MetaData <- DP_ExtractPatientRecordforIndex( PatIndex2017 = PatIndex2017 , PatientCode = PatientID )
ECGs <- DP_LoadReducedECGs( path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile( path , PatientID )

{
  StartBeat <- 1
  numberofBeats <- 500
  
  RegularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  RegularyIrregularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  minImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  minImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  meanImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  meanImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  MulImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  MulImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
  timemat <- rep(mean(PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,1]) , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) )
  
  for(i in 1:ceil(dim(RPeakData$RRCombined)[1]/numberofBeats)){
    #for(i in 1:11 ){
    
    rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
    if(max(rangeofbeats) > dim(RPeakData$RRCombined)[1]){break}
    
    RRtimes <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,3]
    timemat[i] <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats[length(rangeofbeats)],1]
    
    ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , F_x = F_xreg ,f_x = f_xreg ,  MD = MD_Reg , RRtimes = RRtimes , Corr_sdhat = Corr_sdhat2 , imthreshold = ImThreshold2[2])
    ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,x = x ,F_x =  F_x_ReIre, f_x = f_x_ReIre,  MD = MD_ReIre, RRtimes = RRtimes , Corr_sdhat = Corr_sdhat1, imthreshold  = ImThreshold1[2] )
    
    #ReHmOutput <- FM_HistoryMatchRRDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , f_x = f_xreg , MD = MD_Reg , RRtimes = RRtimes , imthreshold = 0.5)
    #ReIreHmOutput <- FM_HistoryMatchRRDensity(PriorNonImplausibleSetRegularyIreRegular,x , f_x_ReIre, MD = MD_ReIre, RRtimes, imthreshold = 0.5)
    
    minImReg[i] <- min(ReHmOutput$Implausability[,1])
    minImRegIre[i] <- min(ReIreHmOutput$Implausability[,1])
    
    meanImReg[i] <- min(ReHmOutput$Implausability[,2])
    meanImRegIre[i] <- min(ReIreHmOutput$Implausability[,2])
    
    MulImReg[i] <- min(ReHmOutput$Implausability[,3])
    MulImRegIre[i] <- min(ReIreHmOutput$Implausability[,3])
    
    if( ( ( meanImReg[i] < ImThreshold2[1] )*( minImReg[i] <  ImThreshold2[2] )*( MulImReg[i] <  ImThreshold2[3] )) ==1 ){
      print('Regular')
      RegularLogical[i] <- 1
    }
    if( ( ( meanImRegIre[i] < ImThreshold1[1] )*( minImRegIre[i] <  ImThreshold1[2] )*( MulImRegIre[i] <  ImThreshold1[3] ) ) ==1 ){
      print('Regularly Irregular')
      RegularyIrregularLogical[i] <- 1
      print(StartBeat)
    }
    if( RegularLogical[i] ==0 && RegularyIrregularLogical[i] == 0 ){
      print('Irregular')
      RegularyIrregularLogical[i] <- 1
      print(StartBeat)
    }
   
    
    StartBeat <- StartBeat + numberofBeats
    DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
  }
}

{
  RegularyIrregularLogical <- ( (meanImRegIre < ImThreshold1[1] )*(minImRegIre <  ImThreshold1[2])*( MulImRegIre <  ImThreshold1[3] ) ) ==1
  RegularLogical <- ( ( (meanImReg < ImThreshold2[1] )*( minImReg <  ImThreshold2[2])*( MulImReg <  ImThreshold2[3] ) ) ==1)
  
  
  RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularLogical  , minutethreshold = 1)
  RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularyIrregularLogical  , minutethreshold = 1)
  IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = ((RegularyIrregularLogical ==0)*(RegularLogical==0)) == 1 , minutethreshold = 1)
  
  RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( PatientID , ' RRTimes' ) )
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
  x11(40 , 25)
  RRPlot + geom_line(data = data.frame( t = timemat , P = ((1/MulImReg) / ((1/MulImReg) + (1/MulImRegIre) ))) , aes(t, P) , col ='blue')+
    geom_hline(yintercept = 0.6)
}



for(PatientID in listAllPatients[180:length(listAllPatients)] ){

MetaData <- DP_ExtractPatientRecordforIndex( PatIndex2017 = PatIndex2017 , PatientCode = PatientID )
ECGs <- DP_LoadReducedECGs( path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile( path , PatientID )

{
  StartBeat <- 1
  numberofBeats <- 500
  
  RegularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  RegularyIrregularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  minImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  minImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  meanImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  meanImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  MulImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  MulImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  
  timemat <- rep(mean(PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,1]) , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) )
  
  for(i in 1:ceil(dim(RPeakData$RRCombined)[1]/numberofBeats)){
    #for(i in 1:11 ){
    
    rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
    if(max(rangeofbeats) > dim(RPeakData$RRCombined)[1]){break}
    
    RRtimes <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,3]
    timemat[i] <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats[length(rangeofbeats)],1]
    
    if(length(RRtimes) < 100 ){
      next
    }
    
    ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , F_x = F_xreg ,f_x = f_xreg ,  MD = MD_Reg , RRtimes = RRtimes , Corr_sdhat = Corr_sdhat2 , imthreshold = ImThreshold2[2])
    ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,x = x ,F_x =  F_x_ReIre, f_x = f_x_ReIre,  MD = MD_ReIre, RRtimes = RRtimes , Corr_sdhat = Corr_sdhat1, imthreshold  = ImThreshold1[2] )
    
    #ReHmOutput <- FM_HistoryMatchRRDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , f_x = f_xreg , MD = MD_Reg , RRtimes = RRtimes , imthreshold = 0.5)
    #ReIreHmOutput <- FM_HistoryMatchRRDensity(PriorNonImplausibleSetRegularyIreRegular,x , f_x_ReIre, MD = MD_ReIre, RRtimes, imthreshold = 0.5)
    
    minImReg[i] <- min(ReHmOutput$Implausability[,1])
    minImRegIre[i] <- min(ReIreHmOutput$Implausability[,1])
    
    meanImReg[i] <- min(ReHmOutput$Implausability[,2])
    meanImRegIre[i] <- min(ReIreHmOutput$Implausability[,2])
    
    MulImReg[i] <- min(ReHmOutput$Implausability[,3])
    MulImRegIre[i] <- min(ReIreHmOutput$Implausability[,3])
    
    if( ( ( meanImReg[i] < ImThreshold2[1] )*( minImReg[i] <  ImThreshold2[2] )*( MulImReg[i] <  ImThreshold2[3] )) ==1 ){
      print('Regular')
      RegularLogical[i] <- 1
    }
    if( ( ( meanImRegIre[i] < ImThreshold1[1] )*( minImRegIre[i] <  ImThreshold1[2] )*( MulImRegIre[i] <  ImThreshold1[3] ) ) ==1 ){
      print('Regularly Irregular')
      RegularyIrregularLogical[i] <- 1
      print(StartBeat)
    }
    if( RegularLogical[i] ==0 && RegularyIrregularLogical[i] == 0 ){
      print('Irregular')
      RegularyIrregularLogical[i] <- 1
      print(StartBeat)
    }
    
    
    StartBeat <- StartBeat + numberofBeats
    DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
  }
}

{
  RegularyIrregularLogical <- ( (meanImRegIre < ImThreshold1[1] )*(minImRegIre <  ImThreshold1[2])*( MulImRegIre <  ImThreshold1[3] ) ) ==1
  RegularLogical <- ( ( (meanImReg < ImThreshold2[1] )*( minImReg <  ImThreshold2[2])*( MulImReg <  ImThreshold2[3] ) ) ==1)
  
  
  RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularLogical  , minutethreshold = 1)
  RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularyIrregularLogical  , minutethreshold = 1)
  IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = ((RegularyIrregularLogical ==0)*(RegularLogical==0)) == 1 , minutethreshold = 1)
  
  RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( PatientID , ' RRTimes' ) )
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
  
  pdf(file = paste0("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\Annotation" , PatientID , '.pdf') )
  print(RRPlot + geom_line(data = data.frame( t = timemat , P = ((1/MulImReg) / ((1/MulImReg) + (1/MulImRegIre) ))) , aes(t, P) , col ='blue')+
    geom_hline(yintercept = 0.6))
  dev.off()
  #x11(40 , 25)
  #print(RRPlot + geom_line(data = data.frame( t = timemat , P = ((1/MulImReg) / ((1/MulImReg) + (1/MulImRegIre) ))) , aes(t, P) , col ='blue')+
  #  geom_hline(yintercept = 0.6))
  
}
}
