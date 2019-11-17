
for(PatientID in listAllPatients[26:length(listAllPatients)] ){
  {
    MetaData <- DP_ExtractPatientRecordforIndex( PatIndex2017 = PatIndex2017 , PatientCode = PatientID )
    if(DP_CheckIfAFPatient(MetaData)){next}
    if(!DP_CheckECGreducedfilesprocessed(path , PatientID , Filestoprocess = 'ECGII_reduced')){next}
    ECGs <- DP_LoadReducedECGs( path ,  PatientID , FilestoProcess = FilestoProcess  )
    RPeakData <- DP_LoadRpeaksfile( path , PatientID )
    
    dovalidationplots = 0
    
    StartBeat <- 1
    numberofBeats <- 500
    
    RegularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    RegularyIrregularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    RegularLogical2 <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    RegularyIrregularLogical2 <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    SecondWaveLogical <-  matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    ECGAbscence <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    RpeaksFail <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    
    StartBeatMat <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
    timemat <- rep( RPeakData$RRCombined[1,1], ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) )
    
    outputstruct <- list()
    
    for(i in 1:ceil(dim(RPeakData$RRCombined)[1]/numberofBeats)){
      #for(i in 1:11 ){
      
      rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
      if(max(rangeofbeats) > dim(RPeakData$RRCombined)[1]){break}
      if(length(rangeofbeats) < 100){
        StartBeat <- StartBeat + numberofBeats
        next}
      
      RRtimes <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,3]
      
      t_tmp <- length(rangeofbeats)/as.numeric(abs(diff(range(RPeakData$RRCombined[rangeofbeats,1])))) / (60/median(RRtimes))
      timemat[i] <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats[length(rangeofbeats)],1]
      StartBeatMat[i] <- StartBeat
      
      if(abs(1-t_tmp)>0.1){
        disp(t_tmp)
        RpeaksFail[i] <- 1
        StartBeat <- StartBeat + numberofBeats
        next
      }
      if(var( RRtimes)==0){next}
      
      # HM Heart-rhythm
      ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, 
                                                       x = xreg ,
                                                       F_x = F_xreg ,
                                                       f_x = f_xreg ,  
                                                       MD = MD_Reg , 
                                                       RRtimes = RRtimes , 
                                                       Corr_sdhat = Corr_sdhat2 , 
                                                       #imthreshold = ImThreshold2[2],
                                                       #imthreshold2 = ImThreshold2[1])
                                                       imthreshold = ImThresholdMaxRegular,
                                                       imthreshold2 = ImThresholdMeanRegular )
      ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,
                                                          x = xIrIreg ,
                                                          F_x =  F_x_ReIre,
                                                          f_x = f_x_ReIre,
                                                          MD = MD_ReIre,
                                                          RRtimes = RRtimes ,
                                                          Corr_sdhat = Corr_sdhat1, 
                                                          #imthreshold = ImThreshold1[2],
                                                          #imthreshold2 = ImThreshold1[1]     )
                                                          imthreshold  = ImThresholdMaxIrregularlyIrregular,
                                                          imthreshold2 = ImThresholdMeanIrregularlyIrregular )
      TotalHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetTotal,
                                                          x = xTotal ,
                                                          F_x =  F_total, 
                                                          f_x = f_total,  
                                                          MD = MD_Total , 
                                                          RRtimes = RRtimes , 
                                                          Corr_sdhat = Corr_sdhat1, 
                                                          imthreshold  = ImThresholdMaxTotal , 
                                                          imthreshold2 = ImThresholdMeanTotal )
      
      if(sum((ReIreHmOutput$Implausability[,1] < max(ImThresholdMaxIrregularlyIrregular))&(ReIreHmOutput$Implausability[,2] < max(ImThresholdMeanIrregularlyIrregular))) > 1){
        RegularyIrregularLogical2[i] <- T
        disp('Irregularly-irregular Heart-rhythm')
      }else{
        RegularyIrregularLogical2[i] <- F
      }
      
      if(sum((ReHmOutput$Implausability[,1] < max(ImThresholdMaxRegular))&(ReHmOutput$Implausability[,2] < max(ImThresholdMeanIrregularlyIrregular))) > 1){
        RegularLogical2[i] <- T
        disp('Regular Heart-rhythm')
      }else{
        RegularLogical2[i] <- F
      }
      
      
      # HM P-waves
      source('FM_HistoryMatchMeanPWave.R')
      
      if(length(NonImplausibleX) > 5){
        if(sum(NonImplausibleX[,5] == 0) > 0 ){
          RegularyIrregularLogical[i] = 1
          disp('Pwave Abscent')
        }
        if(sum(NonImplausibleX[,5] != 0) > 0 ){
          RegularLogical[i] = 1
          disp('Pwave Present')
        }
        
        # second wave analysis
        if(Wave2Pwave ){
          RegularLogical[i] <- 0
          RegularyIrregularLogical[i] <- 1
        }else{
          RegularLogical[i] <- 1
          RegularyIrregularLogical[i] <- 0
        }
        
        
        # Pwave both nonimplausible logic.
        if( (MeanIm[length(MeanIm)] - min(MeanIm[1:(length(MeanIm) -2)])) < -0.5  ){
          RegularLogical[i] <- 0
          RegularyIrregularLogical[i] <- 1
        }
        if( (MeanIm[length(MeanIm)] - min(MeanIm[1:(length(MeanIm) -2)])) > 0.5  ){
          RegularLogical[i] <- 1
          RegularyIrregularLogical[i] <- 0
        }
      }
      
      if((RegularLogical[i] == 0) & (RegularLogical2[i] ==1) & (RegularyIrregularLogical2[i] == 1) & (RegularyIrregularLogical[i] == 1)){
        SecondWaveLogical[i] <- FM_DubiousCaseLogic(RRtimes = RRtimes[!is.na(RRtimes)],
                                                    ReIreHmOutput = ReIreHmOutput,
                                                    ReHmOutput = ReHmOutput,
                                                    ObservedIm_xx = ObservedIm_xx,
                                                    EmulatorParametersCDFMean = EmulatorParametersCDFMean
                                                    , EmulatorParametersCDFMax = EmulatorParametersCDFMax ,
                                                    RegularLogical)
        if(SecondWaveLogical[i] == 1){
          disp('AFib Clasified in Dubious Case')
        }else{
          disp('AFib Not Clasified in Dubious Case')
        }
      }
      StartBeat <- StartBeat + numberofBeats
      DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
      outputstruct[[i]] <- (list(NonImplausibleX , ReHmOutput , ReIreHmOutput ,TotalHmOutput , timemat[i]) )
      
    }
    
    
  }
  
  
  
  {
    RegularLogicalTotal <- (RegularLogical ==1) & (RegularLogical2 ==1) & (RegularyIrregularLogical == 0)&(RegularyIrregularLogical2 == 0)
    IrregularlyIrregularLogical <- (RegularyIrregularLogical == 1)&(RegularyIrregularLogical2 == 1)&(RegularLogical ==0) & (RegularLogical2 ==0)
    IrregularlyIrregularLogical[SecondWaveLogical == 1] <- 1
    RegularLogicalTotal[SecondWaveLogical == 1] <- 0 
    BadDataLogical <- RpeaksFail ==1 | ECGAbscence ==1
    
    RegularlyIrregularLogical <- (RegularLogicalTotal ==0) &(IrregularlyIrregularLogical ==0)
    
    RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = IrregularlyIrregularLogical , minutethreshold = 1)
    RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularLogicalTotal , minutethreshold = 1)
    IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularlyIrregularLogical, minutethreshold = 1)
    BadDataTimePoints <-  ASWF_GetStartEndAF(timemat , logicaltimeseries = BadDataLogical, minutethreshold = 1)
    #UndecidedtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = Undecided, minutethreshold = 1)
    
    
    RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( PatientID , ' RRTimes' ) )
    RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
    RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
    RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
    RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = BadDataTimePoints , fillcolor = 'black')
    pdf(file = paste0("D:\\AllPatientAnnotations\\Annotation" , PatientID , '.pdf') )
    print(RRPlot)
    dev.off()
    
    outputstruct[[i+1]] <- cbind(RegularLogical, RegularLogical2, RegularyIrregularLogical,RegularyIrregularLogical2,SecondWaveLogical)
    
    save(outputstruct , file = paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData'))
  }
}
