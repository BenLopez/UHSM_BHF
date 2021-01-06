{
  
  
  MetaData <- DP_ExtractPatientRecordforIndex( PatIndex2017 = PatIndex2017 , PatientCode = PatientID )
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
  ECGAbscence <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  RpeaksFail <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  
  StartBeatMat <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
  timemat <- rep(mean(PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,1]) , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) )
  
  outputstruct <- list()
  
  for(i in 1:ceil(dim(RPeakData$RRCombined)[1]/numberofBeats)){
    #for(i in 1:11 ){
    rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
    if(max(rangeofbeats) > dim(RPeakData$RRCombined)[1]){break}
    if(length(rangeofbeats) < 100){
      StartBeat <- StartBeat + numberofBeats
      next}
    
    RRtimes <- RPeakData$RRCombined[rangeofbeats,3]
    t_tmp <- length(rangeofbeats)/as.numeric(abs(diff(range(RPeakData$RRCombined[rangeofbeats,1])))) / (60/median(RRtimes))
    timemat[i] <-  RPeakData$RRCombined[rangeofbeats[length(rangeofbeats)],1]
    StartBeatMat[i] <- StartBeat
    
    if(median(RRtimes) > 0.6){
      if(abs(1-t_tmp)>0.1){
        disp(t_tmp)
        RpeaksFail[i] <- 1
        StartBeat <- StartBeat + numberofBeats
        next
      }
    }
    # HM Heart-rhythm
    TotalHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetTotal,
                                                        x = xTotal ,
                                                        F_x =  F_total, 
                                                        f_x = f_total,  
                                                        MD = MD_Total , 
                                                        RRtimes = RRtimes , 
                                                        Corr_sdhat = Corr_sdhat1, 
                                                        imthreshold  = ImThresholdMaxTotal , 
                                                        imthreshold2 = ImThresholdMeanTotal )
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
                                                     imthreshold2 = ImThresholdMeanRegular)
    ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,
                                                        x = xIrIreg ,
                                                        F_x =  F_x_ReIre,
                                                        f_x = f_x_ReIre,
                                                        MD = MD_ReIre,
                                                        RRtimes = RRtimes ,
                                                        Corr_sdhat = Corr_sdhat1, 
                                                        # imthreshold = ImThreshold1[2],
                                                        # imthreshold2 = ImThreshold1[1])
                                                        imthreshold = ImThresholdMaxIrregularlyIrregular,
                                                        imthreshold2 = ImThresholdMaxIrregularlyIrregular)
    
    if(median(RRtimes) > 0.6 & median(RRtimes) < 1 & (quantile(RRtimes , 0.95)  - quantile(RRtimes , 0.05)) < 0.1){
      RegularLogical2[i] <- T
      RegularyIrregularLogical2[i] <- F
      disp('Regular Heart-rhythm by IQR alone.')
    }else{
      if(!is.null(ReIreHmOutput$NonImplausibleSets) ){
        RegularyIrregularLogical2[i] <- T
        disp('Irregularly-irregular Heart-rhythm')
      }else{
        RegularyIrregularLogical2[i] <- F
      }
      if(!is.null(ReHmOutput$NonImplausibleSets) ){
        RegularLogical2[i] <- T
        disp('Regular Heart-rhythm')
      }else{
        RegularLogical2[i] <- F
      }
    }
    
    if(RegularLogical2[i] ==1 & RegularyIrregularLogical2[i] ==1 ){
      SecondWaveResult <- FM_DubiousCaseLogic(RRtimes = RRtimes[!is.na(RRtimes)],
                                              ReIreHmOutput = ReIreHmOutput,
                                              ReHmOutput = ReHmOutput,
                                              ObservedIm_xx = ObservedIm_xx,
                                              EmulatorParametersCDFMean = EmulatorParametersCDFMean
                                              , EmulatorParametersCDFMax = EmulatorParametersCDFMax ,
                                              RegularLogical)
      RegularyIrregularLogical2[i] <- SecondWaveResult
      RegularLogical2[i] <- !SecondWaveResult 
      if(SecondWaveResult){
        disp('AFib Clasified in Dubious Case')
      }else{
        disp('AFib Not Clasified in Dubious Case')
      }
    }
    
    ##### History Match P-waves #####
    source('FM_HistoryMatchMeanPWave.R')
    if(is.na(RegularLogical2[i])){
      break
    }
    # legacy code  
    # Pwave both nonimplausible logic.
    #if( (MeanIm[length(MeanIm)] - min(MeanIm[1:(length(MeanIm) -2)])) < -0.5  ){
    #  RegularLogical[i] <- 0
    #  RegularyIrregularLogical[i] <- 1
    #}
    #if( (MeanIm[length(MeanIm)] - min(MeanIm[1:(length(MeanIm) -2)])) > 0.5  ){
    #  RegularLogical[i] <- 1
    #  RegularyIrregularLogical[i] <- 0
    #}
    
    
    StartBeat <- StartBeat + numberofBeats
    DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
    if(exists('outputstruct')){
    outputstruct[[i]] <- (list(NonImplausibleX , ReHmOutput , ReIreHmOutput , TotalHmOutput , timemat[i]) )
    }
    }
}