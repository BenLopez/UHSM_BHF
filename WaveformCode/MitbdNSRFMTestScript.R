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

DP_ChooseDataReps()
source('FM_CreateMeanPWavePriors.R')
source('FM_CreateRhythumPriors.R')
source('CTEm_LoadDataandCreateEmulatorStructures.R')


{
  PatientID <- DP_choosepatient( listAllPatients )
  ECGs <- DP_LoadReducedECGs( path ,  PatientID , FilestoProcess = FilestoProcess  )
  ECGs$ECGII$Value <- 100*ECGs$ECGII$Value
  RPeakData <- DP_LoadRpeaksfile( path , PatientID )
  
  dovalidationplots = 0
  
  StartBeat <- 1
  numberofBeats <- 500
  
  RegularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  RegularyIrregularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  RegularLogical2 <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  RegularyIrregularLogical2 <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
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
    timemat[i] <-  RPeakData$RRCombined[rangeofbeats[length(rangeofbeats)],1]
    StartBeatMat[i] <- StartBeat
    
    # HM Heart-rhythm
    ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, 
                                                     x = x ,
                                                     F_x = F_xreg ,
                                                     f_x = f_xreg ,  
                                                     MD = MD_Reg , 
                                                     RRtimes = RRtimes , 
                                                     Corr_sdhat = Corr_sdhat2 , 
                                                     imthreshold = ImThresholdMaxRegular,
                                                     imthreshold2 = ImThresholdMeanRegular)
    ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,
                                                        x = x ,
                                                        F_x =  F_x_ReIre,
                                                        f_x = f_x_ReIre,
                                                        MD = MD_ReIre,
                                                        RRtimes = RRtimes ,
                                                        Corr_sdhat = Corr_sdhat1,
                                                        imthreshold  = ImThresholdMaxIrregularlyIrregular,
                                                        imthreshold2 = ImThresholdMeanIrregularlyIrregular)
    
    if(sum((ReIreHmOutput$Implausability[,1] < ImThreshold1[2])&(ReIreHmOutput$Implausability[,2] < ImThreshold1[1])&(ReIreHmOutput$Implausability[,3] < ImThreshold1[3])) > 1){
      RegularyIrregularLogical2[i] <- T
      disp('Irregularly-irregular Heart-rhythm')
    }else{
      RegularyIrregularLogical2[i] <- F
    }
    
    if(sum((ReHmOutput$Implausability[,1] < ImThreshold2[2])&(ReHmOutput$Implausability[,2] < ImThreshold2[1])&(ReHmOutput$Implausability[,3] < ImThreshold2[3])) > 1){
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
      
      if(sum(NonImplausibleX[,5] != 0) > 0 & sum(NonImplausibleX[,5] == 0) > 0  ){
        if(ImDiff < -0.5){
          RegularyIrregularLogical[i] = 1
          RegularLogical[i] = 0}
      }
      if(sum(NonImplausibleX[,5] != 0) > 0 & sum(NonImplausibleX[,5] == 0) > 0  ){
        if(ImDiff > 0.5){
          
          RegularyIrregularLogical[i] = 0
          RegularLogical[i] = 1
        }
      } 
    }
    if(length(NonImplausibleX) == 5){
      if(sum(NonImplausibleX[5] == 0) > 0 ){
        RegularyIrregularLogical[i] = 1
        disp('Pwave Abscent')
      }
      if(sum(NonImplausibleX[5] != 0) > 0 ){
        RegularLogical[i] = 1
        disp('Pwave Present')
      }
    }
    
    StartBeat <- StartBeat + numberofBeats
    DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
    #outputstruct[[i]] <- (list(NonImplausibleX , ReHmOutput , ReIreHmOutput , cbind(MaxIm[tmplog] , MeanIm[tmplog])) )
    
  }
  
}

{
  #RegularyIrregularLogical2 <- ( (meanImRegIre < ImThreshold1[1] )*(minImRegIre <  ImThreshold1[2])*( MulImRegIre <  ImThreshold1[3] ) ) ==1
  #RegularLogical2 <- ( ( (meanImReg < ImThreshold2[1] )*( minImReg <  ImThreshold2[2])*( MulImReg <  ImThreshold2[3] ) ) ==1)
  MetaData = DP_CreateDummyMetaData(PatIndex2017 = PatIndex2017)  
  
  RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = ((RegularyIrregularLogical*RegularyIrregularLogical2*(RegularLogical ==0)*(RegularLogical2 ==0))==1)  , minutethreshold = 1)
  RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = (RegularLogical*RegularLogical2)==1  , minutethreshold = 1)
  IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = (((RegularyIrregularLogical*RegularyIrregularLogical2) ==0)*((RegularLogical*RegularLogical2)==0)) == 1 , minutethreshold = 1)

  RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( PatientID , ' RRTimes' ) )
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
  x11(40 , 25)
  RRPlot + ggtitle('PWaves')
  
  #pdf(file = paste0("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\Annotation" , PatientID , '.pdf') )
  #print(RRPlot)
  #dev.off()
  
  #save(outputstruct , file = paste0("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\HMOutput" , PatientID , '.RData'))
}




{
  StartBeat <- 1001
  rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
  RRtimes <-  RPeakData$RRCombined[rangeofbeats,3]
  mm <- FM_EmulatorEstimate( Y = RRtimes )
  m <- median(RRtimes)
  #m <- 0
  #mm <- 0
  tmp <- RRtimes + m - mm + rnorm(length(RRtimes) , 0 , 0.0025 )
  tmp <- tmp[!is.na(tmp)]
  kdmodel <- kde( tmp )
  ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, 
                                                   x = x , 
                                                   F_x = F_xreg ,
                                                   f_x = f_xreg ,  
                                                   MD = MD_Reg , 
                                                   RRtimes = RRtimes , 
                                                   Corr_sdhat = Corr_sdhat2 , 
                                                   imthreshold  = ImThresholdMaxRegular,
                                                   imthreshold2 = ImThresholdMeanRegular)
  ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,
                                                      x = x ,
                                                      F_x =  F_x_ReIre,
                                                      f_x = f_x_ReIre,
                                                      MD = MD_ReIre,
                                                      RRtimes = RRtimes ,
                                                      Corr_sdhat = Corr_sdhat1,
                                                      imthreshold  = ImThresholdMaxIrregularlyIrregular,
                                                      imthreshold2 = ImThresholdMeanIrregularlyIrregular )
  
  x11()
  plot(x , predict( kdmodel,x = x) , type ='l')
  
  if(length(ReHmOutput$Implausability) > 4){
    tmp <- ReHmOutput$f_x[(ReHmOutput$Implausability[,1] < ImThreshold2[2])&(ReHmOutput$Implausability[,2] < ImThreshold2[1])&(ReHmOutput$Implausability[,3] < ImThreshold2[3]),]
    #tmp <- ReHmOutput$f_x[(ReHmOutput$Implausability[,1] < ImThreshold2[2]),]
    #plot(x , predict( kdmodel,x = x) , type ='l')
    for(i in 1:dim(tmp)[1]){
      #lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
      lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
      
    }
    BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[1:500,] , PriorNonImplausibleSetRegular[1:500,] , ReHmOutput$NonImplausibleSets[(ReHmOutput$Implausability[,1] < ImThreshold2[2])&(ReHmOutput$Implausability[,2] < ImThreshold2[1])&(ReHmOutput$Implausability[,3] < ImThreshold2[3]) , ] )  
    
  }
  if((ReIreHmOutput$Implausability[,1] < ImThreshold1[2])&(ReIreHmOutput$Implausability[,2] < ImThreshold1[1])&(ReIreHmOutput$Implausability[,3] < ImThreshold1[3]) >= 1){
    tmp <- as.matrix(ReIreHmOutput$f_x[(ReIreHmOutput$Implausability[,1] < ImThreshold1[2])&(ReIreHmOutput$Implausability[,2] < ImThreshold1[1])&(ReIreHmOutput$Implausability[,3] < ImThreshold1[3]),])
    
    #tmp <- ReHmOutput$f_x[(ReHmOutput$Implausability[,1] < ImThreshold2[2]),]
    #plot(kdmodel$eval.points , kdmodel$estimate , type ='l')
    for(i in 1:dim(tmp)[1]){
      #lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
      lines(x , FM_EvaluateDenistyEstimate(x , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
    }
    BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[1:500,] , PriorNonImplausibleSetRegular[1:500,] , ReIreHmOutput$NonImplausibleSets[(ReIreHmOutput$Implausability[,1] < ImThreshold1[2])&(ReIreHmOutput$Implausability[,2] < ImThreshold1[1])&(ReIreHmOutput$Implausability[,3] < ImThreshold1[3]) , ] )  
    
  }
  
  dovalidationplots = 1
  source('FM_HistoryMatchMeanPWave.R')
  
}


{
  ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetTotal,
                                                      x = x ,
                                                      F_x =  F_total, 
                                                      f_x = f_total,  
                                                      MD = MD_Total , 
                                                      RRtimes = RRtimes , 
                                                      Corr_sdhat = Corr_sdhat1, 
                                                      imthreshold  = ImThresholdMaxTotal,
                                                      imthreshold2 = ImThresholdMeanTotal)
  
  BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[1:500,] , PriorNonImplausibleSetRegular[1:500,] , ReIreHmOutput$NonImplausibleSets[(ReIreHmOutput$Implausability[,1] < ImThreshold1[2])&(ReIreHmOutput$Implausability[,2] < ImThreshold1[1])&(ReIreHmOutput$Implausability[,3] < ImThreshold1[3]) ,] )  
  
  x11()
  plot(x , predict( kdmodel,x = x) , type ='l')
  if(length(ReIreHmOutput$Implausability) > 4){
    tmp <- as.matrix(ReIreHmOutput$f_x[(ReIreHmOutput$Implausability[,3] < ImThreshold1[3]),])
    
    #tmp <- ReHmOutput$f_x[(ReHmOutput$Implausability[,1] < ImThreshold2[2]),]
    #plot(kdmodel$eval.points , kdmodel$estimate , type ='l')
    for(i in 1:dim(tmp)[1]){
      #lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
      lines(x , FM_EvaluateDenistyEstimate(x , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
    }
  }
}


for(PatientID in listAllPatients[1:length(listAllPatients)] ){
  {
    ECGs <- DP_LoadReducedECGs( path ,  PatientID , FilestoProcess = FilestoProcess  )
    ECGs$ECGII$Value <- 100*ECGs$ECGI$Value
    
    RPeakData <- DP_LoadRpeaksfile( path , PatientID )
    
    dovalidationplots = 0
    
    StartBeat <- 1
    numberofBeats <- 500
    
    RegularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    RegularyIrregularLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    RegularLogical2 <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    RegularyIrregularLogical2 <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
    
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
      timemat[i] <-  RPeakData$RRCombined[rangeofbeats[length(rangeofbeats)],1]
      StartBeatMat[i] <- StartBeat# HM Heart-rhythm
      
      ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, 
                                                       x = xreg ,
                                                       F_x = F_xreg ,
                                                       f_x = f_xreg ,  
                                                       MD = MD_Reg , 
                                                       RRtimes = RRtimes , 
                                                       Corr_sdhat = Corr_sdhat2 , 
                                                       imthreshold = ImThresholdMaxRegular,
                                                       imthreshold2 = ImThresholdMeanRegular)
      ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,
                                                          x = xIrIreg ,
                                                          F_x =  F_x_ReIre,
                                                          f_x = f_x_ReIre,
                                                          MD = MD_ReIre,
                                                          RRtimes = RRtimes ,
                                                          Corr_sdhat = Corr_sdhat1,
                                                          imthreshold  = ImThresholdMaxIrregularlyIrregular,
                                                          imthreshold2 = ImThresholdMeanIrregularlyIrregular)
      ReIreHmTotal <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetTotal,
                                                          x = xTotal ,
                                                          F_x =  F_total, 
                                                          f_x = f_total,  
                                                          MD = MD_Total , 
                                                          RRtimes = RRtimes , 
                                                          Corr_sdhat = Corr_sdhat1, 
                                                          imthreshold  = ImThresholdMaxTotal,
                                                          imthreshold2 = ImThresholdMeanTotal)
      
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
        
        if(sum(NonImplausibleX[,5] != 0) > 0 & sum(NonImplausibleX[,5] == 0) > 0  ){
          if(ImDiff < -0.5){
            RegularyIrregularLogical[i] = 1
            RegularLogical[i] = 0}
        }
        if(sum(NonImplausibleX[,5] != 0) > 0 & sum(NonImplausibleX[,5] == 0) > 0  ){
          if(ImDiff > 0.5){
            
            RegularyIrregularLogical[i] = 0
            RegularLogical[i] = 1
          }
        } 
      }
      if(length(NonImplausibleX) == 5){
        if(sum(NonImplausibleX[5] == 0) > 0 ){
          RegularyIrregularLogical[i] = 1
          disp('Pwave Abscent')
        }
        if(sum(NonImplausibleX[5] != 0) > 0 ){
          RegularLogical[i] = 1
          disp('Pwave Present')
        }
      }
      
      StartBeat <- StartBeat + numberofBeats
      DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
      outputstruct[[i]] <- (list(NonImplausibleX , ReIreHmTotal , cbind(MaxIm[tmplog] , MeanIm[tmplog]) , timemat[i]) )
      
    }
    
    
  }
  
  
  
{
    #RegularyIrregularLogical2 <- ( (meanImRegIre < ImThreshold1[1] )*(minImRegIre <  ImThreshold1[2])*( MulImRegIre <  ImThreshold1[3] ) ) ==1
    #RegularLogical2 <- ( ( (meanImReg < ImThreshold2[1] )*( minImReg <  ImThreshold2[2])*( MulImReg <  ImThreshold2[3] ) ) ==1)
    MetaData = DP_CreateDummyMetaData(PatIndex2017 = PatIndex2017)  
    RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = (RegularLogical*RegularLogical2)==1  , minutethreshold = 1)
    RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = (RegularyIrregularLogical*RegularyIrregularLogical2)==1  , minutethreshold = 1)
    IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = (((RegularyIrregularLogical*RegularyIrregularLogical2) ==0)*((RegularLogical*RegularLogical2)==0)) == 1 , minutethreshold = 1)
    
    RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( PatientID , ' RRTimes' ) )
    RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
    RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
    RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
    #x11(40 , 25)
    #RRPlot + ggtitle('PWaves')
    
    pdf(file = paste0("D:\\nsrdb_outputs\\Annotation" , PatientID , '.pdf') )
    print(RRPlot)
    dev.off()
    
    save(outputstruct , file = paste0("D:\\nsrdb_outputs\\HMOutput" , PatientID , '.RData'))
  }
}
