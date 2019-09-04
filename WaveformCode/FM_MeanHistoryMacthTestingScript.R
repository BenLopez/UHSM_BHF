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
  precomputedfolderpath <- DP_SelectPrecomputedFolder()
  }

UseAnnotatedData <- 1
source( 'FM_CreateRhythumPriors.R' )
source( 'CTEm_LoadDataandCreateEmulatorStructures.R'  )
source( 'CTEm_LoadDataandCreateEmulatorStructuresCDF.R' )
source( 'FM_CreateMeanPWavePriors.R' )

{
  PatientID <- DP_choosepatient( listAllPatients )

  MetaData <- DP_ExtractPatientRecordforIndex( PatIndex2017 = PatIndex2017 , PatientCode = PatientID )
  if(!DP_CheckIfAFPatient(MetaData)){next}
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
  minImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  minImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  meanImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  meanImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  MulImReg <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  MulImRegIre <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  
  SecondWaveLogical <- matrix(0 , ceil(dim(RPeakData$RRCombined)[1]/numberofBeats) , 1)
  
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
    
    if((RegularLogical[i] == 0) & (RegularLogical2[i] ==1) & (RegularyIrregularLogical2[i] == 1) & (RegularyIrregularLogical[i] == 1)){
    SecondWaveLogical[i] <- FM_DubiousCaseLogic(RRtimes[!is.na(RRtimes)], ReIreHmOutput,ReHmOutput,ObservedIm_xx,EmulatorParametersCDFMean, EmulatorParametersCDFMax , RegularLogical)
    disp('AFib Clasified in Dubious Case')
    }
    
    StartBeat <- StartBeat + numberofBeats
    DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
    #outputstruct[[i]] <- (list(NonImplausibleX , ReHmOutput , ReIreHmOutput , cbind(MaxIm[tmplog] , MeanIm[tmplog])) )
    
  }

}


{
  RegularLogical <- (RegularLogical ==1) & (RegularLogical2 ==1) & (RegularyIrregularLogical == 0)&(RegularyIrregularLogical2 == 0)
  IrregularlyIrregularLogical <- (RegularyIrregularLogical == 1)&(RegularyIrregularLogical2 == 1)&(RegularLogical ==0) & (RegularLogical2 ==0)
  #Undecided <- (RegularLogical == 0) & (RegularLogical2 ==1) & (RegularyIrregularLogical2 == 1) & (RegularyIrregularLogical == 1)
  IrregularlyIrregularLogical[SecondWaveLogical == 1] <- 1
  RegularLogical[SecondWaveLogical == 1] <- 0
  
  RegularlyIrregularLogical <- (RegularLogical ==0) &(IrregularlyIrregularLogical ==0)&(Undecided ==0)
  
  
  
  RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = IrregularlyIrregularLogical , minutethreshold = 1)
  RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularLogical , minutethreshold = 1)
  IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularlyIrregularLogical, minutethreshold = 1)
  #UndecidedtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = Undecided, minutethreshold = 1)
  
  
  RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( PatientID , ' RRTimes' ) )
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
  #RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = UndecidedtimePoints , fillcolor = 'purple')
  x11(40 , 25)
  RRPlot + ggtitle('PWaves')
  
  #pdf(file = paste0("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\Annotation" , PatientID , '.pdf') )
  #print(RRPlot)
  #dev.off()
  
  #save(outputstruct , file = paste0("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\HMOutput" , PatientID , '.RData'))
}

{

{StartBeat <- 22001
rangeofbeats <- c(StartBeat: min(dim(RPeakData$RRCombined)[1] , (StartBeat + numberofBeats)) )
RRtimes <-  RPeakData$RRCombined[rangeofbeats,3]
mm <- FM_EmulatorEstimate( Y = RRtimes )
m <- median(RRtimes)
#m <- 0
#mm <- 0
tmp <- RRtimes - m + mm + rnorm(length(RRtimes) , 0 , 0.0025 )
tmp <- tmp[!is.na(tmp)]
kdmodel <- kde( tmp )
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
}
  
x11()
plot(kdmodel$eval.points , kdmodel$estimate , type ='l')

if(length(ReHmOutput$Implausability) > 2){
tmp <- ReHmOutput$f_x[,]
#tmp <- ReHmOutput$f_x[(ReHmOutput$Implausability[,1] < ImThreshold2[2]),]
#plot(x , predict( kdmodel,x = x) , type ='l')
for(i in 1:dim(tmp)[1]){
#lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
  lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
  
}
BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[sample(1:dim(PriorNonImplausibleSetRegularyIreRegular)[1] , 1000 ),] ,
                               PriorNonImplausibleSetRegular[sample(100000:dim(PriorNonImplausibleSetRegular)[1] , 1000 ),] , ReHmOutput$NonImplausibleSets[ , ] )  

}
if(length(ReIreHmOutput$Implausability) > 2){
tmp <- as.matrix(ReIreHmOutput$f_x)
#tmp <- ReHmOutput$f_x[(ReHmOutput$Implausability[,1] < ImThreshold2[2]),]
#plot(kdmodel$eval.points , kdmodel$estimate , type ='l')
for(i in 1:dim(tmp)[1]){
  #lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
  lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
}
BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[sample(1:dim(PriorNonImplausibleSetRegularyIreRegular)[1] , 1000 ),] ,
                               PriorNonImplausibleSetRegular[sample(100000:dim(PriorNonImplausibleSetRegular)[1] , 1000 ),] ,
                               ReIreHmOutput$NonImplausibleSets )  

}

dovalidationplots = 1
source('FM_HistoryMatchMeanPWave.R')


if( is.null(ReIreHmOutput$MinImplausibility) & is.null(ReHmOutput$MinImplausibility) ){

ReMeanProb <- FM_EvaluateObservedImplausibilityEmulator(Xstar =  ReHmOutput$NonImplausibleSets[,-3] ,
                                                        Im =  ReHmOutput$Implausability[,2]  ,
                                                        EmulatorParameters = EmulatorParametersCDFMean ,
                                                        ObservedIm_xx = ObservedIm_xx)
IrReMeanProb <- FM_EvaluateObservedImplausibilityEmulator(Xstar =  ReIreHmOutput$NonImplausibleSets[,-3] ,
                                                          Im =  ReIreHmOutput$Implausability[,2]  ,
                                                          EmulatorParameters = EmulatorParametersCDFMean ,
                                                          ObservedIm_xx = ObservedIm_xx)
ReMaxProb <- FM_EvaluateObservedImplausibilityEmulator(Xstar =  ReHmOutput$NonImplausibleSets[,-3] ,
                                                       Im =  ReHmOutput$Implausability[,1]  , 
                                                       EmulatorParameters = EmulatorParametersCDFMax ,
                                                       ObservedIm_xx = ObservedIm_xx)
IrReMaxProb <- FM_EvaluateObservedImplausibilityEmulator(Xstar =  ReIreHmOutput$NonImplausibleSets[,-3] ,
                                                          Im =  ReIreHmOutput$Implausability[,1]  ,
                                                          EmulatorParameters = EmulatorParametersCDFMax ,
                                                          ObservedIm_xx = ObservedIm_xx)

x11()
plot( IrReMaxProb , IrReMeanProb , col = 'red' , xlab = 'Max Im Prob' , ylab = 'Mean Im Prob' , title = 'Second Wave Comparison' , xlim = c(0,1) , ylim =c(0,1))
points( ReMaxProb , ReMeanProb )


}

}


{
ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetTotal,
                                                    x = xTotal ,
                                                    F_x =  F_total, 
                                                    f_x = f_total,  
                                                    MD = MD_Total , 
                                                    RRtimes = RRtimes , 
                                                    Corr_sdhat = Corr_sdhat1, 
                                                    imthreshold  = ImThresholdMaxTotal , 
                                                    imthreshold2 = ImThresholdMeanTotal )
if(is.null(ReIreHmOutput$NonImplausibleSets)){
BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[1:500,] , PriorNonImplausibleSetRegular[1:500,] , t(matrix(ReIreHmOutput$MinImplausiblePoint , 10, 100)) )  
BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[1:500,] , PriorNonImplausibleSetRegular[1:500,] , t(matrix(ReHmOutput$MinImplausiblePoint , 10, 100)) )  
    }else{
BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[sample(1:dim(PriorNonImplausibleSetRegularyIreRegular)[1] , 1000 ),] ,
                               PriorNonImplausibleSetRegular[sample(100000:dim(PriorNonImplausibleSetRegular)[1] , 1000 ),], ReIreHmOutput$NonImplausibleSets )  
    }

x11()
plot(x , predict( kdmodel,x = x) , type ='l')
if(is.null(ReIreHmOutput$NonImplausibleSets)){

if(length(ReIreHmOutput$Implausability > 3)){
  lines(x , FM_EvaluateDenistyEstimate(x , ReIreHmOutput$MinImplausiblePoint) , col = rgb(0,0,1 , alpha = 0.5) )
  lines(x , FM_EvaluateDenistyEstimate(x , ReHmOutput$MinImplausiblePoint) , col = rgb(1,0,0 , alpha = 0.5) )
  
}
  }
  
if( !is.null(ReIreHmOutput$NonImplausibleSets)> 4 ){
  tmp <- as.matrix(ReIreHmOutput$f_x[,])
  
  #tmp <- ReHmOutput$f_x[(ReHmOutput$Implausability[,1] < ImThreshold2[2]),]
  #plot(kdmodel$eval.points , kdmodel$estimate , type ='l')
  for(i in 1:dim(tmp)[1]){
    #lines(kdmodel$eval.points , FM_EvaluateDenistyEstimate(kdmodel$eval.points  , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
    lines(x , FM_EvaluateDenistyEstimate(x , ReIreHmOutput$NonImplausibleSets[i,]) , col = rgb(0,0,1 , alpha = 0.1) )
  }
}
}


for(PatientID in listAllPatients[125:length(listAllPatients)] ){
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
      
      RRtimes <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,3]
      timemat[i] <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats[length(rangeofbeats)],1]
      StartBeatMat[i] <- StartBeat
      
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
      
      if((RegularLogical[i] == 0) & (RegularLogical2[i] ==1) & (RegularyIrregularLogical2[i] == 1) & (RegularyIrregularLogical[i] == 1)){
        SecondWaveLogical[i] <- FM_DubiousCaseLogic(RRtimes[!is.na(RRtimes)], ReIreHmOutput,ReHmOutput,ObservedIm_xx,EmulatorParametersCDFMean, EmulatorParametersCDFMax , RegularLogical)
        disp('AFib Clasified in Dubious Case')
      }
      
      StartBeat <- StartBeat + numberofBeats
      DP_WaitBar(i / ceil(dim(RPeakData$RRCombined)[1]/numberofBeats))
      outputstruct[[i]] <- (list(NonImplausibleX , ReHmOutput , ReIreHmOutput ,TotalHmOutput , timemat[i]) )
      
    }
    
    
  }
  
  
  
{
  RegularLogical <- (RegularLogical ==1) & (RegularLogical2 ==1) & (RegularyIrregularLogical == 0)&(RegularyIrregularLogical2 == 0)
  IrregularlyIrregularLogical <- (RegularyIrregularLogical == 1)&(RegularyIrregularLogical2 == 1)&(RegularLogical ==0) & (RegularLogical2 ==0)
  #Undecided <- (RegularLogical == 0) & (RegularLogical2 ==1) & (RegularyIrregularLogical2 == 1) & (RegularyIrregularLogical == 1)
  IrregularlyIrregularLogical[SecondWaveLogical == 1] <- 1
  RegularLogical[SecondWaveLogical == 1] <- 0
  
  RegularlyIrregularLogical <- (RegularLogical ==0) &(IrregularlyIrregularLogical ==0)
  
  RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = IrregularlyIrregularLogical , minutethreshold = 1)
  RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularLogical , minutethreshold = 1)
  IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularlyIrregularLogical, minutethreshold = 1)
  #UndecidedtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = Undecided, minutethreshold = 1)
  
  
  RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( PatientID , ' RRTimes' ) )
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
  #RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = UndecidedtimePoints , fillcolor = 'purple')
  pdf(file = paste0("D:\\AllPatientAnnotations\\Annotation" , PatientID , '.pdf') )
    print(RRPlot)
    dev.off()
    
  outputstruct[[i+1]] <- cbind(RegularLogical, RegularLogical2, RegularyIrregularLogical,RegularyIrregularLogical2,SecondWaveLogical)
    
  save(outputstruct , file = paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData'))
  }
}
