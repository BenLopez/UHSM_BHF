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


source('FM_CreatePWavePriors.R')
source('FM_CreateRhythumPriors.R')

PatientID <- DP_choosepatient(listAllPatients)

MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = PatientID)
ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)


{
  #rangeofbeats <- 1:500
  rangeofbeats <- 9500:(9500 + 500)

  RRtimes <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,3]
  tmp <- RRtimes + rnorm(length(RRtimes) , 0 , 0.0025) - rollmedian(RRtimes , k = 21 , na.pad = T) + median(RRtimes)
  RPeakKDEEstmate <- kde( as.matrix(tmp[!is.na(tmp)]) )
  
  MaxDensity <- max(RPeakKDEEstmate$estimate)
  modeRate <- RPeakKDEEstmate$eval.points[which.max(RPeakKDEEstmate$estimate)]
  Peakslogical <- PE_FindLocalTurningPoints(RPeakKDEEstmate$estimate/max(RPeakKDEEstmate$estimate) > 0.075 , RPeakKDEEstmate$estimate)
  numberofpeaks <- sum(Peakslogical)
  PeakLocations <- RPeakKDEEstmate$eval.points[Peakslogical]
  PeakHieghts <- RPeakKDEEstmate$estimate[Peakslogical]
  NumberTachicardicBeats <- sum(RRtimes < 0.6)
  
  PrecentageforRegular <- 0.125
  ProportionRegular <- sum((RRtimes > (modeRate - PrecentageforRegular*modeRate))*(RRtimes < (modeRate + PrecentageforRegular*modeRate))) / length(RRtimes) 
  
  RPeakdisPlot <- ggplot(data = data.frame( RRTimes=RPeakKDEEstmate$eval.points , density = RPeakKDEEstmate$estimate) , aes(RRTimes , density)) + 
    geom_line( col ='blue') +
    geom_vline(xintercept =  0.6) +
    geom_vline(xintercept =  1 ) +
    geom_vline(xintercept =  modeRate  , col = 'red' , linetype = "dashed" ) +
    geom_vline(xintercept =  modeRate + PrecentageforRegular*modeRate  , col = 'red' , linetype = "dashed" ) +
    geom_vline(xintercept =  modeRate - PrecentageforRegular*modeRate  , col = 'red' , linetype = "dashed" ) +
    xlim(0.1 , 1.2) 
  
  RPeakdisPlot <- RPeakdisPlot + ggtitle(ggtitle(paste0(PatientID , ' RR-Times Distribution')))
  
  # Rhythm History Matching
  ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , F_x = F_xreg ,f_x = f_xreg ,  MD = MD_Reg , RRtimes = RRtimes , imthreshold = ImThreshold2[2] )
  ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,x = x ,F_x =  F_x_ReIre, f_x = f_x_ReIre, MD = MD_ReIre, RRtimes = RRtimes, imthreshold =  ImThreshold1[2] )
  
  QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = 10)
  EmulatedQS <- FMPWaveHM_EmulateTQSegment( QS_Struct , EmulatorParameters = EmulatorParameters , Xstar )
  
  source('FM_HistoryMatchPWaves.R')
  
  d <- FM_CreateDefaultFeatureTable( )
  
  RegularLogical <- sum(ReHmOutput$Implausability[ , 1] < ImThreshold2[2])* sum(ReHmOutput$Implausability[ , 2] < ImThreshold2[1])*sum(ReHmOutput$Implausability[ , 3] < ImThreshold2[3])>0
  IregularIrregularlogical <- sum((ReIreHmOutput$Implausability[ , 1] < ImThreshold1[2])*sum(ReIreHmOutput$Implausability[ , 2] < ImThreshold1[1])*sum(ReIreHmOutput$Implausability[ , 3] < ImThreshold1[3]))>0
  
  d[1 , 4] <- signif(mean(60/tmp , na.rm =T) , 2)
  d[1 , 5] <- signif(var(60/tmp, na.rm =T) , 2)
  d[2 , 4] <- signif(E_PA, 2)
  d[2 , 5] <- signif(V_PA, 2)
  d[3 , 4] <- signif(E_PWaveDurations, 2)
  d[3 , 5] <- signif(V_PWaveDurations, 2)
  d[4 , 4] <- signif(E_PRInterval, 2)
  d[4 , 5] <- signif(V_PRInterval, 2)
  
  if(abs(E_PA) < 4){
    d[2 , 2] <- 'Abnormal (abscent)'
  }
  d[2 , 3] <- paste0( signif((sum(apply(Implausability , 2 , min) <0.6)/500)*100 , 3) , '%')
  
  if(RegularLogical){
    d[1 , 2] <- 'Regular'
    d[1 , 3] <- 'Yes'
  }
  if(IregularIrregularlogical){
    d[1 , 2] <- 'Irregularly Irregular'
    d[1 , 3] <- 'Yes'
  }
  if(!IregularIrregularlogical && !RegularLogical){
    d[1 , 2] <- 'Irregular'
    d[1 , 3] <- 'Na'
  }
  
  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
  tbl <- tableGrob(d, rows=NULL, theme=tt)
  
  PWavePlot1 <- FM_PlotPWaveAnalysisSingle(ECG = ECGs$ECGI , Beats = RPeakData$RRCombined[rangeofbeats,] , QSwidth  = 10)
  PWavePlot2 <- FM_PlotPWaveAnalysisSingle(ECG = ECGs$ECGII , Beats = RPeakData$RRCombined[rangeofbeats,] , QSwidth  = 10)
  PWavePlot3 <- FM_PlotPWaveAnalysisSingle(ECG = ECGs$ECGIII , Beats = RPeakData$RRCombined[rangeofbeats,] , QSwidth  = 10)
  lay <- rbind(c(1,1 ) ,
               c(1,1 ) ,
               c(1,1 ) ,
               c(2,2 ) , 
               c(3,3 ) ,
               c(4,4 )) 
  lay <- rbind(c(1,1 , 5 , 5) ,
               c(1,1 , 5 , 5) ,
               c(1,1 , 5 , 5) ,
               c(2,2 , 6 , 6) , 
               c(3,3 , 6 , 6) ,
               c(4,4 , 6 , 6))  
  
  {
    x11(25 , 14)
    print(grid.arrange(RPeakdisPlot , PWavePlot1 ,PWavePlot2 , PWavePlot3  ,tbl, layout_matrix = lay))
  }
}


mQS <- apply( EmulatedQS , 2 , mean )
VQS <- apply( EmulatedQS , 2 , var )

HMOutput <- PWaveHM_EmulateEstimatePAmplitude(QS_Struct , EmulatorParameters, Xstar , PriorNonImplausibleSet , Graphics = 0)

{z <- EmulatedQS[22,] 
Implausability2 <- PWaveHM_CalulateImplausabilityTQSegment(Xstar , z , PriorNonImplausibleSet)

H = PWaveHM_CreateDesignMatrix(Xstar , PriorNonImplausibleSet[which.min(Implausability2) , ] , PsimulatorFunction)

alpha <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + diag(V_me))
AdjustedBeta <- FMPWaveHM_BLUBetas( H , z = z  , E_Beta , V_Beta  , V_me  , alpha = 0)

plot(Xstar , z , type = 'l' , col ='red' , ylim = c(-40 , 40))
lines(Xstar , H[ , 1:4]%*%as.matrix(AdjustedBeta$E_z_Beta[1:4,])  )
lines(Xstar ,  H[ , 1:3]%*%as.matrix(AdjustedBeta$E_z_Beta[1:3,]) , type = 'l' , col ='green' )}


# Create Implausability Thresholds for PWaveHistoryMatching

samplex <- PriorNonImplausibleSet[sample(1:dim(PriorNonImplausibleSet)[1] , 1) , ]
H = PWaveHM_CreateDesignMatrix(Xstar , samplex , PsimulatorFunction)
V_md <- ModelDiscrepancy(Xstar , x = samplex , PsimulatorFunction )
C_md <- CF_ExponentialFamily(Xstar , Xstar , 0.2 , 2)

mdSample <- t(rmvnorm(1 , 0*V_md , (sqrt(V_md)%*%t(sqrt(V_md)))*C_md  ))
plot(mdSample  , type ='l')

samplex <- PriorNonImplausibleSet[which.min(Implausability) , ] 

Samplesofz <- matrix(0, dim(PriorNonImplausibleSet)[1] , dim(H )[1]  )
SamplesofIm <- matrix(0, dim(PriorNonImplausibleSet)[1] , 1  )
for( i in 1:dim(PriorNonImplausibleSet)[1]){
  samplex <- PriorNonImplausibleSet[i,]
  H = PWaveHM_CreateDesignMatrix(Xstar , samplex , PsimulatorFunction)
  BetaSample <- t(rmvnorm(1 , E_Beta , V_Beta ))
  V_md <- ModelDiscrepancy(Xstar , x = samplex , PsimulatorFunction )
  samplesofz[i,] <- H[ , 1:4]%*%BetaSample[1:4] +  t(rmvnorm(1 , 0*V_md , (sqrt(4*V_md)%*%t(sqrt(4*V_md)))*C_md  ))  + rnorm(length(Xstar) , 0 , 1 )
  SamplesofIm[i,] <- CalculateImplausability(Xstar , samplex ,  samplesofz[i,])
}
