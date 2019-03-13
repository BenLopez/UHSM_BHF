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

PsimulatorFunction <- function( x , t_observation ){
  
  Pcen <- x[1]
  Pwidth <- x[2]
  Pcen2 <- x[3]
  Pwidth2 <- x[4]
  
  return( ECGSim_Gaussian( x = t_observation , mu = Pcen , sigma = Pwidth ) + ECGSim_Gaussian( x = t_observation , mu = Pcen2 , sigma = Pwidth2 ) )
}
#BC_PlotPairs(PriorNonImplausibleSet)
ModelDiscrepancy <- function(x , t_observation, PsimulatorFunction){
  1/( PsimulatorFunction( x , t_observation ) + 0.5 )
}
CalculateImplausability <- function( t_observation , x ,  z  , mdfun = ModelDiscrepancy){
  H = PWaveHM_CreateDesignMatrix(t_observation , x , PsimulatorFunction)
  Beta = PWaveHM_CalculateBetas(H , z)
  Im = mean( abs(z - H%*%Beta) / sqrt(ModelDiscrepancy(x , t_observation, PsimulatorFunction)))
  return( Im )
} 
PWaveHM_CreateDesignMatrix <- function(t_observation , x , PsimulatorFunction){
  t_observation <- as.matrix(t_observation)
  return( cbind( matrix(1 , dim(t_observation)[1] , 1 ) , t_observation , t_observation^2, as.matrix(PsimulatorFunction(x , t_observation))  ) )
}


EmulatorParameters <- PWaveHM_CreateDefaultEmulationclass()
FilestoProcess = DP_ChooseECGstoProcess()

{Xstar = seq(0.5 ,1 , 0.01)
  PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.001 , 0.95, 0.001 ) , b = c(0.6  , 0.04 , 0.6 , 0.04 ) , numbersamples = 1000000 )
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.1 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.02 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]
  QSwidth = 10 }


PatientID <- DP_choosepatient(listAllPatients)

MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = PatientID)
ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)

RPeakTs <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData)

StartBeat <- 1
for(i in 1:72){
  rangeofbeats <- c(StartBeat:(StartBeat + 500))
  
  RRtimes <- PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats,3]
  RPeakKDEEstmate <- kde( RRtimes )
  
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
  
      RPeakdisPlot <- RPeakdisPlot + ggtitle(ggtitle(paste0( ' RR-Times Distribution ')))
      
  #QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[10000:10500,] , QSwidth = 10)
  
  PWavePlot1 <- FM_PlotPWaveAnalysisSingle(ECG = ECGs$ECGI , Beats = RPeakData$RRCombined[rangeofbeats,] , QSwidth  = 10) 
  PWavePlot2 <- FM_PlotPWaveAnalysisSingle(ECG = ECGs$ECGII , Beats = RPeakData$RRCombined[rangeofbeats,] , QSwidth  = 10)
  PWavePlot3 <- FM_PlotPWaveAnalysisSingle(ECG = ECGs$ECGIII , Beats = RPeakData$RRCombined[rangeofbeats,] , QSwidth  = 10)
  lay <- rbind(c(1,1 ) ,
               c(1,1 ) ,
               c(1,1 ) ,
               c(2,2 ) , 
               c(3,3 ) ,
               c(4,4 )) 
  lay <- rbind( c( 1 , 1 , 1 , 1) ,
                c( 1 , 1 , 1 , 1) ,
                c( 1 , 1 , 1 , 1) ,
                c( 2 , 2 , 3 , 3) , 
                c( 2 , 2 , 4 , 4) ,
                c( 2 , 2 , 5 , 5) )  
  
  x11(40 , 25)
  tmpRPeakTs <- RPeakTs + geom_vline(xintercept = PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats[1],1] , linetype = "dashed") + geom_vline(xintercept = PE_CleanRpeaks( RPeakData$RRCombined )[rangeofbeats[length(rangeofbeats)],1] , linetype = "dashed") 
  print(grid.arrange(tmpRPeakTs , RPeakdisPlot , PWavePlot1 +ggtitle('P-Waves') ,PWavePlot2 , PWavePlot3 , layout_matrix = lay))
  dev.copy(png , paste0('C:\\Users\\Ben\\Documents\\Output Images\\Videos\\image' , rangeofbeats[1] , '.png'))
  dev.off()
  StartBeat = StartBeat + 500 +1
}


