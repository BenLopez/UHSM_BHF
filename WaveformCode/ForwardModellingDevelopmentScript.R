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

# prior specification regular 

source('FM_CreateRhythumPriors.R')

{Goodbad <- matrix(T , dim(PriorNonImplausibleSetRegularyIreRegular)[1])

for(samplenum in 1:100){

x11()
print(plot(  x  ,  F_x_ReIre[ samplenum , ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )) 
#print(lines( x  ,  F_x_ReIre[ samplenum , ] + 2*(MD_ReIre[ samplenum , ])  , col ='red' , lty=2 ) )
#print(lines( x  ,  F_x_ReIre[ samplenum , ] - 2*(MD_ReIre[ samplenum , ])  , col ='red' , lty=2 ) )
print(abline(0.6 , 0))
Sys.sleep(0)
PriorNonImplausibleSetRegularyIreRegular[samplenum, ]

UserResponse <- winDialog(type = c("yesnocancel") , message = 'Is this sample non-implausible?')

if(UserResponse == 'YES'){
  Goodbad[samplenum] <- T
}
if(UserResponse == 'NO'){
  Goodbad[samplenum] <- F
}
if(UserResponse == 'CANCEL'){
  break
}
dev.off()
}
}


BC_PlotCompareTwoHists(PriorNonImplausibleSetRegularyIreRegular[1:1000 , 1:9] , PriorNonImplausibleSetRegular[1:1000 , 1:9])
BC_PlotPairsFromTwoVariables(PriorNonImplausibleSetRegularyIreRegular[1:1000 , 1:9] , PriorNonImplausibleSetRegular[1:1000 , 1:9] , alpha = 0.1)

PatientID <- DP_choosepatient(listAllPatients)

MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = PatientID)
ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)

{
#rangeofbeats <- 1:500
rangeofbeats <- 30000:(30000 + 500)
#rangeofbeats <- 25000:25500
  
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

#QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[rangeofbeats,] , QSwidth = 10)
#EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct , EmulatorParameters = EmulatorParameters , Xstar )

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

{x11(20 , 14)
print(grid.arrange(RPeakdisPlot , PWavePlot1 ,PWavePlot2 , PWavePlot3 , layout_matrix = lay))
}
}

ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , F_x = F_xreg ,f_x = f_xreg ,  MD = MD_Reg , RRtimes = RRtimes , imthreshold = ImThreshold2[2] )
ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,x = x ,F_x =  F_x_ReIre, f_x = f_x_ReIre, MD = MD_ReIre, RRtimes = RRtimes, imthreshold =  ImThreshold1[2] )

if( length(ReHmOutput$Implausability) > 3 ){
  
samplenum <- which.min(ReHmOutput$Implausability[,1])
plot(  x , ReHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density' , xlim =c(0.5,0.75) )
lines( x , ReHmOutput$y , type ='l' , col = 'black' )
lines( x , ReHmOutput$f_x[samplenum , ] + 4*(MDFunction(ReHmOutput$f_x[ samplenum , ] , ReHmOutput$NonImplausibleSets[samplenum , ] ))  , col ='red' , lty=2 )
lines( x , ReHmOutput$f_x[samplenum , ] - 4*(MDFunction(ReHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  ))  , col ='red' , lty=2 )

x11()
p1 <- ggplot(data.frame(x = x , y = ReHmOutput$f_x[samplenum, ]) , aes(x , y) ) +
  geom_line(col ='blue') +
  xlim(0.5,0.75) +
  geom_line(data = data.frame(x = x , y = ReHmOutput$y) , aes(x , y)  , col ='black')+
  geom_line(data = data.frame(x = x , y = ReHmOutput$f_x[samplenum , ] + 4*(MDFunction(ReHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  ))) , aes(x , y)  , col ='red')+
  geom_line(data = data.frame(x = x , y = ReHmOutput$f_x[samplenum , ] - 4*(MDFunction(ReHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  ))) , aes(x , y)  , col ='red')+
  ggtitle('Acceptable Match') +
  xlab('RR')+
  ylab('Culmulative Probability')
print(p1)


}else{
  
plot(  x , ReHmOutput$f_x , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
lines( x , ReHmOutput$y , type ='l' , col = 'black' )
lines( x , ReHmOutput$f_x + 3*(MDFunction(ReHmOutput$f_x , ReHmOutput$NonImplausibleSets ) + 0.0045 )  , col ='red' , lty=2 )
lines( x , ReHmOutput$f_x - 3*(MDFunction(ReHmOutput$f_x  , ReHmOutput$NonImplausibleSets  ) + 0.0045 )  , col ='red' , lty=2 )

}

if( length( ReIreHmOutput$Implausability ) > 3 ){
  
  samplenum <- which.min(ReIreHmOutput$Implausability[,1])
  plot(  x , ReIreHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReIreHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReIreHmOutput$f_x[samplenum , ] + ImThreshold1[2]*(MDFunction(ReIreHmOutput$f_x[ samplenum , ] , ReIreHmOutput$NonImplausibleSets[samplenum , ] ) + 0.0045)  , col ='red' , lty=2 )
  lines( x , ReIreHmOutput$f_x[samplenum , ] - ImThreshold1[2]*(MDFunction(ReIreHmOutput$f_x[ samplenum , ]  ,ReIreHmOutput$NonImplausibleSets[samplenum , ]  )+ 0.0045)  , col ='red' , lty=2 )

  x11()
  p1 <- ggplot(data.frame(x = x , y = ReIreHmOutput$f_x[samplenum, ]) , aes(x , y) ) +
    geom_line(col ='blue') +
    xlim(0.1,0.75)+
    geom_line(data = data.frame(x = x , y = ReIreHmOutput$y) , aes(x , y)  , col ='black')+
    geom_line(data = data.frame(x = x , y = ReIreHmOutput$f_x[samplenum , ] + 4*(MDFunction(ReIreHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  ))) , aes(x , y)  , col ='red')+
    geom_line(data = data.frame(x = x , y = ReIreHmOutput$f_x[samplenum , ] - 4*(MDFunction(ReIreHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  ))) , aes(x , y)  , col ='red')+
    ggtitle('Acceptable Match') +
    xlab('RR')+
    ylab('Culmulative Probability')
  print(p1)  
  
}else{
  
  plot(  x , ReIreHmOutput$f_x , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReIreHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReIreHmOutput$f_x +  ImThreshold1[2]*( MDFunction2(ReIreHmOutput$f_x ,ReHmOutput$NonImplausibleSets ) )  , col ='red' , lty=2 )
  lines( x , ReIreHmOutput$f_x -  ImThreshold1[2]*( MDFunction2(ReIreHmOutput$f_x  ,ReHmOutput$NonImplausibleSets  ) )  , col ='red' , lty=2 )

  
}

BC_PlotPairsFromTwoVariables( PriorNonImplausibleSetRegularyIreRegular[1:5000,1:9],PriorNonImplausibleSetRegular[1:1000,1:9] , alpha = 0.02)
BC_PlotPairsFromTwoVariables( PriorNonImplausibleSetRegularyIreRegular[1:1000,1:9],ReIreHmOutput$NonImplausibleSets[,1:9] , alpha = 0.03)


###### Test with known ground truth ######

RRtimes <- FM_SampleGMM(X = PriorNonImplausibleSetRegularyIreRegular[1 , ] ,  N = 231)

ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , F_x = F_xreg ,f_x = f_xreg ,  MD = MD_Reg , RRtimes = RRtimes , imthreshold = 3 )
ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,x = x ,F_x =  F_x_ReIre, f_x = f_x_ReIre, MD = MD_ReIre, RRtimes = RRtimes, imthreshold = 3 )

if(length(ReHmOutput$Implausability) > 4){
  samplenum <- which.min(ReHmOutput$Implausability[,4])
  plot(  x , ReHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReHmOutput$f_x[samplenum , ] + 2*(MDFunction(ReHmOutput$f_x[ samplenum , ] ,ReHmOutput$NonImplausibleSets[samplenum , ] ))  , col ='red' , lty=2 )
  lines( x , ReHmOutput$f_x[samplenum , ] - 2*(MDFunction(ReHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  ))  , col ='red' , lty=2 )
}else{
  plot(  x , ReHmOutput$f_x , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReHmOutput$f_x + 2*(MDFunction(ReHmOutput$f_x ,ReHmOutput$NonImplausibleSets ))  , col ='red' , lty=2 )
  lines( x , ReHmOutput$f_x - 2*(MDFunction(ReHmOutput$f_x  ,ReHmOutput$NonImplausibleSets  ))  , col ='red' , lty=2 )
}

if( length( ReIreHmOutput$Implausability[,1] ) > 3 ){
  
  samplenum <- which.max(ReIreHmOutput$Implausability[,4])
  plot(  x , ReIreHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReIreHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReIreHmOutput$f_x[samplenum , ] + 3*(MDFunction(ReIreHmOutput$f_x[ samplenum , ] , ReIreHmOutput$NonImplausibleSets[samplenum , ] ))  , col ='red' , lty=2 )
  lines( x , ReIreHmOutput$f_x[samplenum , ] - 3*(MDFunction(ReIreHmOutput$f_x[ samplenum , ]  ,ReIreHmOutput$NonImplausibleSets[samplenum , ]  ))  , col ='red' , lty=2 )
  
}else{
  
  plot(  x , ReIreHmOutput$f_x , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReIreHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReIreHmOutput$f_x + 3*(MDFunction2(ReIreHmOutput$f_x ,ReHmOutput$NonImplausibleSets ))  , col ='red' , lty=2 )
  lines( x , ReIreHmOutput$f_x - 3*(MDFunction2(ReIreHmOutput$f_x  ,ReHmOutput$NonImplausibleSets  ))  , col ='red' , lty=2 )
 
}

FM_ExtractActiveOutputs <- function(y , x){
  
  NonZeroLogical <- c( ((y > 0.03)*(y < 0.97))==1 ) ==1
  
  if(sum(NonZeroLogical) <= 2){
    NonZeroLogical <- c( ((y > 0.01)*(y < 0.99))==1 ) ==1
  }
  if(sum(NonZeroLogical) <= 2){
    NonZeroLogical <- c( ((y > 0.001)*(y < 0.999))==1 ) ==1
  }
  
  seq(min(x[NonZeroLogical]) , max(x[NonZeroLogical]) , (max(x[NonZeroLogical]) - min(x[NonZeroLogical]))/length(x) )
  return(seq(min(x[NonZeroLogical]) , max(x[NonZeroLogical]) , (max(x[NonZeroLogical]) - min(x[NonZeroLogical]))/length(x) )[1:201])
  
}

# Cabage code for fitting with pdf rather than pdf
#ReHmOutput <- FM_HistoryMatchRRDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , f_x = f_xreg , MD = MD_Reg , RRtimes = RRtimes , imthreshold = 0.5)
#ReIreHmOutput <- FM_HistoryMatchRRDensity(PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,x = x ,f_x =  f_x_ReIre, MD = MD_ReIre, RRtimes = RRtimes, imthreshold = 0.5)

#samplenum <- which.min(ReHmOutput$Implausability)
#plot(  x , ReHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
#lines( x , ReHmOutput$y , type ='l' , col = 'black' )
#lines( x , ReHmOutput$f_x[samplenum , ] + 2*(MDFunction(ReHmOutput$f_x[ samplenum , ] ,ReHmOutput$NonImplausibleSets[samplenum , ] ))  , col ='red' , lty=2 )
#lines( x , ReHmOutput$f_x[samplenum , ] - 2*(MDFunction(ReHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  ))  , col ='red' , lty=2 )

#samplenum <- which.min(ReIreHmOutput$Implausability)
#plot( x ,  ReIreHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'   , ylim = c(0 , 5))   
#lines( x , ReIreHmOutput$y , type ='l' , col = 'black' )
#lines( x , ReIreHmOutput$f_x[samplenum , ] + 3*(MDFunction2(ReIreHmOutput$f_x[ samplenum , ] , ReIreHmOutput$NonImplausibleSets[samplenum , ] ))  , col ='red' , lty=2 )
#lines( x , ReIreHmOutput$f_x[samplenum , ] - 3*(MDFunction2(ReIreHmOutput$f_x[ samplenum , ] , ReIreHmOutput$NonImplausibleSets[samplenum , ]  ))  , col ='red' , lty=2 )

#X <- PriorNonImplausibleSetRegularyIreRegular[1:2000,]
#Z <- ReIreHmOutput$NonImplausibleSets[ReIreHmOutput$Implausability[,2]<ImThreshold1[1] , ]
#Y <- PriorNonImplausibleSetRegular[1:1000,]

#x11(20 , 14)
#pairs( rbind(X ,  Y , Z) ,  col = rbind(as.matrix(rep(rgb(1,0,0 , alpha = 0.05) , size(X)[1])) , as.matrix(rep(rgb(0,0,1 , alpha = 0.05) , size(Y)[1]) ) ,  as.matrix(rep(rgb(0,1,0 , alpha = 0.15) , size(Z)[1]) )  ) , pch = 16) 
