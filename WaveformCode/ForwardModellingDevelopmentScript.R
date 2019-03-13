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
{
  {
    numberofsamples <- 100000
    PriorNonImplausibleSetRegular <- matrix(0 , numberofsamples , 10)
    
    PriorNonImplausibleSetRegular[, 5] <- runif(numberofsamples , 0.6 , 1)
    PriorNonImplausibleSetRegular[, 4] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] - 0.15*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] - 0.05*PriorNonImplausibleSetRegular[, 5])
    PriorNonImplausibleSetRegular[, 6] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] + 0.05*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] + 0.15*PriorNonImplausibleSetRegular[, 5] )
    
    PriorNonImplausibleSetRegular[, 1] <- runif(numberofsamples , 0 , 0.1)
    PriorNonImplausibleSetRegular[, 3] <- runif(numberofsamples , 0 , 0.1)
    PriorNonImplausibleSetRegular[, 2] <- 1 - (PriorNonImplausibleSetRegular[, 1] + PriorNonImplausibleSetRegular[, 3]) 
    
    PriorNonImplausibleSetRegular[, 7] <- runif(numberofsamples , 0.004 , (0.15*PriorNonImplausibleSetRegular[, 4]) / 2.5 )
    PriorNonImplausibleSetRegular[, 8] <- runif(numberofsamples , 0.001 , (0.15*PriorNonImplausibleSetRegular[, 5]) / 2.5 )
    PriorNonImplausibleSetRegular[, 9] <- runif(numberofsamples ,  apply(as.matrix(PriorNonImplausibleSetRegular[, 7]- 0.15*PriorNonImplausibleSetRegular[, 7]) , 1 , function(X){max(0,X)} ) , PriorNonImplausibleSetRegular[, 7] + 0.15*PriorNonImplausibleSetRegular[, 7] )
    
    PriorNonImplausibleSetRegular[, 10]<- runif(numberofsamples ,0,1)
    
    MDFunction <- function(F_x , X){
      output <- sqrt( (F_x*(1-F_x))/231 + 0.000000001 )
      output[output <= 0.0001] <- quantile(output[output > 0] , 0.05, na.rm = T)
      return(output)
    }
    
    #MDFunction <- function(f_x , X){
    #  return( 0.3  )
    #}
    
    x <- seq(0 , 2, 0.01)
    f_xreg <- t(apply(PriorNonImplausibleSetRegular , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
    F_xreg <- t(apply(f_xreg , 1 , function(X){cumsum( c(0,diff(x))*(X/sum(c(0,diff(x))*(X))) ) }))
    
    MD_Reg <- matrix(0 , numberofsamples , length(x))
    for( i in 1:numberofsamples){
      MD_Reg[i , ] <- MDFunction( F_xreg[i , ] , PriorNonImplausibleSetRegular[i , ] )
    }
  }
  
  # Prior specification regularly-irregular
  {
    PriorNonImplausibleSetRegularyIreRegular <- matrix(0 , numberofsamples , 10)
    
    PriorNonImplausibleSetRegularyIreRegular[, 5] <- runif(numberofsamples , 0.3 , 1)
    PriorNonImplausibleSetRegularyIreRegular[, 4] <- runif(numberofsamples , 0.25 ,  apply(as.matrix(0.85*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 , function(X){max(X , 0.26)}) )
    PriorNonImplausibleSetRegularyIreRegular[, 6] <- runif(numberofsamples , apply(as.matrix(PriorNonImplausibleSetRegularyIreRegular[, 5] + 0.1*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 ,function(X){max(0.1 , X)} ) , 2)
    
    
    PriorNonImplausibleSetRegularyIreRegular[, 2] <- runif(numberofsamples , 0.4 , 0.8)
    PriorNonImplausibleSetRegularyIreRegular[, 1] <- runif(numberofsamples , rep(0 , numberofsamples ) ,  apply( cbind(PriorNonImplausibleSetRegularyIreRegular[, 2] , (1 - PriorNonImplausibleSetRegularyIreRegular[,2])) , 1 , min) )
    PriorNonImplausibleSetRegularyIreRegular[, 3] <- 1 - PriorNonImplausibleSetRegularyIreRegular[, 1] - PriorNonImplausibleSetRegularyIreRegular[, 2]
    
    PriorNonImplausibleSetRegularyIreRegular[, 8] <- runif(numberofsamples , 0.05 , 0.2 )
    PriorNonImplausibleSetRegularyIreRegular[, 7] <- runif(numberofsamples , 0.8*PriorNonImplausibleSetRegularyIreRegular[, 8] , PriorNonImplausibleSetRegularyIreRegular[, 8] )
    PriorNonImplausibleSetRegularyIreRegular[, 9] <- runif(numberofsamples , 0.8*PriorNonImplausibleSetRegularyIreRegular[, 8] , PriorNonImplausibleSetRegularyIreRegular[, 8] )
    PriorNonImplausibleSetRegularyIreRegular[, 10] <- 1
    
    f_x_ReIre <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
    
    F_x_ReIre <- t(apply(f_x_ReIre , 1 , function(X){cumsum( c(0,diff(x))*(X/sum(c(0,diff(x))*(X))) ) }))
      
    PeakHeights <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvaluateDenistyEstimate(X[4:6], X) }))
    
    # max height must be central peak
    Valid1 = (PeakHeights[,1] < PeakHeights[,2])*(PeakHeights[,3] < PeakHeights[,2]) == 1
    
    # No daylight between peaks
    MinBetweenPeaks = apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){ min(c( min(FM_EvaluateDenistyEstimate(seq(X[4] , X[5] , 0.01), X) ) , min(FM_EvaluateDenistyEstimate(seq(X[5] , X[6] , 0.01), X) ))) } ) 
    Valid2 = MinBetweenPeaks > 0.2
    Valid = (Valid1*Valid2) == 1
    
    PriorNonImplausibleSetRegularyIreRegular <- PriorNonImplausibleSetRegularyIreRegular[Valid , ]
    f_x_ReIre <- f_x_ReIre[Valid , ]
    F_x_ReIre <- F_x_ReIre[Valid , ]

    MDFunction2 <- function(F_x , X){
      return(MDFunction(F_x , X))
    }
    
    MD_ReIre <- matrix(0 , dim(PriorNonImplausibleSetRegularyIreRegular)[1] , length(x))
    for( i in 1:dim(PriorNonImplausibleSetRegularyIreRegular)[1]){
      MD_ReIre[i , ] <- MDFunction2(F_x_ReIre[i , ] , PriorNonImplausibleSetRegularyIreRegular[i ,])
    }
  }
  
 # Calulate approximate correlation matrix.
  SampleMatrix  <-  FM_SampleRealisationsSet( PriorNonImplausibleSetRegularyIreRegular , N = (500 - 19) , x  )
  ME_matrix     <-  SampleMatrix - F_x_ReIre
  Cov_sdhat1    <-  cov( ME_matrix )
  Corr_sdhat1   <-  cov2cor( Cov_sdhat1 )
  
  ImplausabilityMatrix <- FM_CalulateImForGroundTruth(x = x , F_x = F_x_ReIre , PriorNonImplausibleSet= PriorNonImplausibleSetRegularyIreRegular , MD = MD_ReIre , Corr_sdhat1 , N = (500 - 19)  )
  ImThreshold1 <- apply(ImplausabilityMatrix  , 2 , function(X){quantile(X , 0.98)}) 
  
  SampleMatrix <- FM_SampleRealisationsSet( PriorNonImplausibleSet =  PriorNonImplausibleSetRegular[1:10000,] , N = (500 - 19) , x =  x  )
  ME_matrix <- SampleMatrix - F_xreg[1:10000,]
  Cov_sdhat2 <- cov( ME_matrix )
  Corr_sdhat2 <- cov2cor(DP_AddNugget(Cov_sdhat2 , 0.0045) )
  
  ImplausabilityMatrix <- FM_CalulateImForGroundTruth(x = x , F_x = F_xreg[1:10000,] , PriorNonImplausibleSet= PriorNonImplausibleSetRegular[1:10000,] , MD = MD_Reg[1:10000,], Corr_sdhat2 , N = (500 - 19)  )
  ImThreshold2 <- apply(ImplausabilityMatrix  , 2 , function(X){quantile(X , 0.98)}) 
  
}

Goodbad <- matrix(T , dim(PriorNonImplausibleSetRegularyIreRegular)[1])

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

BC_PlotCompareTwoHists(PriorNonImplausibleSetRegularyIreRegular[1:1000 , 1:9] , PriorNonImplausibleSetRegular[1:1000 , 1:9])
BC_PlotPairsFromTwoVariables(PriorNonImplausibleSetRegularyIreRegular[1:1000 , 1:9] , PriorNonImplausibleSetRegular[1:1000 , 1:9] , alpha = 0.1)

PatientID <- DP_choosepatient(listAllPatients)

MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = PatientID)
ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)

{
rangeofbeats <- 1:500
#rangeofbeats <- 30501:(30501 + 500)
#rangeofbeats <- 22000:22500
  
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

x11(20 , 14)
grid.arrange(RPeakdisPlot , PWavePlot1 ,PWavePlot2 , PWavePlot3 , layout_matrix = lay)
}

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

ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , F_x = F_xreg , MD = MD_Reg , RRtimes = RRtimes , Corr_sdhat = Corr_sdhat2, imthreshold = ImThreshold2[2] )
ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular , x = x , F_x =  F_x_ReIre , MD = MD_ReIre , RRtimes = RRtimes , Corr_sdhat = Corr_sdhat1 , imthreshold = ImThreshold1[2] )

if(length(ReHmOutput$Implausability) > 3){
  
samplenum <- which.min(ReHmOutput$Implausability[,3])
plot(  x , ReHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
lines( x , ReHmOutput$y , type ='l' , col = 'black' )
lines( x , ReHmOutput$f_x[samplenum , ] + 2*(MDFunction(ReHmOutput$f_x[ samplenum , ] ,ReHmOutput$NonImplausibleSets[samplenum , ] )+ 0.0045)  , col ='red' , lty=2 )
lines( x , ReHmOutput$f_x[samplenum , ] - 2*(MDFunction(ReHmOutput$f_x[ samplenum , ]  ,ReHmOutput$NonImplausibleSets[samplenum , ]  )+ 0.0045)  , col ='red' , lty=2 )

}else{
  
plot(  x , ReHmOutput$f_x , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
lines( x , ReHmOutput$y , type ='l' , col = 'black' )
lines( x , ReHmOutput$f_x + 2*(MDFunction(ReHmOutput$f_x ,ReHmOutput$NonImplausibleSets ) + 0.0045)  , col ='red' , lty=2 )
lines( x , ReHmOutput$f_x - 2*(MDFunction(ReHmOutput$f_x  ,ReHmOutput$NonImplausibleSets  )+ 0.0045)  , col ='red' , lty=2 )

}

if( length( ReIreHmOutput$Implausability ) > 3 ){
  
  samplenum <- which.min(ReIreHmOutput$Implausability[,1])
  plot(  x , ReIreHmOutput$f_x[samplenum, ] , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReIreHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReIreHmOutput$f_x[samplenum , ] + 3.82*(MDFunction(ReIreHmOutput$f_x[ samplenum , ] , ReIreHmOutput$NonImplausibleSets[samplenum , ] ) + 0.0045)  , col ='red' , lty=2 )
  lines( x , ReIreHmOutput$f_x[samplenum , ] - 3.82*(MDFunction(ReIreHmOutput$f_x[ samplenum , ]  ,ReIreHmOutput$NonImplausibleSets[samplenum , ]  )+ 0.0045)  , col ='red' , lty=2 )
  
  MDTmp <- MDFunction(ReIreHmOutput$f_x[ samplenum , ]  ,ReIreHmOutput$NonImplausibleSets[samplenum , ]  )+ 0.0045 
  
  covMattmp <- DP_AddNugget( ((MDTmp) %*% t(MDTmp)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[ii , NonZeroLogical]+ 0.00000001) )
  mahalanobis(x = y[NonZeroLogical] , center = F_x[ii , NonZeroLogical] , cov = covMattmp  ) 
  
}else{
  
  plot(  x , ReIreHmOutput$f_x , col = 'purple', type = 'l' , xlab = 'RRTimes' , ylab = 'Density'  )
  lines( x , ReIreHmOutput$y , type ='l' , col = 'black' )
  lines( x , ReIreHmOutput$f_x + 3.82*(MDFunction2(ReIreHmOutput$f_x ,ReHmOutput$NonImplausibleSets ))  , col ='red' , lty=2 )
  lines( x , ReIreHmOutput$f_x - 3.82*(MDFunction2(ReIreHmOutput$f_x  ,ReHmOutput$NonImplausibleSets  ))  , col ='red' , lty=2 )

  
  MDTmp <- MDFunction(ReIreHmOutput$f_x  ,ReIreHmOutput$NonImplausibleSets  )+ 0.0045 
  
  covMattmp <- DP_AddNugget( ((MDTmp) %*% t(MDTmp)) * Corr_sdhat , 0.0001*diag(MD[ii , NonZeroLogical]+ 0.00000001) )
  mahalanobis(x = y[NonZeroLogical] , center = F_x[ii , NonZeroLogical] , cov = covMattmp  ) 
  
  
}

BC_PlotPairsFromTwoVariables( PriorNonImplausibleSetRegularyIreRegular[1:5000,1:9],PriorNonImplausibleSetRegular[1:1000,1:9] , alpha = 0.02)
BC_PlotPairsFromTwoVariables( PriorNonImplausibleSetRegularyIreRegular[1:1000,1:9],ReIreHmOutput$NonImplausibleSets[,1:9] , alpha = 0.03)



###### Test with known ground truth ######

RRtimes <- FM_SampleGMM(X = PriorNonImplausibleSetRegularyIreRegular[1 , ] ,  N = 231)

ReHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegular, x = x , F_x = F_xreg , MD = MD_Reg , RRtimes = RRtimes , imthreshold = 3 )
ReIreHmOutput <- FM_HistoryMatchRRCulmativeDensity( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular,x = x ,F_x =  F_x_ReIre, MD = MD_ReIre, RRtimes = RRtimes, imthreshold = 3 )

if(length(ReHmOutput$Implausability) > 2){
  samplenum <- which.min(ReHmOutput$Implausability[,1])
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

if( length( ReIreHmOutput$Implausability[,1] ) > 2 ){
  
  samplenum <- which.min(ReIreHmOutput$Implausability[,1])
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





