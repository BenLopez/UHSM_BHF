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

load( 'ForwardModelling.RData')
source("LibrariesAndSettings.R" , print.eval  = TRUE )

UseAnnotatedData <- 1
precomputedfolderpath <- DP_SelectPrecomputedFolder()
source( 'FM_CreateRhythumPriors.R' )
source( 'CTEm_LoadDataandCreateEmulatorStructures.R'  )
source( 'CTEm_LoadDataandCreateEmulatorStructuresCDF.R' )
source( 'FM_CreateMeanPWavePriors.R' )

save.image(file = 'ForwardModelling.RData')

{
PatientID <- DP_choosepatient( listAllPatients )
source('FM_DecisionSupportSinglePatient.R')
}

{
  RegularLogicalTotal <- (RegularLogical ==1)&(RegularLogical2 ==1)&(RegularyIrregularLogical == 0)&(RegularyIrregularLogical2 == 0)
  IrregularlyIrregularLogical <- (RegularyIrregularLogical == 1)&(RegularyIrregularLogical2 == 1)&(RegularLogical ==0) & (RegularLogical2 ==0)
  BadDataLogical <- RpeaksFail ==1 | ECGAbscence ==1
  
  RegularlyIrregularLogical <- (RegularLogicalTotal ==0) &(IrregularlyIrregularLogical ==0)
  
  RegIrRegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = IrregularlyIrregularLogical , minutethreshold = 1)
  RegulartimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularLogicalTotal , minutethreshold = 1)
  IrtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = RegularlyIrregularLogical, minutethreshold = 1)
  BadDataTimePoints <-  ASWF_GetStartEndAF(timemat , logicaltimeseries = BadDataLogical, minutethreshold = 1)
  #UndecidedtimePoints <- ASWF_GetStartEndAF(timemat , logicaltimeseries = Undecided, minutethreshold = 1)
  
  
  t <- RPeakData$RRCombined$t[seq(1 , length(RPeakData$RRCombined$t) , 3)]   
     RR <- RPeakData$RRCombined$RR[seq(1 ,  length(RPeakData$RRCombined$t)  , 3)] 
     
       plotstruct <- ggplot(data.frame( t = t , RR=RR ), aes(t , RR) ) +
           geom_point( colour = 'blue' ,  alpha=0.03 ) +
           ggtitle('RR times') +
           xlab( "t" ) +
           ylab( "RR" ) + 
           coord_cartesian(ylim = c(0, 1.5))
       RRPlot <- BC_PlotAddAFLines(plotstruct = plotstruct , MetaData = MetaData)  
    
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegulartimePoints , fillcolor = 'green')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = IrtimePoints , fillcolor = 'orange')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = RegIrRegulartimePoints , fillcolor = 'red')
  RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = BadDataTimePoints , fillcolor = 'black')
  x11(40 , 25)
  RRPlot + ggtitle('PWaves')
  
  #pdf(file = paste0("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\Annotation" , PatientID , '.pdf') )
  #print(RRPlot)
  #dev.off()
  
  #save(outputstruct , file = paste0("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\HMOutput" , PatientID , '.RData'))
}

{

{
StartBeat <-  5501
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
  
  
SecondWaveLogicaltest <- FM_DubiousCaseLogic(RRtimes = RRtimes[!is.na(RRtimes)],
                                              ReIreHmOutput = ReIreHmOutput,
                                              ReHmOutput = ReHmOutput,
                                              ObservedIm_xx = ObservedIm_xx,
                                              EmulatorParametersCDFMean = EmulatorParametersCDFMean
                                              , EmulatorParametersCDFMax = EmulatorParametersCDFMax ,
                                                  RegularLogical)  
  
x11()
plot(kdmodel$eval.points , kdmodel$estimate , type ='l')

if(length(ReHmOutput$Implausability) > 2){
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


if( (mean((c(ReMeanProb[maxRegIndex],ReMaxProb[maxRegIndex])/c(IrReMeanProb[maxIrIrIndex],IrReMaxProb[maxIrIrIndex]))^2)>2 ||  sqrt(sum((c(ReMeanProb[maxRegIndex],ReMaxProb[maxRegIndex]) - c(IrReMeanProb[maxIrIrIndex],IrReMaxProb[maxIrIrIndex]))^2)) > 0.1) & sum(RegularLogical)>1 ){
  print(TRUE)
}else{
  print(FALSE)
}


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



