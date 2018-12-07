Xstar = seq(0.5 ,1 , 0.01)
PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 1, 0.001 ) ,b = c(0.6  , 0.1 ) , numbersamples = 10000)

AFBeats <- RPeaksStruct$RRCombined[BC_CreateAnnotationFromInference(t = RPeaksStruct$RRCombined$t , AFLocations = AFLocations[ incidencedetected , ] )== 1,]
ECGBeats <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = AFBeats[500:min(1000 , dim(AFBeats)[1]),] , QSwidth  = 10)
PAmplitudeAF <- PWaveHM_EmulateEstimatePAmplitude(QS_Struct = ECGBeats , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar = Xstar ,PriorNonImplausibleSet =  PriorNonImplausibleSet)

NAFBeats <- RPeaksStruct$RRCombined[AnnotatedAFInference== 0 , ]
ECGBeats <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = NAFBeats[1000:min(1500 , dim(NAFBeats)[1]),] , QSwidth  = 10)
PAmplitudeNAF <- PWaveHM_EmulateEstimatePAmplitude(QS_Struct = ECGBeats , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar , PriorNonImplausibleSet)



AFNAFPlot1 <- BC_PlotPWaveAnalysis(ECG = ECGs$ECGI , Beats =  AFBeats , Beats2 =  NAFBeats , QSwidth  = 10 ) +
             ggtitle('ECGI') + ylim(-40,40)

AFNAFPlot2 <- BC_PlotPWaveAnalysis(ECG = ECGs$ECGII , Beats =  AFBeats , Beats2 =  NAFBeats, QSwidth  = 10 ) +
           ggtitle('ECGII') + ylim(-40,40)


AFNAFPlot3 <- BC_PlotPWaveAnalysis(ECG = ECGs$ECGIII , Beats =  AFBeats , Beats2 =  NAFBeats, QSwidth  = 10 ) +
  ggtitle('ECGIII') + ylim(-40,40)



x11(15,12)
print(grid.arrange( AFNAFPlot1,
                    AFNAFPlot2,
                    AFNAFPlot3,
                    nrow = 3 ,
                    ncol = 1  , 
                    top = paste0('P-Wave Analysis. Red: suspected AFib. Blue: not AFib. ECGII Amplitude % Difference ' , round(( (abs(PAmplitudeNAF) - abs(PAmplitudeAF)) /PAmplitudeNAF )*100 ) , '%' )))

