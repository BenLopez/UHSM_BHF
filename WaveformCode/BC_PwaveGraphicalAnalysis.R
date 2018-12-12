
if(!exists('Xstar')){
Xstar = seq(0.5 ,1 , 0.01)}
if(!exists('PriorNonImplausibleSet')){
PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95, 0.001 ) ,b = c(0.6  , 0.1 ) , numbersamples = 10000)
}

AFBeats <- RPeaksStruct$RRCombined[BC_CreateAnnotationFromInference(t = RPeaksStruct$RRCombined$t , AFLocations = AFLocations[ incidencedetected , ] )== 1,]
regiontocheck <- 500:min(1000 , dim(AFBeats)[1])
lengthcheck <- 0
while(lengthcheck < 400){
  ECGBeats <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = AFBeats[regiontocheck,] , QSwidth  = 10)
  lengthcheck <- length(ECGBeats$numvalues)
  if(lengthcheck < 400 ){regiontocheck <- regiontocheck + 500}
}

PAmplitudeAF <- PWaveHM_EmulateEstimatePAmplitude(QS_Struct = ECGBeats , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar = Xstar ,PriorNonImplausibleSet =  PriorNonImplausibleSet)

NAFBeats <- RPeaksStruct$RRCombined[AnnotatedAFInference== 0 , ]

regiontocheck <- 1000:min(1500 , dim(NAFBeats)[1])
lengthcheck <- 0
while(lengthcheck < 400){
ECGBeats <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = NAFBeats[regiontocheck,] , QSwidth  = 10)
lengthcheck <- length(ECGBeats$numvalues)
if(lengthcheck < 400 ){regiontocheck <- regiontocheck + 500}
}

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
                    top = paste0('P-Wave Analysis. Red: suspected AFib. Blue: not AFib. ECGII Amplitude % Difference ' , round(( (abs(PAmplitudeNAF - PAmplitudeAF)) / abs(PAmplitudeNAF) )*100 ) , '%' )))
UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to P-Wave Inference Diagnostics?')
if(UserResponse == 'YES'){
PWaveHM_EmulateEstimatePAmplitude(QS_Struct = ECGBeats , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar = Xstar ,PriorNonImplausibleSet =  PriorNonImplausibleSet , Graphics = 1)
ECGBeats <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = AFBeats[500:min(1000 , dim(AFBeats)[1]),] , QSwidth  = 10)
PWaveHM_EmulateEstimatePAmplitude(QS_Struct = ECGBeats , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar , PriorNonImplausibleSet, Graphics = 1)
}
