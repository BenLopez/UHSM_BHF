

PAmplitudes <- BC_CalculatePAmplitudes(AFLocations = AFLocations[incidencedetected,] ,RPeaksStruct =  RPeaksStruct ,ECGs =  ECGs , AnnotatedAFInference = AnnotatedAFInference)


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
                    top = paste0('P-Wave Analysis. Red: suspected AFib. Blue: not AFib. ECGII Amplitude % Difference ' , round(( (abs(PAmplitudes$NAFPAmplitude - PAmplitudes$AFPAmplitude)) / abs(PAmplitudes$NAFPAmplitude) )*100 ) , '%' )))
UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to P-Wave Inference Diagnostics?')
if(UserResponse == 'YES'){
  PAmplitudes <- BC_CalculatePAmplitudes(AFLocations = AFLocations[incidencedetected,] ,RPeaksStruct =  RPeaksStruct ,ECGs =  ECGs , AnnotatedAFInference = AnnotatedAFInference , Graphics = 1)
  }
