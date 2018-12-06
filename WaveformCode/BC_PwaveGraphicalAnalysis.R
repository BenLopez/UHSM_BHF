

AFBeats <- RPeaksStruct$RRCombined[BC_CreateAnnotationFromInference(t = RPeaksStruct$RRCombined$t , AFLocations = AFLocations[ incidencedetected , ] )== 1,]
NAFBeats <- RPeaksStruct$RRCombined[AnnotatedAFInference== 0 , ]

AFNAFPlot1 <- BC_PlotPWaveAnalysis(ECG = ECGs$ECGI , Beats =  AFBeats , Beats2 =  NAFBeats, QSwidth  = 10 ) +
             ggtitle('ECGI')
AFNAFPlot2 <- BC_PlotPWaveAnalysis(ECG = ECGs$ECGII , Beats =  AFBeats , Beats2 =  NAFBeats, QSwidth  = 10 ) +
  ggtitle('ECGII')
AFNAFPlot3 <- BC_PlotPWaveAnalysis(ECG = ECGs$ECGIII , Beats =  AFBeats , Beats2 =  NAFBeats, QSwidth  = 10 ) +
  ggtitle('ECGIII')


x11(15,12)
print(grid.arrange( AFNAFPlot1,
                    AFNAFPlot2,
                    AFNAFPlot3,
                    nrow = 3 ,
                    ncol = 1  , 
                    top = paste0('P-Wave Analysis. Red: suspected AFib. Blue: not AFib ') ))

