pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))

source("LibrariesAndSettings.R" , print.eval  = TRUE )

DP_LoadPatientIndex()
DP_ChooseDataReps()

DP_choosepatient(listAllPatients)
sub_pat = subset( PatIndex2017, PseudoId %in% subList )

print('Loading ECG.')
WaveData <- DP_LoadECGReduced(path , subList , numberrep , ECGNum = 1)
print('ECG Loaded.')

print('Calulating Rpeaks')
RWaveExtractedData <- CleanRpeaks(RPeakExtractionWavelet( WaveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8) , 2)
print('Rpeaks Calulated')

AFScore <- ExtractIHVAFScore(RWaveExtractedData ,  binlims <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3))
VBM <- AFScore$IHAVFScore

p <- ggplot(RWaveExtractedData , aes(t , RR)) +
  geom_point(colour="blue", alpha=0.01) +
  xlab("t") +
  ylab("RR") + coord_cartesian(ylim = c(0, 1.2)) +
  ggtitle( paste0(subList , " RR-times" ) ) 


binMatrix <- AFD_CalulateBinMatrixKernelDensityEstimated(RWaveExtractedData , n = 100)

#output <-  ExtractNumberofModes( RWaveExtractedData  , densitythresh = 0.025 )
NumberModes <- AFD_Calculatemodalmode( RWaveExtractedData , binlims=  c(0, seq(from = 0  , to = 1.8  , 0.025  )) , n =100 , nn = 250  , densitythresh = 0.01 )


NumberModes <- AFD_Calculatemodalmode( RWaveExtractedData , 
                        binlims = SettingsAFDetection[['BinlimsMM']] , 
                        n = SettingsAFDetection[['BandWidthScore']] ,
                        densitythresh = SettingsAFDetection[['DensityThresholdMM']],
                        nn = SettingsAFDetection[['BadnWidthMM']])

AFScore<- AFD_ExtractIHVAFScore(RWaveExtractedData ,  binlims = SettingsAFDetection[['BinlimsScore']] , n = SettingsAFDetection[['BandWidthScore']] )

indexOI <- 25000
p2 <- ggplot(data.frame(x = c(1:dim(binMatrix)[2])  , y = binMatrix[indexOI , ]) , aes(x , y)) +
      geom_line() + geom_hline(yintercept = 0.01)
p3 <- ggplot( data.frame(x = NumberModes$t , y = NumberModes$NumModes) , aes(x , y) ) +
      geom_line( ) +
      ggtitle('Number of Modes')

grid.arrange( p + geom_vline(xintercept = as.numeric( NumberModes$t[indexOI]))  , p2, p3  , nrow = 3 , ncol =1 )




