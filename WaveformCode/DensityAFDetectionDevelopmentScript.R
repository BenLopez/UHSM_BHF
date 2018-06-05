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

NumberModes <- AF_Calculatemodalmode( RWaveExtractedData )

p3 <- ggplot( data.frame(x = NumberModes$t , y = NumberModes$NumModes) , aes(x , y) ) +
      geom_line( ) +
      ggtitle('Number of Modes')

grid.arrange( p  , p3  , nrow = 2 , ncol =1 )




