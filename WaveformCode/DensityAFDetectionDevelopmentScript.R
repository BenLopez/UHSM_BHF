pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))

source("LibrariesAndSettings.R" , print.eval  = TRUE )

load(choose.files(caption = "Select patientindexmaster.RData"))
DP_ChooseDataReps()

DP_choosepatient(listAllPatients)
sub_pat = subset( PatIndex2017, PseudoId %in% subList )

print('Loading ECG.')
WaveData <- DP_LoadECG(path , subList , numberrep , ECGNum = 1)
print('ECG Loaded.')

interestingtimepoint <- DP_SelectInterestingTimePoint(WaveData[ seq(from = 1 , to = length(WaveData[ , 1]), by = 1000) , ] , sub_pat)
interestingindex <- which.min( abs(WaveData[ , 1] - as.POSIXct(interestingtimepoint)))
WaveData <- DP_CropWaveData(WaveData , interestingindex , DP_SelectHoursBeforeandAfter() )

print('Calulating Rpeaks')
RWaveExtractedData <- CleanRpeaks(RPeakExtractionWavelet( WaveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8) , 2)
print('Rpeaks Calulated')


AFScore <- ExtractIHVAFScore(RWaveExtractedData ,  binlims <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3))
VBM <- AFScore$IHAVFScore

p <- ggplot(RWaveExtractedData , aes(t , RR)) +
  geom_point(colour="blue", alpha=0.01) +
  xlab("t") +
  ylab("RR") + coord_cartesian(ylim = c(0, 1.2)) +
  geom_vline( xintercept = as.numeric( as.POSIXct( DP_StripTime(sub_pat$ConfirmedFirstNewAF)) ) , linetype="dashed" , color = "purple" )+
  ggtitle( paste0(subList , " RR-times" ) ) 
p2 <- ggplot(data.frame(x=AFScore$t , y=VBM), aes(x,y)) + geom_line(colour = "blue") + ggtitle('AF Score') + ylim(c(0,600))

tmp <- ASWF_GetStartEndAF(RWaveExtractedData$t , logicaltimeseries = (VBM > 150) )

for( i in ( 1:length(tmp$Start) ) )
{
p <- p + annotate("rect" , xmin = tmp$Start[i], xmax = tmp$End[i], ymin = -1000, ymax= 1000 , fill = 'pink' , alpha = 0.5)
}

print(grid.arrange( p  + 
                      geom_vline( xintercept = as.numeric( as.POSIXct(sub_pat$FirstNewAF[1] ) ) ,linetype="dashed" , color = "black" ) ,
                    p2 + geom_hline( yintercept = 150 ,linetype="dashed" , color = "black" ) , nrow = 2 , ncol = 1))


binmatrix <- CalulateBinMatrix(RWaveExtractedData ) 

par(mfow = c(2 , 1))
plot(RWaveExtractedData$t , RWaveExtractedData$RR)
abline(v = as.numeric( as.POSIXct(sub_pat$FirstNewAF[1] )))
plot(binmatrix[10000 , ])


