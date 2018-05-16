pathFiles <- choose.dir(caption="Select folder with source code.")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

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

p <- ggplot(RWaveExtractedData , aes(t , RR)) +
  geom_point(colour="blue", alpha=0.01) +
  xlab("t") +
  ylab("RR") + coord_cartesian(ylim = c(0, 1.2)) +
  geom_vline( xintercept = as.numeric( as.POSIXct( sub_pat$FirstNewAF) ) , linetype="dashed" , color = "purple" )+
  ggtitle( paste0(" RR-times" ) )
p

timeintervalforhist <- 5 # time interval in minutes

starttime <- RWaveExtractedData[1 , 1]


timelogical <- ( (RWaveExtractedData[ , 1] > starttime)*(RWaveExtractedData[ , 1] < (starttime + (timeintervalforhist*60)) )*(RWaveExtractedData$RR < 2) ) == 1
p2 <- ggplot(RWaveExtractedData[ timelogical , ] , aes( RR)) +
    geom_histogram(colour="blue")

#dev.off()
grid.arrange( p  + 
                geom_vline( xintercept = as.numeric( as.POSIXct(starttime ) ) ,linetype="dashed" , color = "black" ) + 
                geom_vline( xintercept = as.numeric( as.POSIXct(starttime + (timeintervalforhist*60) )) ,linetype="dashed" , color = "black" ) ,
               p2 , nrow = 2 , ncol = 1)

starttime <- starttime + (timeintervalforhist*60)


Beta <- cumsum(as.numeric(RWaveExtractedData$RR))/cumsum(rep(1,length(as.numeric(RWaveExtractedData$RR))))
W_t <- smth(RWaveExtractedData$RR - Beta , method = 'sma' , n = 500 )

tmp <- as.numeric(RWaveExtractedData$RR - Beta - W_t)
tmp <- tmp[!is.na(tmp)]
#tmp <- tmp[tmp < quantile(tmp , probs = 0.999 , na.rm = 'TRUE') ]
#tmp <- tmp[tmp > quantile(tmp , probs = 0.001 , na.rm = 'TRUE') ]

V_t <- smth(as.numeric(tmp)^2 , method = 'sma' , n = 500 ) - smth((tmp) , method = 'sma' , n = 500 )^2
Beta_V <- cumsum(V_t[!is.na(V_t)])/cumsum(rep(1,length(as.numeric(V_t[!is.na(V_t)]))))

#Sigma_W <- cumsum(V_t[!is.na(V_t)] - Beta_V)/
W_V_t <- smth(V_t[!is.na(V_t)] - Beta_V , method = 'sma' , n = 100 )
plot(W_V_t)

