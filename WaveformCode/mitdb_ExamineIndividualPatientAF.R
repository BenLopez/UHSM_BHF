files <- choose.files(multi = FALSE)


tmp <- readMat(files)
waveData <- setNames(  data.frame(tmp$Data[[2]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2)] ,
                                  tmp$Data[[1]][seq(1 , length(tmp$Data[[1]][ , 1]) , 2) , 1] ) , c('Date' , 'Value'))
waveData <-ReturnWaveformwithPositiveOrientation( waveData )
Annotation <- setNames(data.frame( round(tmp$Data[[3]]/2) ,  tmp$Data[[4]])  , c('Index' , 'Code'))
stdthresh <- 2.5
RData <- CleanRpeaks(RPeakExtractionWavelet( waveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = stdthresh) , 2) 
Resultsrpeaks <- mitdb_ComputePeakresults( waveData , Annotation , stdthresh = stdthresh )

AFScore <- ExtractIHVAFScore(RData ,  binlims <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3))
#StartEndTimesAF <- ASWF_GetStartEndAF(t = AFScore$t , logicaltimeseries = (AFScore$IHAVFScore > 145)  , minutethreshold = 9)

#p1 <- ggplot(RData , aes(t , RA)) + 
# geom_point(colour="blue", alpha=0.1) +
# ggtitle('R-amplitudes') +
#  xlab("t") +
#  ylab("RA") + coord_cartesian(ylim = c(0.5, 2)) 


p2 <- ggplot() + 
  geom_line( data = AFScore , aes(x = t , y = IHAVFScore/300) , colour ='red' , alpha = 0.25)  + 
  geom_point(data = RData  , aes(x = t , y = RR) , colour="blue", alpha=0.1)+
  scale_y_continuous(sec.axis = sec_axis(~.*300, name = "AF Score")) +
  xlab("t") +
  ylab("RR") + coord_cartesian(ylim = c(0.2, 1.2))
  
p2
