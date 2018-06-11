pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )

DP_LoadPatientIndex()
DP_ChooseDataReps()
DP_choosepatient(listAllPatients)

ECGI <- DP_LoadECGReduced(path , subList   , numberrep , 1)
ECGII <- DP_LoadECGReduced(path , subList , numberrep , 2)
ECGIII <-  DP_LoadECG(path , subList , numberrep , 3)
sub_pat <- subset(PatIndex2017, PseudoId %in% subList)


timeindex <- which.min( abs(difftime( ECGIII$Date , ECGI$Date[1] , units = 'secs')) )
numberhoursbefore = 0
numberhoursafter = abs( difftime( ECGI$Date[length(ECGI$Date)] , ECGI$Date[1]  , units ='hours' ))
HoursBeforeAndAFter = data.frame(numberhoursbefore , numberhoursafter)
ECGIII <- ReturnWaveformwithPositiveOrientation(DP_CropWaveData(ECGIII , timeindex , HoursBeforeAndAFter))

outputdata<- list()
outputdata[[1]] <- CleanRpeaks(RPeakExtractionWavelet( ECGI   , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
outputdata[[2]] <- CleanRpeaks(RPeakExtractionWavelet( ECGII  , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
outputdata[[3]] <- CleanRpeaks(RPeakExtractionWavelet( ECGIII , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5) , 2)
outputdata[[4]] <- sub_pat

outputdata <- setNames(outputdata , c('ECGI' , 'ECGII' , 'ECGIII' , 'Meta_Data'))


regionofinterestI   <-  DP_ChooseRegionofInterest( ECGI )
regionofinterestII  <-  DP_AlignRegionofInterests( ECGI , ECGII , regionofinterestI )
regionofinterestIII <-  DP_AlignRegionofInterests( ECGI , ECGIII , regionofinterestI )

p1 <- ggplot(ECGI[regionofinterestI , ] , aes(Date , Value)) +
  geom_line(colour="blue") + ylab('Hz') +
  geom_point(data = outputdata$ECGI[ ((outputdata$ECGI$t > ECGI$Date[regionofinterestI[1]])*(outputdata$ECGI$t < ECGI$Date[regionofinterestI[length(regionofinterestI)]]) == 1) , ] , aes(t , RA) ) +
  xlim(ECGI[regionofinterestI[1] , 1] , ECGI[regionofinterestI[length(regionofinterestI)] , 1] )  +
  ggtitle('ECGI')

p2 <- ggplot(ECGII[regionofinterestII , ] , aes(Date , Value)) +
  geom_line(colour="red") + ylab('Hz') +
  geom_point(data = outputdata$ECGII[ ((outputdata$ECGII$t > ECGII$Date[regionofinterestII[1]])*(outputdata$ECGII$t < ECGII$Date[regionofinterestII[length(regionofinterestII)]]) == 1) , ] , aes(t , RA) ) +
  xlim(ECGI[regionofinterestI[1] , 1] , ECGI[regionofinterestI[length(regionofinterestI)] , 1] )  +
  ggtitle('ECGII')  

p3 <- ggplot(ECGIII[regionofinterestIII , ] , aes(Date , Value)) +
  geom_line(colour="green") + ylab('Hz') +
  geom_point(data = outputdata$ECGIII[ ((outputdata$ECGIII$t > ECGIII$Date[regionofinterestIII[1]])*(outputdata$ECGIII$t < ECGIII$Date[regionofinterestIII[length(regionofinterestIII)]]) == 1) , ] , aes(t , RA) ) +
  xlim(ECGI[regionofinterestI[1] , 1] , ECGI[regionofinterestI[length(regionofinterestI)] , 1] )  +
  ggtitle('ECGIII')


p <- ggplot( ) + 
  geom_point( data = outputdata$ECGI    , aes(x = t , y = RR) , colour="blue"  , alpha=0.01 ) +
  geom_point( data = outputdata$ECGII   , aes(x = t , y = RR) , colour="red"   , alpha=0.01 ) +
  geom_point( data = outputdata$ECGIII  , aes(x = t , y = RR) , colour="green" , alpha=0.01 ) +
  scale_y_continuous( sec.axis = sec_axis(~.*150, name = "AF Score") ) +
  xlab( "t" ) +
  ylab( "RR" ) + coord_cartesian(ylim = c(0, 1.2)) + 
  ggtitle(paste0(subList , ' RR-Times')) + 
  geom_vline( xintercept = as.numeric(ECGI[regionofinterestI[1] , 1]) , linetype="dashed" , color = "black" ) + 
  geom_vline( xintercept = as.numeric(ECGI[regionofinterestI[length(regionofinterestI)] , 1]) , linetype="dashed" , color = "black" )
x11(15,10)
grid.arrange(p1,p2,p3 , p , nrow = 4, ncol = 1)

output1 <- PE_MultipleECGRPeaks(outputdata , thresh = 0.0075)

p5 <- ggplot( ) + 
  geom_point( data = output1  , aes(x = t , y = RR) , colour="blue" , alpha=0.01 ) +
  scale_y_continuous( sec.axis = sec_axis(~.*150, name = "AF Score") ) +
  xlab( "t" ) +
  ylab( "RR" ) + coord_cartesian(ylim = c(0, 1.2)) + 
  ggtitle(paste0(subList , ' RR-Times Combined')) + 
  geom_vline( xintercept = as.numeric(ECGI[regionofinterestI[1] , 1]) , linetype="dashed" , color = "black" ) + 
  geom_vline( xintercept = as.numeric(ECGI[regionofinterestI[length(regionofinterestI)] , 1]) , linetype="dashed" , color = "black" )

x11(15,10)
grid.arrange(p1 , p2 , p3 , p , p5 , nrow = 5 , ncol = 1)
