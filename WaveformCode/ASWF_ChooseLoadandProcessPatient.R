source('ASWF_ChoosePatient.R')

# Plot discrete data
plot(1)
dev.off()
myplot <- ggplot( data.frame(x<- DataSet$Data$tt , y <- DataSet$Data$HeartRate ) , aes(x,y))  + geom_point(colour="blue", alpha=0.009) + 
  ggtitle(DataSet$MetaData$PseudoId) +
  xlab("Time") + ylab("Heart Rate") +
  geom_hline( yintercept = 130 , linetype="dashed" , color = "red" ) + 
  geom_hline( yintercept = 100 , linetype="dashed" , color = "black" ) +
  geom_hline( yintercept = 60  , linetype="dashed" , color = "blue" )  +
  geom_vline( xintercept = as.numeric(as.POSIXct(DataSet$MetaData$FirstNewAF[1])) , linetype="dashed" , color = "black" ) +
  geom_vline( xintercept = as.numeric(as.POSIXct(DataSet$MetaData$LastITUEntry[1])) , linetype="dashed" , color = "red" ) 
x11(12 , 7)
print(myplot)

interestingtimepoint <- which( as.vector(as.character(round.POSIXt(DataSet$Data$tt , units = c('hours')))) 
  == select.list(as.vector(as.character(unique(round.POSIXt(DataSet$Data$tt , units = c('hours')))))
   , preselect = DataSet$MetaData$FirstNewAF[1]
   , multiple = FALSE
   , title = 'Select Interest Time Point'
   , graphics = TRUE ))


# Load wave form data 
for(i in 1:(numberrep+1))
{
  if(i > (numberrep))
  {
    stop('Error: No ECGI data processed.') 
    break
  }
  if(file.exists(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGI_' , subList , '.RData' )))
  {
    load(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGI_' , subList , '.RData' ))
    break 
  }  
}

timeindex <- which.min( abs(difftime( WaveData$Date ,  DataSet$Data$tt[interestingtimepoint] , units = 'secs')) )

numberhours <- c('1' , '2' , '3' , '4' , '5' , '6' , '7' , '8')

numberhoursbefore = as.numeric(select.list(numberhours
                                           , preselect = '5'
                                           , multiple = FALSE
                                           ,title = 'Select number of hours before.'
                                           , graphics = TRUE ))

numberhoursafter = as.numeric(select.list(numberhours
                                          , preselect = '1'
                                          , multiple = FALSE
                                          ,title = 'Select number of hours after.'
                                          , graphics = TRUE ))

print('Loading ECGI.')
ECGI <- WaveData[ max( 1 , timeindex - (numberhoursbefore*((60^2)/0.005)) ) : min(length(WaveData[ , 1]) , timeindex + (numberhoursafter*((60^2)/0.005)) ) , ]
ECGI <- ReturnWaveformwithPositiveOrientation(ECGI)
print('ECGI Loaded.')

# Load wave form data 
for(i in 1:(numberrep+1))
{
  if(i > (numberrep))
  {
    stop('Error: No ECGII data processed.') 
    break
  }
  if(file.exists(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGII_' , subList , '.RData' )))
  {
    load(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGII_' , subList , '.RData' ))
    break 
  }  
}

timeindex <- which.min( abs(difftime( WaveData$Date ,  DataSet$Data$tt[interestingtimepoint] , units = 'secs')) )

print('Load ECGII.')
ECGII <- WaveData[ max( 1 , timeindex - (numberhoursbefore*((60^2)/0.005)) ) : min(length(WaveData[ , 1]) , timeindex + (numberhoursafter*((60^2)/0.005)) ) , ]
ECGII <- ReturnWaveformwithPositiveOrientation( ECGII )
print('ECGII loaded.')
rm( WaveData )

DataSet$Data <-  DataSet$Data[ ((DataSet$Data$tt > ECGI$Date[1])*(DataSet$Data$tt < ECGI$Date[length( ECGI$Date)]) == 1), ]

print('Extracting RA and R-R times.')
RWaveExtractedDataI <- CleanRpeaks(RPeakExtractionWavelet( ECGI , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8))
#RWaveExtractedDataII <- RPeakExtractionWavelet( WaveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8)
print('RA and R-R times Extracted.')

AFScore <- ExtractIHVAFScore(RWaveExtractedDataI ,  binlims <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3))
StartEndTimesAF <- ASWF_GetStartEndAF(t = AFScore$t , logicaltimeseries = (AFScore$IHAVFScore > 150)  , minutethreshold = 9)

timelist <- as.vector(as.character(round.POSIXt(ECGI[seq(from = 1 , to = length(ECGI[ , 1]) , by = 1000), 1] , units = 'mins')))

startindex = which(timelist == (select.list(unique(timelist)
                                            , preselect = NULL
                                            , multiple = FALSE
                                            , title = 'Select time to view ECGI'
                                            , graphics = TRUE ) ))
startindex <- startindex[1]

timeintervaloptions <- as.character(c(1:100))
timeinterval <- as.numeric(select.list(unique(timeintervaloptions)
                                       , preselect = '10'
                                       , multiple = FALSE
                                       , title = 'Select number of seconds of full ECG to be viewed'
                                       , graphics = TRUE ))

endindex <- startindex + (round(timeinterval / as.numeric(abs(ECGI$Date[1]- ECGI$Date[2]))))
regionofinterest <- startindex:endindex

regionofinterest2  <- ASWF_AlignRegionofInterests(ECGI , ECGII , regionofinterest)

startindex <- which.min( abs(ECGII[ , 1] - ECGI[ regionofinterest[1] , 1]) )
endindex <- startindex + length(regionofinterest)
regionofinterest2 <- startindex:endindex

p1 <- ggplot(RWaveExtractedDataI , aes(t , RA)) + 
  geom_point(colour="blue", alpha=0.01) +
  ggtitle('R-amplitudes') +
  xlab("t") +
  ylab("RA") + coord_cartesian(ylim = c(50, 200)) 


p2 <- ggplot() + 
      geom_line( data = AFScore , aes(x = t , y = IHAVFScore/300) , colour ='red' , alpha = 0.25)  + 
      geom_point(data = RWaveExtractedDataI  , aes(x = t , y = RR) , colour="blue", alpha=0.01)+
      scale_y_continuous(sec.axis = sec_axis(~.*300, name = "AF Score")) +
      xlab("t") +
      ylab("RR") + coord_cartesian(ylim = c(0.2, 1.2))


if(length(StartEndTimesAF$Start) > 0)
{
for( i in ( 1:length(StartEndTimesAF$Start) ) )
{
  p2 <- p2 + annotate("rect" , xmin = StartEndTimesAF$Start[i], xmax = StartEndTimesAF$End[i], ymin = -1000, ymax= 1000 , fill = 'pink' , alpha = 0.25)
}
}

p4 <- ggplot(DataSet$Data , aes(tt , HeartRate)) +
  geom_point(colour="blue", alpha=0.1) +
  ggtitle('Discrete Heart Rate') +
  xlab("t") +
  ylab("HR") 
