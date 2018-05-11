source('ASWF_ChoosePatient.R')


# Plot discrete data
plot(1)
dev.off()
df<- data.frame(x<- DataSet$Data$tt , y <- DataSet$Data$HeartRate )
myplot <- ggplot( df , aes(x,y))  + geom_point(colour="blue", alpha=0.009) + 
  ggtitle(DataSet$MetaData$PseudoId) +
  xlab("Time") + ylab("Heart Rate") +
  geom_hline( yintercept = 130 , linetype="dashed" , color = "red" ) + 
  geom_hline( yintercept = 100 , linetype="dashed" , color = "black" ) +
  geom_hline( yintercept = 60  , linetype="dashed" , color = "blue" )  +
  geom_vline( xintercept = as.numeric(as.POSIXct(DataSet$MetaData$FirstNewAF)) , linetype="dashed" , color = "black" ) 
x11()
print(myplot)


interestingtimepoint <- which( as.vector(as.character(round.POSIXt(DataSet$Data$tt , units = c('hours')))) 
  == select.list(as.vector(as.character(unique(round.POSIXt(DataSet$Data$tt , units = c('hours')))))
   , preselect = DataSet$MetaData$FirstNewAF
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

WaveData <- WaveData[ max( 1 , timeindex - (numberhoursbefore*((60^2)/0.005)) ) : min(length(WaveData[ , 1]) , timeindex + (numberhoursafter*((60^2)/0.005)) ) , ]
WaveData <- ReturnWaveformwithPositiveOrientation(WaveData)

DataSet$Data <-  DataSet$Data[ ((DataSet$Data$tt > WaveData$Date[1])*(DataSet$Data$tt < WaveData$Date[length( WaveData$Date)]) == 1), ]

print('Extracing RA and R-R times')
RWaveExtractedData <- RPeakExtractionWavelet( WaveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8)
print('RA and R-R times Extracted')

timelist <- as.vector(as.character(round.POSIXt(WaveData[, 1] , units = 'mins')))

startindex = which(timelist == (select.list(unique(timelist)
                                            , preselect = NULL
                                            , multiple = FALSE
                                            , title = 'Select time to view ECGI'
                                            , graphics = TRUE ) ))
startindex <- startindex[1]

timeintervaloptions <- as.character(c(1:100))
timeinterval <- as.numeric(select.list(unique(timeintervaloptions)
                                       , preselect = '25'
                                       , multiple = FALSE
                                       , title = 'Select number of seconds of full ECG to be viewed'
                                       , graphics = TRUE ))

endindex <- startindex + (round(timeinterval)/0.005)
regionofinterest <- startindex:endindex


p1 <- ggplot(RWaveExtractedData , aes(t , RA)) + 
  geom_point(colour="blue", alpha=0.01) +
  ggtitle('R-amplitdes') +
  xlab("t") +
  ylab("RA") + coord_cartesian(ylim = c(50, 200)) 

p2 <- ggplot(RWaveExtractedData , aes(t , RR)) +
  geom_point(colour="blue", alpha=0.01) +
  ggtitle('R-R times') +
  xlab("t") +
  ylab("RR") + coord_cartesian(ylim = c(0.4, 1.2)) 

p4 <- ggplot(DataSet$Data , aes(tt , HeartRate)) +
  geom_point(colour="blue", alpha=0.1) +
  ggtitle('Discrete Heart Rate') +
  xlab("t") +
  ylab("HR") 
