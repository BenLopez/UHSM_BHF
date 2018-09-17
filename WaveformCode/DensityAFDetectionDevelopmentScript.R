{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )
DP_LoadPatientIndex()
DP_ChooseDataReps()}


DP_choosepatient(listAllPatients)
sub_pat = subset( PatIndex2017, PseudoId %in% subList )

print('Loading ECG.')
WaveData <- DP_LoadECGReduced(path , subList , numberrep , ECGNum = 1)
print('ECG Loaded.')

print('Calulating Rpeaks')
outputdata <- DP_LoadRpeaksfile(path , subList) 
print('Rpeaks Calulated')

AFScore <- ExtractIHVAFScore(outputdata$RRCombined ,  binlims <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3))
VBM <- AFScore$IHAVFScore

p <- ggplot(outputdata$RRCombined , aes(t , RR)) +
  geom_point(colour="blue", alpha=0.01) +
  xlab("t") +
  ylab("RR") + coord_cartesian(ylim = c(0, 1.2)) +
  ggtitle( paste0(subList , " RR-times" ) ) 


binMatrix <- AFD_CalulateBinMatrixKernelDensityEstimated(outputdata$RRCombined , n = 100)

SettingsAFDetection <- AFD_CreateDefaultSettings()
#output <-  ExtractNumberofModes( outputdata$RRCombined  , densitythresh = 0.025 )
NumberModes <- AFD_Calculatemodalmode( outputdata$RRCombined , 
                        binlims = SettingsAFDetection[['BinlimsMM']] , 
                        n = SettingsAFDetection[['BandWidthScore']] ,
                        densitythresh = SettingsAFDetection[['DensityThresholdMM']],
                        nn = SettingsAFDetection[['BadnWidthMM']])

AFScore<- AFD_ExtractIHVAFScore(outputdata$RRCombined ,  binlims = SettingsAFDetection[['BinlimsScore']] , n = SettingsAFDetection[['BandWidthScore']] )

x11()
indexOI <- 23500
p2 <- ggplot(data.frame(x = (c(1:dim(binMatrix)[2])*2.75)/100  , y = binMatrix[indexOI , ]) , aes(x , y)) +
  geom_line() + ylab('Density') + xlab('RR Times')+ ggtitle('Black')
indexOI <- 30000
p3 <- ggplot(data.frame( x = (c(1:dim(binMatrix)[2])*2.75)/100   , y = binMatrix[indexOI , ]) , aes(x , y)) +
  geom_line() + ylab('Density') + xlab('RR Times') + ggtitle('Green')
p4 <- ggplot(data.frame( x = (c(1:dim(binMatrix)[2])*2.75)/100   , y = binMatrix[5000 , ]) , aes(x , y)) +
  geom_line() + ylab('Density') + xlab('RR Times') +ggtitle('Blue')
grid.arrange( p + geom_vline(xintercept = as.numeric( NumberModes$t[indexOI]) , col = 'green') + geom_vline(xintercept = as.numeric( NumberModes$t[5000]) , col = 'blue')  , p3 , p2, p4, nrow = 2 , ncol =2 )



binlimits = (c(1:dim(binMatrix)[2])/73)*1.8

plot(binlimits , binMatrix[indexOI , ] , type ='l')
d1 <- c(0,diff(binMatrix[indexOI , ]))
d2 <- c(0,0,diff(diff(binMatrix[indexOI , ])))
peakslogical<- ((d1<0.025)*(d2<0)) ==1
peakslogical[is.na(peakslogical)]=FALSE
peakslogical <- FindLocalTurningPoints(peakslogical , binMatrix[indexOI , ])
points( binlimits[peakslogical] , binMatrix[indexOI , peakslogical])
