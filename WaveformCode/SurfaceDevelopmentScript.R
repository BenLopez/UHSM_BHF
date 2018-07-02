{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )
DP_GetDirectories()}

PatientID <- DP_choosepatient( listAllPatients )

ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = DP_ChooseECGstoProcess()  )
PeakData <- DP_LoadRpeaksfile(path , PatientID)

QSwidth = 0
ylim <- c( -60 , 200 )
x11(20,20)
{dev.off()
  x11(20,20)
  par(mfrow = c(1 , 1))
QS <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = PeakData$RRCombined[30000:31000,] , QSwidth = QSwidth)
plot( 0  , 0  , xlim = c(-0,1) , ylim = ylim , xlab='t' , ylab ='Hz'  , type ='l')
title('T-P Waves Distribution in AF')
goodbeatslogical <- (QS$numvalues < as.numeric(quantile( QS$numvalues  , 0.95)))*(QS$numvalues > as.numeric(quantile( QS$numvalues  , (1-0.95)))) == 1

tmp <- QS$Value[ goodbeatslogical , ]
tmp<- apply(tmp , 1 , function(X){X - mean(X , na.rm = TRUE)})

m <- apply(tmp , 1 , function(X){mean(X , na.rm = TRUE)})
v <- apply(tmp , 1 , function(X){var(X , na.rm = TRUE)})
  
for(i in 1:dim(QS$Date)[1]){
  if(QS$numvalues[i ] > as.numeric(quantile( QS$numvalues  , 0.95)) ){next}
  if(QS$numvalues[i ] < as.numeric(quantile( QS$numvalues  , 0.5)) ){next}
  tmp = 1:(min(dim(QS$Date)[2] , QS$numvalues[i ]) -1)
  lines((1:length(QS$Date[i,tmp]))/length(QS$Date[i,tmp]), QS$Value[i,tmp] - mean(QS$Value[i,tmp]) , col= rgb(1 , 0 , 0 , alpha = 0.01))
}

#lines( (1:length(m))/length(m) , m , col = 'blue')
#lines( (1:length(m))/length(m), m + 2*sqrt(v) , col = 'black')
#lines( (1:length(m))/length(m), m - 2*sqrt(v) , col = 'black')

#plot( 0  , 0  , xlim = c(-0,1) , ylim = ylim , xlab='t' , ylab ='Hz'  , type ='l')

#title('T-P Second order Statistics')
#lines( (1:length(m))/length(m) , m , col = 'blue')
#lines( (1:length(m))/length(m), m + 2*sqrt(v) , col = 'black')
#lines( (1:length(m))/length(m), m - 2*sqrt(v) , col = 'black')

QS <- AFD_ExtractAllSQ(ECGs$ECGII , PeakData$RRCombined[26000:228000,] , QSwidth = QSwidth)

goodbeatslogical <- (QS$numvalues < as.numeric(quantile( QS$numvalues  , 0.95)))*(QS$numvalues > as.numeric(quantile( QS$numvalues  , (1-0.95)))) == 1

tmp <- QS$Value[ goodbeatslogical , ]
tmp<- apply(tmp , 1 , function(X){X - mean(X , na.rm = TRUE)})

m <- apply(tmp , 1 , function(X){mean(X , na.rm = TRUE)})
v <- apply(tmp , 1 , function(X){var(X , na.rm = TRUE)})

#lines( (1:length(m))/length(m) , m , col = 'yellow')
#lines( (1:length(m))/length(m), m + 2*sqrt(v) , col = 'red')
#lines( (1:length(m))/length(m), m - 2*sqrt(v) , col = 'red')

plot( 0  , 0  , xlim = c(-0,1) , ylim = ylim , xlab='t' , ylab ='Hz'  , type ='l')
title('T-P Waves Distribution')
for(i in 1:dim(QS$Date)[1]){
  if(QS$numvalues[i ] > as.numeric(quantile( QS$numvalues  , 0.8)) ){next}
  if(QS$numvalues[i ] < as.numeric(quantile( QS$numvalues  , 0.2)) ){next}
  tmp = 1:(min(dim(QS$Date)[2] , QS$numvalues[i ]) -1)
  lines((1:length(QS$Date[i,tmp]))/length(QS$Date[i,tmp]), QS$Value[i,tmp] - mean(QS$Value[i,tmp]) , col= rgb(i/dim(QS$Date)[1] , 0 , 0 , alpha = 0.01))
}
}

{dev.off()
  x11(20,20)
  par(mfrow = c(1 , 1))
  QS <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = PeakData$RRCombined[30000:31000,] , QSwidth = QSwidth)
  plot( 0  , 0  , xlim = c(-0,1) , ylim = ylim , xlab='t' , ylab ='Hz'  , type ='l')
  title('T-P Waves Distribution in AF')
  goodbeatslogical <- (QS$numvalues < as.numeric(quantile( QS$numvalues  , 0.95)))*(QS$numvalues > as.numeric(quantile( QS$numvalues  , (1-0.95)))) == 1
  
  tmp <- QS$Value[ goodbeatslogical , ]
  tmp<- apply(tmp , 1 , function(X){X - mean(X , na.rm = TRUE)})
  
  m <- apply(tmp , 1 , function(X){mean(X , na.rm = TRUE)})
  v <- apply(tmp , 1 , function(X){var(X , na.rm = TRUE)})
  
  for(i in 1:dim(QS$Date)[1]){
    if(QS$numvalues[i ] > as.numeric(quantile( QS$numvalues  , 0.95)) ){next}
    if(QS$numvalues[i ] < as.numeric(quantile( QS$numvalues  , 0.5)) ){next}
    tmp = 1:(min(dim(QS$Date)[2] , QS$numvalues[i ]) -1)
    lines((1:length(QS$Date[i,tmp]))/length(QS$Date[i,tmp]), QS$Value[i,tmp] - mean(QS$Value[i,tmp]) , col= rgb(1 , 0 , 0 , alpha = 0.01))
  }
  
  #lines( (1:length(m))/length(m) , m , col = 'blue')
  #lines( (1:length(m))/length(m), m + 2*sqrt(v) , col = 'black')
  #lines( (1:length(m))/length(m), m - 2*sqrt(v) , col = 'black')
  
  #plot( 0  , 0  , xlim = c(-0,1) , ylim = ylim , xlab='t' , ylab ='Hz'  , type ='l')
  
  #title('T-P Second order Statistics')
  #lines( (1:length(m))/length(m) , m , col = 'blue')
  #lines( (1:length(m))/length(m), m + 2*sqrt(v) , col = 'black')
  #lines( (1:length(m))/length(m), m - 2*sqrt(v) , col = 'black')
  
  QS <- AFD_ExtractAllSQ(ECGs$ECGII , PeakData$RRCombined[1:1000,] , QSwidth = QSwidth)
  
  goodbeatslogical <- (QS$numvalues < as.numeric(quantile( QS$numvalues  , 0.95)))*(QS$numvalues > as.numeric(quantile( QS$numvalues  , (1-0.95)))) == 1
  
  tmp <- QS$Value[ goodbeatslogical , ]
  tmp<- apply(tmp , 1 , function(X){X - mean(X , na.rm = TRUE)})
  
  m <- apply(tmp , 1 , function(X){mean(X , na.rm = TRUE)})
  v <- apply(tmp , 1 , function(X){var(X , na.rm = TRUE)})
  
  #lines( (1:length(m))/length(m) , m , col = 'yellow')
  #lines( (1:length(m))/length(m), m + 2*sqrt(v) , col = 'red')
  #lines( (1:length(m))/length(m), m - 2*sqrt(v) , col = 'red')
  
 # plot( 0  , 0  , xlim = c(-0,1) , ylim = ylim , xlab='t' , ylab ='Hz'  , type ='l')
#  title('T-P Waves Distribution')
  for(i in 1:dim(QS$Date)[1]){
    if(QS$numvalues[i ] > as.numeric(quantile( QS$numvalues  , 0.8)) ){next}
    if(QS$numvalues[i ] < as.numeric(quantile( QS$numvalues  , 0.2)) ){next}
    tmp = 1:(min(dim(QS$Date)[2] , QS$numvalues[i ]) -1)
    lines((1:length(QS$Date[i,tmp]))/length(QS$Date[i,tmp]), QS$Value[i,tmp] - mean(QS$Value[i,tmp]) , col= rgb(0 , 0 , 1 , alpha = 0.01))
  }
}
