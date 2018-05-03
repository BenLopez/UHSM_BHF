pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

path_PatIndex = choose.files(caption="Select 2017 PatientIndex.csv file")

if(length(path_PatIndex)>0){
  PatIndex2017 = read.csv(file=path_PatIndex, stringsAsFactors = FALSE)
  # PatIndex2017$FirstITUEntry=as.POSIXct(PatIndex2017$FirstITUEntry, format="%d/%m/%Y %H:%M")
  # PatIndex2017$LastITUEntry=as.POSIXct(PatIndex2017$LastITUEntry, format="%d/%m/%Y %H:%M")
} else {
  warning("No Patient Info provided")
  sub_pat = list()
}

path = choose.dir(caption="Select folder containing data repository")
listAllPatients = list.dirs(path = path, full.names = FALSE, recursive = FALSE)
subList = select.list(listAllPatients, preselect = NULL, multiple = TRUE, title = NULL, graphics = TRUE )


j <- 1
plot(DiscreteData[[j]]$tt - DiscreteData[[j]]$tt[1] , c((cumsum(DiscreteData[[j]]$HeartRate > 130))) , type ='l' , ylim = c(0,500) , xlab ='Time (seconds) from Beginning of Recording' , ylab ='Number of Intervals With HeartRate over 130')
title('Discrete Data Analysis Heart Rate Over 130')
for(j in c(1:length(DiscreteData)) )
{
  sub_pat = subset(PatIndex2017, PseudoId %in% DiscreteDataPatientCode[[j]])
  if(sub_pat$Pre_OperativeHeartRhythm != "Sinus Rhythm"){next}
  if(is.na(TimeNAF[[j]])){lines(DiscreteData[[j]]$tt - DiscreteData[[j]]$tt[1] , c((cumsum(DiscreteData[[j]]$HeartRate > 130))), col = 'black')}
  
  if(!is.na(TimeNAF[[j]])){
    lines(DiscreteData[[j]]$tt[DiscreteData[[j]]$tt<TimeNAF[[j]]] - DiscreteData[[j]]$tt[1] , cumsum(DiscreteData[[j]]$HeartRate > 130 )[DiscreteData[[j]]$tt<TimeNAF[[j]]], col = 'red')
    lines(difftime(DiscreteData[[j]]$tt[DiscreteData[[j]]$tt>TimeNAF[[j]]] , DiscreteData[[j]]$tt[1] , units = 'secs'), cumsum(DiscreteData[[j]]$HeartRate > 130 )[DiscreteData[[j]]$tt>TimeNAF[[j]]], col = 'blue')
    }

  #title(DiscreteDataPatientCode[[j]])
#abline(200,0 , col = 'blue')
#abline(100,0 , col = 'red')
}
legend( 250000, 450, legend=c("No Diagnosed AF", "Before AF Diagnosed" , "After AF Diagnosed"),
       col=c("Black","red", "blue"), lty=rep(1 , 3), cex=0.8)


num_over130 <- matrix(0 , length(DiscreteData) , 1 )
lengthoftimeseries <- num_over130

for(i in c(1:length(DiscreteData)) )
{
num_over130[i] <- sum( DiscreteData[[i]]$HeartRate > 130 )/length(DiscreteData[[i]]$HeartRate)
lengthoftimeseries[i] <- length(DiscreteData[[i]]$HeartRate)
}

setforgoodpopulation <- num_over130 < mean(num_over130[ !is.na(TimeNAF)])

MatrixforSecondOrderCalulations <- matrix(0 , length(lengthoftimeseries[setforgoodpopulation]) , max(lengthoftimeseries[setforgoodpopulation]))
MatrixforSecondOrderCalulations2 <- MatrixforSecondOrderCalulations
counter <- 1

for(i in 1:length(DiscreteData) )
{
  if(setforgoodpopulation[i] == FALSE){next}
  if(!is.na(TimeNAF[[i]])){next}
  MatrixforSecondOrderCalulations2[ counter , 1:length( DiscreteData[[i]]$tt )] <- ( ( cumsum(DiscreteData[[i]]$HeartRate)/cumsum(DiscreteData[[i]]$HeartRate>0) ))
  MatrixforSecondOrderCalulations[ counter , 1:length( DiscreteData[[i]]$tt )] <- DiscreteData[[i]]$HeartRate - MatrixforSecondOrderCalulations2[ counter , 1:length( DiscreteData[[i]]$tt )]
  counter <- counter +1
}

n <- apply(MatrixforSecondOrderCalulations2 != 0 , 2 , sum )
MatrixforSecondOrderCalulations <- MatrixforSecondOrderCalulations[ , n > 2]
MatrixforSecondOrderCalulations2 <- MatrixforSecondOrderCalulations2[ , n > 2]
n <- n[n>2]

Sumf <- apply(MatrixforSecondOrderCalulations , 2 , sum )
SumSquaresf <- apply(MatrixforSecondOrderCalulations^2 , 2 , sum )

Sumf2 <- apply(MatrixforSecondOrderCalulations2 , 2 , sum )
SumSquaresf2 <- apply(MatrixforSecondOrderCalulations2^2 , 2 , sum )

Sumf[is.na(Sumf)]<- 0
SumSquaresf[is.na(SumSquaresf)]<-0
n[n < 2]<- 2
n[is.na(n)]<-2

fbar <- imfilter1D( Sumf/(n) , rep(1/1000 , 1000 ) )
Sigma_f <- imfilter1D( n/(n-1)*((SumSquaresf/n) - fbar^2) , rep(1/1000 , 1000 ) )

fbar2 <- Sumf2/n
Sigma_f2 <- n/(n-1)*( (SumSquaresf2/n) - (fbar2^2) )

par(mfrow = c(2 , 1))
plot( fbar , type ='l' , ylim = c(-100,100) , xlab = 'time from start' , ylab = 'Heart-rate' )
lines( fbar + 2*sqrt(Sigma_f + Sigma_f/n) , col = 'blueviolet' )
lines( fbar - 2*sqrt(Sigma_f + Sigma_f/n) , col = 'blueviolet' )
title('Heart-rate Local Behaviour')
j <- 25
Temp <- DiscreteData[[j]]$HeartRate - ( cumsum(DiscreteData[[j]]$HeartRate)/cumsum(DiscreteData[[j]]$HeartRate>0) )
if(!is.na(TimeNAF[[j]])){
lines(1:length( DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt < TimeNAF[[j]]]) ,  Temp[DiscreteData[[j]]$tt<TimeNAF[[j]]] , col = 'red')
lines((length(DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt<TimeNAF[[j]]]) + 1):length(DiscreteData[[j]]$HeartRate) , Temp[DiscreteData[[j]]$tt>TimeNAF[[j]]], col = 'blue')
}
if(is.na(TimeNAF[[j]])){ lines(  Temp , col = 'red')}
  
legend( 80000, 100, legend=c("Mean Pop with No Diagnosed AF", "2Sigma Credible Interval" , "Before AF Diagnosed" , "After AF Diagnosed"),
        col=c("Black" , "Yellow","red", "blue"), lty=rep(1 , 3), cex=0.8)
plot( fbar2 , type ='l' , ylim = c(0,200) , xlab = 'time from start' , ylab = 'Heart-rate' )
lines( fbar2 + 2*sqrt(Sigma_f2 + Sigma_f2/n) , col = 'blueviolet' )
lines( fbar2 - 2*sqrt(Sigma_f2 + Sigma_f2/n) , col = 'blueviolet' )
Temp <- ( cumsum(DiscreteData[[j]]$HeartRate)/cumsum(DiscreteData[[j]]$HeartRate>0) )
if(!is.na(TimeNAF[[j]])){
  lines(1:length( DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt < TimeNAF[[j]]]) ,  Temp[DiscreteData[[j]]$tt<TimeNAF[[j]]] , col = 'red')
  lines((length(DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt<TimeNAF[[j]]]) + 1):length(DiscreteData[[j]]$HeartRate) , Temp[DiscreteData[[j]]$tt>TimeNAF[[j]]], col = 'blue')
}
if(is.na(TimeNAF[[j]])){ lines(  Temp , col = 'red')}
title('Heart-rate Global Behaviour')
legend( 80000,200, legend=c("Mean Pop with No Diagnosed AF", "2Sigma Credible Interval" , "Before AF Diagnosed" , "After AF Diagnosed"),
        col=c("Black" , "Yellow","red", "blue"), lty=rep(1 , 3), cex=0.8)

dev.off()
par(mfrow = c(2 , 1))
j <- 1
Temp <- DiscreteData[[j]]$HeartRate - ( cumsum(DiscreteData[[j]]$HeartRate)/cumsum(DiscreteData[[j]]$HeartRate>0) )
stdresid<- c(abs(Temp - fbar[1:length(Temp)])/sqrt(Sigma_f[1:length(Temp)]))
plot(stdresid , type ='l' , col ='red' , ylim = c(0,15))
if(!is.na(TimeNAF[[j]])){
  lines((length(DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt<TimeNAF[[j]]]) + 1):length(DiscreteData[[j]]$HeartRate) , stdresid[DiscreteData[[j]]$tt>TimeNAF[[j]]], col = 'black')
  lines(1:length( DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt < TimeNAF[[j]]]) ,  stdresid[DiscreteData[[j]]$tt<TimeNAF[[j]]] , col = 'red')
  }
abline(2,0 )
title(paste0( TimeNAF[[j]] , ' ' , DiscreteDataPatientCode[[j]]))
Temp <- ( cumsum(DiscreteData[[j]]$HeartRate)/cumsum(DiscreteData[[j]]$HeartRate>0) )
stdresid<- c(abs(Temp - fbar2[1:length(Temp)])/sqrt(Sigma_f2[1:length(Temp)]))
plot( stdresid   , type ='l', col ='red', ylim = c(0,15))
if(!is.na(TimeNAF[[j]])){
  lines((length(DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt<TimeNAF[[j]]]) + 1):length(DiscreteData[[j]]$HeartRate) , stdresid[DiscreteData[[j]]$tt>TimeNAF[[j]]], col = 'black')
  lines(1:length( DiscreteData[[j]]$HeartRate[DiscreteData[[j]]$tt < TimeNAF[[j]]]) ,  stdresid[DiscreteData[[j]]$tt<TimeNAF[[j]]] , col = 'red')
}
abline(2,0 )



j <- 1
par(mfrow = c(1 , 1))
Temp <- DiscreteData[[j]]$HeartRate - ( cumsum(DiscreteData[[j]]$HeartRate)/cumsum(DiscreteData[[j]]$HeartRate>0) )
stdresid<- c(abs(Temp - fbar[1:length(Temp)])/sqrt(Sigma_f[1:length(Temp)]))
stdresid[is.na(stdresid)]<-0
plot(imfilter1D(stdresid>2, c(rep(1,100) , rep(0,100))) )
title(paste0( TimeNAF[[j]] , ' ' , DiscreteDataPatientCode[[j]]))
j<-j+1