{
  if(file.exists('CheckforDefaultsScript.R')){
    source('CheckforDefaultsScript.R')
  }else{
    pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
    source("LibrariesAndSettings.R" , print.eval  = TRUE )
    DP_LoadPatientIndex()
    DP_ChooseDataReps()
    FilestoProcess <- DP_ChooseECGstoProcess() 
    HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
  }
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}


PatientID <- DP_choosepatient( listAllPatients )

ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = DP_ChooseECGstoProcess()  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)
QSwidth = 10

QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[1000:1500,] , QSwidth = QSwidth)

plot( 0  , 0  , xlim = c(0,1) , ylim = c(-40,40) , xlab='t' , ylab ='Hz'  , type ='l')
for(i in 1:dim(QS_Struct$Date)[1]){
  if(QS_Struct$numvalues[i ] > as.numeric(quantile( QS_Struct$numvalues  , 0.95)) ){next}
  if(QS_Struct$numvalues[i ] < as.numeric(quantile( QS_Struct$numvalues  , 0.5)) ){next}
  tmp = 1:(min(dim(QS_Struct$Date)[2] , QS_Struct$numvalues[i ]) -1)
  lines((1:length(QS_Struct$Date[i,tmp]))/length(QS_Struct$Date[i,tmp]), QS_Struct$Value[i,tmp] - mean(QS_Struct$Value[i,tmp]) , col= rgb(1 , 0 , 0 , alpha = 0.01))
}


Xstar = seq(0.5 ,1 , 0.01)
PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 1, 0.001 ) ,b = c(0.6  , 0.1 ) , numbersamples = 10000)

PAmplitude <- PWaveHM_EmulateEstimatePAmplitude(QS_Struct , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar , PriorNonImplausibleSet)





BC_PlotPairsFromTwoVariables(X = PriorNonImplausibleSet[which(Implausability > 1) , ] , Y = PriorNonImplausibleSet[which(Implausability < 1) , ]  , alpha = 0.1)


# Set up a P-wave simulator



PWaveHM_PlotMatch(Xstar , H%*%Beta , mQS , 2)
title(paste0( 'Non-ImplausibleMatch. P-Amplitude = ' , E_Y[which.min(abs(Xstar - XminIm[1])) ]))
abline(E_Y[which.min(abs(Xstar - XminIm[1])) ] , 0)

