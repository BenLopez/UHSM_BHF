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

PsimulatorFunction <- function( x , t_observation ){
  
  Pcen <- x[1]
  Pwidth <- x[2]
  Pcen2 <- x[3]
  Pwidth2 <- x[4]
  
  return( ECGSim_Gaussian( x = t_observation , mu = Pcen , sigma = Pwidth ) + ECGSim_Gaussian( x = t_observation , mu = Pcen2 , sigma = Pwidth2 ) )
}

#BC_PlotPairs(PriorNonImplausibleSet)

ModelDiscrepancy <- function(x , t_observation, PsimulatorFunction){
  1/( PsimulatorFunction( x , t_observation ) + 0.5 )
}

CalculateImplausability <- function( t_observation , x ,  z  , mdfun = ModelDiscrepancy){
  H = PWaveHM_CreateDesignMatrix(t_observation , x , PsimulatorFunction)
  Beta = PWaveHM_CalculateBetas(H , z)
  Im = mean( abs(z - H%*%Beta) / sqrt(ModelDiscrepancy(x , t_observation, PsimulatorFunction)))
  return( Im )
} 

PWaveHM_CreateDesignMatrix <- function(t_observation , x , PsimulatorFunction){
  t_observation <- as.matrix(t_observation)
  return( cbind( matrix(1 , dim(t_observation)[1] , 1 ) , t_observation , t_observation^2, as.matrix(PsimulatorFunction(x , t_observation))  ) )
}
FilestoProcess = DP_ChooseECGstoProcess()

Xstar = seq(0.5 ,1 , 0.01)
PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.001 , 0.95, 0.001 ) , b = c(0.6  , 0.04 , 0.6 , 0.04 ) , numbersamples = 1000000 )
PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.1 ,  ]
PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.02 ,  ]
PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]


RLDiff <-matrix(0 , length(listAllPatients) , 1)
RLWidthDiff <-matrix(0 , length(listAllPatients) , 1)
PWidth <- matrix(0 , length(listAllPatients) , 1)
PAmplitudeAF <- matrix(0 , length(listAllPatients) , 1)
ImplausabilityStruct <- matrix(0 , length(listAllPatients) , 1)
BetaStruct <- matrix(0 , length(listAllPatients) , 4)
AFLogicalStruct <- matrix(0 , length(listAllPatients) , 1)
XminStruct <- matrix(0 , length(listAllPatients) , 4)
mQSStruct <- matrix(0 , length(Xstar) , length(listAllPatients))
vQSStruct <- matrix(0 , length(Xstar) , length(listAllPatients))
PAmplitudeAFHat <- PAmplitudeAF
Ponoff <- matrix(0 , length(listAllPatients) , 2)
QSwidth = 10

EmulatorParameters <- PWaveHM_CreateDefaultEmulationclass()
for( k in 1:length(listAllPatients) ){
  
PatientID <- listAllPatients[k]

if(DP_CheckIfAFPatient( DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientID))){
  AFLogicalStruct[k] = 1
}

if(!DP_CheckECGreducedfilesprocessed(path = path , PatientsId = PatientID , Filestoprocess = 'ECGI_reduced')){next}
ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)

if(dim(RPeakData$RRCombined)[1] < 2000){next}

QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[1000:1500,] , QSwidth = QSwidth)

if(dim(as.matrix(QS_Struct$Date))[1] < 250){next}
if(length(QS_Struct$numvalues) == 0){next}
if(length(QS_Struct$Date) == 0){next}

EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct = QS_Struct , EmulatorParameters = EmulatorParameters ,Xstar =  Xstar )
mQS <- apply(EmulatedQS , 2 , median) 
mQS <- mQS - mean( mQS )
vQS <- apply(EmulatedQS , 2 , IQR)

Implausability <- PWaveHM_CalulateImplausabilityTQSegment( Xstar , mQS ,   PriorNonImplausibleSet , ImplausabilityFunction = CalculateImplausability )

#BC_PlotPairsFromTwoVariables(PriorNonImplausibleSet[sample(1:dim(PriorNonImplausibleSet)[1] , 10000) , ] , PriorNonImplausibleSet[Implausability < 50 , ])

ImplausabilityStruct[k] <- min(Implausability)
Xmin <- PriorNonImplausibleSet[which.min(Implausability) , ]  

H <- PWaveHM_CreateDesignMatrix( Xstar , Xmin , PsimulatorFunction )
Beta <- PWaveHM_CalculateBetas( H , mQS )
E_Y <- H%*%Beta

mQSStruct[, k] <- mQS
vQSStruct[,k] <- vQS
XminStruct[k,] <- Xmin
BetaStruct[k, ] <- Beta
RLDiff[k] <- abs(Xmin[1] - Xmin[3])
RLWidthDiff[k] <- Xmin[2] - Xmin[4]
PWidth[k] <- sum(PsimulatorFunction(Xmin , Xstar) > 0.01)
PAmplitudeAF[k] <- E_Y[which.max( PsimulatorFunction(Xmin , Xstar))]
PAmplitudeAFHat[k] <- BetaStruct[k , 4]*PsimulatorFunction(XminStruct[k,] , Xstar)[which.max(PsimulatorFunction(XminStruct[k,] , Xstar) )]
Ponoff[k,] <- which(abs(diff(PsimulatorFunction(XminStruct[2,] , Xstar) > 0.1)) ==1 )[1:2]

DP_WaitBar(k/length(listAllPatients))

}

VPAmplitude <- PAmplitudeAF
for(k in 1:length(PAmplitudeAF)){
if(BetaStruct[k , 4] != 0){
}
  VPAmplitude[k] <- max(vQSStruct[1:which.max(PsimulatorFunction(XminStruct[k,] , Xstar)) , k]) 
}

# Set up a P-wave simulator
BC_PlotCompareSingleHists(VPAmplitude[((AFLogicalStruct == 0)*(VPAmplitude!=0)) == 1 ] , VPAmplitude[((AFLogicalStruct == 1)*(VPAmplitude!=0)) == 1 ] , main = 'Observed P-Amplitude')

BC_PlotCompareSingleHists(X= PAmplitudeAFHat[((AFLogicalStruct == 0)*(PAmplitudeAFHat!=0)) == 1 ] , Y = PAmplitudeAFHat[((AFLogicalStruct == 1)*(PAmplitudeAFHat!=0)) == 1 ] , main = 'Predicted P-Amplitude')
BC_PlotCompareSingleHists(PAmplitudeAF[((AFLogicalStruct == 0)*(PAmplitudeAF!=0)) == 1 ] , PAmplitudeAF[((AFLogicalStruct == 1)*(PAmplitudeAF!=0)) == 1 ] , main = 'Observed P-Amplitude')
BC_PlotCompareSingleHists(PWidth[((AFLogicalStruct == 0)*(PWidth!=0)) == 1 ] , PWidth[((AFLogicalStruct == 1)*(PWidth!=0)) == 1 ] , main = 'P Width')
BC_PlotCompareSingleHists(XminStruct[((AFLogicalStruct == 0)*(XminStruct[,1]!=0)) == 1  , 1] , XminStruct[((AFLogicalStruct == 1)*(XminStruct[,1]!=0)) == 1  , 1] , main = 'X1')
BC_PlotCompareSingleHists(XminStruct[((AFLogicalStruct == 0)*(XminStruct[,1]!=0)) == 1  , 2] , XminStruct[((AFLogicalStruct == 1)*(XminStruct[,1]!=0)) == 1  , 2] , main = 'X2')
BC_PlotCompareSingleHists(XminStruct[((AFLogicalStruct == 0)*(XminStruct[,1]!=0)) == 1  , 3] , XminStruct[((AFLogicalStruct == 1)*(XminStruct[,1]!=0)) == 1  , 3] , main = 'X3')
BC_PlotCompareSingleHists(XminStruct[((AFLogicalStruct == 0)*(XminStruct[,1]!=0)) == 1  , 4] , XminStruct[((AFLogicalStruct == 1)*(XminStruct[,1]!=0)) == 1  , 4] , main = 'X4')


FullDataMatrix <-  cbind(XminStruct[ ,] , as.matrix(PAmplitudeAF) , as.matrix(PWidth) , as.matrix(PAmplitudeAFHat) )
rownames(FullDataMatrix) <- listAllPatients
tmp <- FullDataMatrix[(AFLogicalStruct) == 0 ,]
tmp <- tmp[tmp[,1] >0 ,]

tmp2 <- FullDataMatrix[(AFLogicalStruct) == 1 ,]
tmp2 <- tmp2[tmp2[,1] >0 ,]

BC_PlotPairsFromTwoVariables( tmp2 , tmp[sample(1:dim(tmp)[1] , 53),] , alpha = 0.2)



for(PatientID in  listAllPatients[1]){
#  PatientID <- NAFPatients[sample(1:length(NAFPatients) , 1)]
if(!DP_CheckECGreducedfilesprocessed(path = path , PatientsId = PatientID , Filestoprocess = 'ECGI_reduced')){next}
ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)

if(dim(RPeakData$RRCombined)[1] < 2000){next}

QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[1000:1500,] , QSwidth = QSwidth)

if(dim(QS_Struct$Date)[1] < 250){next}
if(length(QS_Struct$numvalues) == 0){next}
if(length(QS_Struct$Date) == 0){next}

EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct = QS_Struct , EmulatorParameters = EmulatorParameters ,Xstar =  Xstar )
mQS <- apply(EmulatedQS , 2 , median) 
mQS <- mQS - mean(mQS)

Implausability <- PWaveHM_CalulateImplausabilityTQSegment( Xstar , mQS ,   PriorNonImplausibleSet , ImplausabilityFunction = CalculateImplausability )
Xmin <- PriorNonImplausibleSet[which.min(Implausability) , ]  

H <- PWaveHM_CreateDesignMatrix( Xstar , Xmin , PsimulatorFunction )
Beta <- PWaveHM_CalculateBetas( H , mQS )
E_Y <- H%*%Beta

PWaveHM_PlotMatch(Xstar , as.matrix(H[,dim(H)[2]])%*%Beta[dim(H)[2]] , mQS - (H[,1:(dim(H)[2] -1 )]%*%Beta[1:(dim(H)[2] -1 ),]) , ModelDiscrepancy(Xmin , Xstar, PsimulatorFunction))
lines(Xstar ,  Beta[(dim(H)[2])]*ECGSim_Gaussian( x = Xstar , mu = Xmin[1]  , sigma = Xmin[2] ) , col = 'black')
lines(Xstar ,  Beta[(dim(H)[2])]*ECGSim_Gaussian( x = Xstar , mu = Xmin[3] , sigma = Xmin[4] ) , col = 'black')
abline(v = Xstar[which.max( PsimulatorFunction(Xmin , Xstar))])

title(paste0(PatientID , ' Non-ImplausibleMatch Im = ' , min(Implausability)))
}

{
PatientID <- DP_choosepatient(listAllPatients = listAllPatients)
ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)
QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[1000:1500,] , QSwidth = QSwidth)
EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct = QS_Struct , EmulatorParameters = EmulatorParameters ,Xstar =  Xstar )
mQS <- apply(EmulatedQS , 2 , median) 
mQS <- mQS - mean( mQS )
vQS <- apply(EmulatedQS , 2 , IQR)

Implausability <- PWaveHM_CalulateImplausabilityTQSegment( Xstar , mQS ,   PriorNonImplausibleSet , ImplausabilityFunction = CalculateImplausability )
Xmin <- PriorNonImplausibleSet[which.min(Implausability) , ]  

H <- PWaveHM_CreateDesignMatrix( Xstar , Xmin , PsimulatorFunction )
Beta <- PWaveHM_CalculateBetas( H , mQS )
E_Y <- H%*%Beta

{patientnumber = 26
plot(DP_NormaliseData(mQSStruct[,patientnumber]) , type ='l' , ylim = c(-3,3))
points(DP_NormaliseData(vQSStruct[,patientnumber])  , type ='l' , col = 'red')
listAllPatients[patientnumber]}

PWaveHM_PlotMatch(Xstar , as.matrix(H[,dim(H)[2]])%*%Beta[dim(H)[2]] , mQS - (H[,1:(dim(H)[2] -1 )]%*%Beta[1:(dim(H)[2] -1 ),]) , ModelDiscrepancy(Xmin , Xstar, PsimulatorFunction))
lines(Xstar ,  Beta[(dim(H)[2])]*ECGSim_Gaussian( x = Xstar , mu = Xmin[1]  , sigma = Xmin[2] ) , col = 'black')
lines(Xstar ,  Beta[(dim(H)[2])]*ECGSim_Gaussian( x = Xstar , mu = Xmin[3] , sigma = Xmin[4] ) , col = 'black')
abline(v = Xstar[which.max( PsimulatorFunction(Xmin , Xstar))])

title(paste0(PatientID , ' Non-ImplausibleMatch Im = ' , min(Implausability)))
}


# Second order classification for 



FullDataMatrix <-  cbind(XminStruct[ ,] , as.matrix(PAmplitudeAF) , as.matrix(PWidth) )
tmp <- FullDataMatrix[which((apply(FullDataMatrix ,1 , function(X){sum(X == 0) }) != dim(FullDataMatrix)[2])) ,  ] 
tmpAFLogical <- AFLogicalStruct[apply(FullDataMatrix[,] ,1 , function(X){sum(X == 0) }) !=dim(FullDataMatrix)[2]]
m_PWaveAF <- as.matrix( apply(tmp[ tmpAFLogical ==1 , ], 2 , mean))
C_PWaveAF <- cov((tmp[ tmpAFLogical ==1 , ]))
m_PWaveNAF <- as.matrix(apply(tmp[ tmpAFLogical==0 , ] , 2 , mean))
C_PWaveNAF <- cov((tmp[ tmpAFLogical==0 , ]))

{par(mfrow = c(1,2))
Im1 <- apply(FullDataMatrix[which((apply(FullDataMatrix ,1 , function(X){sum(X == 0) }) != dim(FullDataMatrix)[2])) ,  ]  , 1 , function(X){ t(X - m_PWaveAF)%*%solve(DP_AddNugget(C_PWaveAF , 0.0000001) )%*%(X - m_PWaveAF) }  )

plot(tmpAFLogical , log(Im1)  , pch = 16 , col = rgb(1,0,0, alpha = 0.05) , xlab = c('AF Indicator') , ylab = c('Implausability of not being AF'))
abline(2.2 , 0)
title('Log Multivariate Implausability')

Im2 <- apply(FullDataMatrix[which((apply(FullDataMatrix ,1 , function(X){sum(X == 0) }) != dim(FullDataMatrix)[2])) ,  ]  , 1 , function(X){ t(X - m_PWaveNAF)%*%solve(DP_AddNugget(C_PWaveNAF) )%*%(X - m_PWaveNAF) }  )

plot(tmpAFLogical , log(Im2)  , pch = 16 , col = rgb(1,0,0, alpha = 0.05) , xlab = c('AF Indicator') , ylab = c('Implausability of being AF'))
title('Log Multivariate Implausability')}


{x11()
par(mfrow = c(2 , 2))
plot(Xstar,apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) , type = 'l'  , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 0] , 1 , mean ) , type = 'l'  , col = rgb(1,0,0 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
title('Mean P-Waves')

plot(Xstar,apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , var ) , type = 'l'  , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 0] , 1 , var ) , type = 'l'  , col = rgb(1,0,0 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
title('Variance Mean P-Waves')

plot(Xstar,apply(vQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) , type = 'l'  , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(vQSStruct[ , AFLogicalStruct == 0] , 1 , mean ) , type = 'l'  , col = rgb(1,0,0 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
title('Mean Variance P-Waves')

plot(Xstar,apply(vQSStruct[ , AFLogicalStruct == 1] , 1 , var ) , type = 'l'  , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(vQSStruct[ , AFLogicalStruct == 0] , 1 , var ) , type = 'l'  , col = rgb(1,0,0 , alpha = 0.5) , xlab = 'time' , ylab = 'V')
title('Variance Variance P-Waves')}

{x11()
par(mfrow = c(1 , 1))
plot(Xstar,0*apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) , type = 'l'  , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'time' , ylab = 'V' , ylim = c(-40,40))
for(i in 1:dim(vQSStruct)[2]){
  if(AFLogicalStruct[i] == 1){
    lines(Xstar,mQSStruct[ , i]  , type = 'l'  , col = rgb(0,0,1 , alpha = 0.2) , xlab = 'time' , ylab = 'V')
  }
  if(AFLogicalStruct[i] == 0){
    lines(Xstar,mQSStruct[ , i]  , type = 'l'  , col = rgb(1,0,0 , alpha = 0.05) , xlab = 'time' , ylab = 'V')
  }
}
lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) , type = 'l'  , col = rgb(0,0,1 , alpha = 1) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) + 3*apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , function(X){sqrt(var(X))} ) , type = 'l'  , col = rgb(0,1,0 , alpha = 1) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) - 3*apply(mQSStruct[ , AFLogicalStruct == 1] , 1 , function(X){sqrt(var(X))} ), type = 'l'  , col = rgb(0,1,0 , alpha = 1) , xlab = 'time' , ylab = 'V')

lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 0] , 1 , mean ) , type = 'l'  , col = rgb(1,0,0 , alpha = 1) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 0] , 1 , mean ) + 3*apply(mQSStruct[ , AFLogicalStruct == 0] , 1 , function(X){sqrt(var(X))} ) , type = 'l'  , col = rgb(0,0,0 , alpha = 1) , xlab = 'time' , ylab = 'V')
lines(Xstar,apply(mQSStruct[ , AFLogicalStruct == 0] , 1 , mean ) - 3*apply(mQSStruct[ , AFLogicalStruct == 0] , 1 , function(X){sqrt(var(X))} ), type = 'l'  , col = rgb(0,0,0 , alpha = 1) , xlab = 'time' , ylab = 'V')

title('Pwaves AF (Blue) , NAF (Red)')
}


{x11()
  par(mfrow = c(1 , 1))
  plot(Xstar,0*apply(vQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) , type = 'l'  , col = rgb(0,0,1 , alpha = 0.5) , xlab = 'time' , ylab = 'V' , ylim = c(0,100))
  for(i in 1:dim(vQSStruct)[2]){
    if(AFLogicalStruct[i] == 1){
      lines(Xstar,vQSStruct[ , i]  , type = 'l'  , col = rgb(0,0,1 , alpha = 0.2) , xlab = 'time' , ylab = 'V')
    }
    if(AFLogicalStruct[i] == 0){
      lines(Xstar,vQSStruct[ , i]  , type = 'l'  , col = rgb(1,0,0 , alpha = 0.05) , xlab = 'time' , ylab = 'V')
    }
  }
  lines(Xstar,apply(vQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) , type = 'l'  , col = rgb(0,0,1 , alpha = 1) , xlab = 'time' , ylab = 'V')
  lines(Xstar,apply(vQSStruct[ , AFLogicalStruct == 1] , 1 , mean ) + 3*apply(vQSStruct[ , AFLogicalStruct == 1] , 1 , function(X){sqrt(var(X))} ) , type = 'l'  , col = rgb(0,1,0 , alpha = 1) , xlab = 'time' , ylab = 'V')

  lines(Xstar,apply(vQSStruct[ , AFLogicalStruct == 0] , 1 , mean ) , type = 'l'  , col = rgb(1,0,0 , alpha = 1) , xlab = 'time' , ylab = 'V')
  lines(Xstar,apply(vQSStruct[ , AFLogicalStruct == 0] , 1 , mean ) + 3*apply(vQSStruct[ , AFLogicalStruct == 0] , 1 , function(X){sqrt(var(X))} ) , type = 'l'  , col = rgb(0,0,0 , alpha = 1) , xlab = 'time' , ylab = 'V')

  title('Pwaves AF (Blue) , NAF (Red)')
}



x11(20,14)
parcoord(t(mQSStruct[ , which(AFLogicalStruct == 0)[sample(1:sum(AFLogicalStruct == 0) , 53)]] )  ,  rgb(1,0,0 , alpha = 0.2) )
parcoord(t(mQSStruct[ , AFLogicalStruct == 1])  ,  rgb(0,0,1 , alpha = 0.2),add = T)
title('Parallel Corrdinates Mean')

x11(20,14)
parcoord(t(vQSStruct[ , which(AFLogicalStruct == 0)[sample(1:sum(AFLogicalStruct == 0) , 53)]] )  ,  rgb(1,0,0 , alpha = 0.2) )
parcoord(t(vQSStruct[ , AFLogicalStruct == 1])  ,  rgb(0,0,1 , alpha = 0.2),add = T)
title('Parallel Corrdinates Variance')

variables = 20:25
BC_PlotPairsFromTwoVariables(t(mQSStruct[variables , which(AFLogicalStruct == 0)[sample(1:sum(AFLogicalStruct == 0) , 53)]]) , t(mQSStruct[variables , AFLogicalStruct == 1]) , alpha = 0.1 )

variables <- 1:(dim(mQSStruct)[2])
tmp <- mQSStruct[ , which((apply(mQSStruct ,2 , function(X){sum(X == 0) }) != length(Xstar))[variables]) ] 
tmpAFLogical <- AFLogicalStruct[apply(mQSStruct[,] ,2 , function(X){sum(X == 0) }) != length(Xstar)]
m_PWaveAF <- as.matrix(apply(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 1], 1 , mean))
C_PWaveAF <- cov(t(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 1]))
m_PWaveNAF <- as.matrix(apply(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 0] , 1 , mean))
C_PWaveNAF <- cov(t(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 0]))

tmp <- vQSStruct[ , which((apply(vQSStruct ,2 , function(X){sum(X == 0) }) != length(Xstar))[variables]) ] 
tmpAFLogical <- AFLogicalStruct[apply(vQSStruct[,] ,2 , function(X){sum(X == 0) }) != length(Xstar)]
mV_PWaveAF <- as.matrix(apply(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 1], 1 , mean))
CV_PWaveAF <- cov(t(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 1]))
mV_PWaveNAF <- as.matrix(apply(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 0] , 1 , mean))
CV_PWaveNAF <- cov(t(tmp[ , tmpAFLogical[1:dim(tmp)[2]] == 0]))

plot(mV_PWaveAF + BE_SampleGP(C_PWaveAF) , type ='l' , col = rgb(0,0,1 , alpha = 0.1) )
for(i in 1:100){
lines(mV_PWaveAF + BE_SampleGP(C_PWaveAF) , type ='l' , col = rgb(0,0,1 , alpha = 0.1))
}
for(i in 1:100){
  lines(mV_PWaveNAF + BE_SampleGP(C_PWaveNAF) , type ='l' , col = rgb(1,0,0 , alpha = 0.1))
}

{x11()
par(mfrow = c(2,2))
Im1 <- apply(mQSStruct[ , apply(mQSStruct ,2 , function(X){sum(X == 0) }) != length(Xstar)]  , 2 , function(X){ t(X - m_PWaveAF)%*%solve(DP_AddNugget(C_PWaveAF , 0.0000001) )%*%(X - m_PWaveAF) }  )

plot(tmpAFLogical , log(Im1)  , pch = 16 , col = rgb(1,0,0, alpha = 0.05) , xlab = c('AF Indicator') , ylab = c('Implausability of not being AF'))
abline(5 , 0)
title('Log Multivariate Implausability')

Im2 <- apply(mQSStruct[ , apply(mQSStruct ,2 , function(X){sum(X == 0) }) != length(Xstar)]  , 2 , function(X){ t(X - m_PWaveNAF)%*%solve(DP_AddNugget(C_PWaveNAF) )%*%(X - m_PWaveNAF) }  )

plot(tmpAFLogical , log(Im2)  , pch = 16 , col = rgb(1,0,0, alpha = 0.05) , xlab = c('AF Indicator') , ylab = c('Implausability of being AF'))
title('Log Multivariate Implausability')

Im3 <- apply(vQSStruct[ , apply(vQSStruct ,2 , function(X){sum(X == 0) }) != length(Xstar)]  , 2 , function(X){ t(X - mV_PWaveAF)%*%solve(DP_AddNugget(CV_PWaveAF , 0.0000001) )%*%(X - mV_PWaveAF) }  )

plot(tmpAFLogical , log(Im3)  , pch = 16 , col = rgb(1,0,0, alpha = 0.05) , xlab = c('AF Indicator') , ylab = c('Implausability of not being AF'))
abline(5 , 0)
title('Log Multivariate Implausability')

Im4 <- apply(vQSStruct[ , apply(vQSStruct ,2 , function(X){sum(X == 0) }) != length(Xstar)]  , 2 , function(X){ t(X - mV_PWaveNAF)%*%solve(DP_AddNugget(CV_PWaveNAF , 0.0000001) )%*%(X - mV_PWaveNAF) }  )
plot(tmpAFLogical , log(Im4)  , pch = 16 , col = rgb(1,0,0, alpha = 0.05) , xlab = c('AF Indicator') , ylab = c('Implausability of not being AF'))
title('Log Multivariate Implausability')
}

