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
  Im = max( abs(z - H%*%Beta) / sqrt(ModelDiscrepancy(x , t_observation, PsimulatorFunction)))
  return( Im )
} 

FilestoProcess = DP_ChooseECGstoProcess()

Xstar = seq(0.5 ,1 , 0.01)
PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.001 , 1, 0.001 ) , b = c(0.6  , 0.1 , 0.6 , 0.1 ) , numbersamples = 1000000 )
PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.2 ,  ]
PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.02 ,  ]
PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]

RLDiff <-matrix(0 , length(listAllPatients) , 1)
RLWidthDiff <-matrix(0 , length(listAllPatients) , 1)
PWidth <- matrix(0 , length(listAllPatients) , 1)
PAmplitudeAF <- matrix(0 , length(listAllPatients) , 1)
ImplausabilityStruct <- matrix(0 , length(listAllPatients) , 1)
BetaStruct <- matrix(0 , length(listAllPatients) , 3)
AFLogicalStruct <- matrix(0 , length(listAllPatients) , 1)
XminStruct <- matrix(0 , length(listAllPatients) , 4)

EmulatorParameters <- PWaveHM_CreateDefaultEmulationclass()
for( k in 1:length(listAllPatients) ){
  
PatientID <- listAllPatients[k]

if(DP_CheckIfAFPatient( DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientID))){
  AFLogicalStruct[k] = 1
}


ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = FilestoProcess  )
RPeakData <- DP_LoadRpeaksfile(path , PatientID)
QSwidth = 10


QS_Struct <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[1000:1500,] , QSwidth = QSwidth)

if(length(QS_Struct$numvalues) == 0){next}
if(length(QS_Struct$Date) == 0){next}

EmulatedQS <- PWaveHM_EmulateTQSegment( QS_Struct , EmulatorParameters = EmulatorParameters , Xstar )
mQS <- apply(EmulatedQS , 2 , median) 
mQS <- mQS - mean( mQS )

Implausability <- PWaveHM_CalulateImplausabilityTQSegment( Xstar , mQS ,   PriorNonImplausibleSet , ImplausabilityFunction = CalculateImplausability )

#BC_PlotPairsFromTwoVariables(PriorNonImplausibleSet[sample(1:dim(PriorNonImplausibleSet)[1] , 10000) , ] , PriorNonImplausibleSet[Implausability < 50 , ])

ImplausabilityStruct[k] <- min(Implausability)
Xmin <- PriorNonImplausibleSet[which.min(Implausability) , ]  

H <- PWaveHM_CreateDesignMatrix( Xstar , Xmin , PsimulatorFunction )
Beta <- PWaveHM_CalculateBetas( H , mQS )
E_Y <- H%*%Beta

XminStruct[k,] <- Xmin
BetaStruct[k, ] <- Beta
RLDiff[k] <- abs(Xmin[1] - Xmin[3])
RLWidthDiff[k] <- Xmin[2] - Xmin[4]
PWidth[k] <- sum(PsimulatorFunction(Xmin , Xstar) > 0.01)
PAmplitudeAF[k] <- E_Y[which.max( PsimulatorFunction(Xmin , Xstar))]

print(RLDiff[k])
print(PWidth[k])
DP_WaitBar(k/length(listAllPatients))

}
# Set up a P-wave simulator

PWaveHM_PlotMatch(Xstar , H%*%Beta , mQS , ModelDiscrepancy(Xmin , Xstar, PsimulatorFunction))
lines(Xstar , H[,1:2]%*%Beta[1:2,]+ Beta[3]*ECGSim_Gaussian( x = Xstar , mu = Xmin[1]  , sigma = Xmin[2] ) , col = 'black')
lines(Xstar , H[,1:2]%*%Beta[1:2,]+ Beta[3]*ECGSim_Gaussian( x = Xstar , mu = Xmin[3] , sigma = Xmin[4] ) , col = 'black')
abline(PAmplitudeAF , 0)
abline(v = Xstar[which.max( PsimulatorFunction(Xmin , Xstar))])

title(paste0( 'Non-ImplausibleMatch.'))
abline(E_Y[which.min(abs(Xstar - XminIm[1])) ] , 0)

