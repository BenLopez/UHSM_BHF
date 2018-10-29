DataSetPriorProbabilities <- BC_ExtractValidationPriorsLocalPrediction(Priorprobabilities , DataBaseMaster , DataBase)

DataSetPriorProbabilities$B =  0.0581065165409433
DataSetPriorProbabilities$`B^c` = 0.941893483459057

ProbabilityStuct <- matrix( 0 , 1 , 3)
AFMeanStruct <- matrix(0 , 1 , 2)
NAFMeanStruct <- matrix(0 , 1 , 2)
counter <-1


for(ii in 1:length(listAllPatients)){
  # Check file exists
  if(DP_CheckDistributionSummariesExists(path , listAllPatients[[ii]]) == FALSE){next}
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017  , listAllPatients[[ii]])
  t <- DP_LoadRpeaksfile(path , listAllPatients[ii])$RRCombined$t
  # Some logic to avoid dodgy data
  if(length(t) < 5000){next}
  
  DistributionSummaries <- DP_LoadDistributionSummaries(path , listAllPatients[[ii]])
  AdjustedBeliefs <- BC_EstimateGlobalandLocalParametersDisSum( DistributionSummaries )
  
  if( DP_CheckIfAFPatient( MetaData ) ){
    
    AdjustedBeliefs$W <- AdjustedBeliefs$W[BC_CreateAFAnnotationFomMetaData(t , MetaData) == 0, ]
    Probabilities <- BC_CalulateDenistiesGMM( W = DP_RemoveNaRows(AdjustedBeliefs$W) , LocalDistributionStruct = LocalDistributionStruct )[ , 2:3]
    #Probabilities <- (DataSetPriorProbabilities$B*Probabilities[ , 2]) / (DataSetPriorProbabilities$B*Probabilities[ , 2] + DataSetPriorProbabilities$`B^c`*Probabilities[ , 1] ) 
    tmp <- cbind(as.matrix(Probabilities) , matrix(1 , dim(Probabilities)[1] , 1))
    ProbabilityStuct <- rbind(ProbabilityStuct , tmp)
    }
  if( DP_CheckIfAFPatient( MetaData ) == FALSE ){
    
    AdjustedBeliefs$W <- AdjustedBeliefs$W
    
    Probabilities <- BC_CalulateDenistiesGMM( W = DP_RemoveNaRows(AdjustedBeliefs$W) , LocalDistributionStruct = LocalDistributionStruct )[ , 2:3]
    #Probabilities <- (DataSetPriorProbabilities$B*Probabilities[ , 2]) / (DataSetPriorProbabilities$B*Probabilities[ , 2] + DataSetPriorProbabilities$`B^c`*Probabilities[ , 1] ) 
    tmp <- cbind(as.matrix(Probabilities) , matrix(0 , dim(Probabilities)[1] , 1))
    ProbabilityStuct <- rbind(ProbabilityStuct , tmp)
    }
  if(ii == 1){ProbabilityStuct <- ProbabilityStuct[-1,]}
  DP_WaitBar(ii/length(listAllPatients))
}

{EmulatorX <- matrix(0, 1 ,  2 )
EmulatorY <- matrix(0, 1 ,  1 )
EmulatorSigma_n <- matrix(0, 1 ,  1 )
setofbeta = randomLHS(100, 1)
for(beta in setofbeta){

PosteriorProbabilityStructure <- rbind(ProbabilityStuct[sample(which(ProbabilityStuct[ , 3] == 1) , round(beta*1000000)) , ] , ProbabilityStuct[sample(which(ProbabilityStuct[ , 3] == 0) , round((1-beta)*1000000)) , ])
PosteriorProbabilityStructure <- cbind(beta*PosteriorProbabilityStructure[,2] /(beta*PosteriorProbabilityStructure[,2] + (1-beta)*PosteriorProbabilityStructure[,1]) , PosteriorProbabilityStructure[,3])
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput( BC_CreateCalibrationStructure(PosteriorProbabilityStructure , runif(1, 0.01 , 0.2)) )
EmulatorX <- rbind(EmulatorX , cbind(beta +0*EmpericalProbabilityStructure$x , EmpericalProbabilityStructure$x))
EmulatorY <- as.matrix(rbind(EmulatorY , as.matrix(EmpericalProbabilityStructure$y - EmpericalProbabilityStructure$x) ))
EmulatorSigma_n <- as.matrix(rbind(EmulatorSigma_n ,as.matrix(EmpericalProbabilityStructure$sd^2)))
DP_WaitBar(which(setofbeta == beta)/length(setofbeta))
}
EmulatorX <- EmulatorX[-1,]
EmulatorY <- as.matrix(EmulatorY[-1,])
EmulatorSigma_n <-EmulatorSigma_n[-1,]
# Fit emulator
BC_PlotPairs(cbind(EmulatorX , EmulatorY) , alpha=0.1)}
plot_ly(data = data.frame(alpha = EmulatorX[ , 1],
                          predicted_probability = EmulatorX[ , 2],
                          epsilon = EmulatorY), x = ~alpha, y = ~predicted_probability, z = ~epsilon, type = 'scatter3d' )

TrainingSamples <- sample(1:dim(EmulatorX)[1] , 1000)
LocalEmulationClass = BE_CreateDefaultEmulationClass()
LocalEmulationClass$X = EmulatorX[TrainingSamples,]
LocalEmulationClass$Y = EmulatorY[TrainingSamples]
LocalEmulationClass$MeanFunction = function(X){
  X <- as.matrix(X)
  H = cbind( as.matrix(1 + 0*X[,1]) , X  )
  return(H)
}
LocalEmulationClass$CorrelationLength <- function(X , n){
  return(c(0.1 , 0.1))
}
LocalEmulationClass$w <- function(X){
  return( diag(EmulatorSigma_n[TrainingSamples]) )
}
xstar = EmulatorX[-TrainingSamples,]

LocalEmulationClass <- BE_PerformPreCalulationForLSEEmulator( LocalEmulationClass )


emulatoroutput <- BE_BayesLinearEmulatorLSEstimates(xstar = xstar , EmulatorSettings = LocalEmulationClass)

BE_PlotStdResiduals( zstar = EmulatorY[-TrainingSamples] , emulatoroutput , e= EmulatorSigma_n[-TrainingSamples])

LocalEmulationClass <- BE_PerformPreCalulationForLSEEmulator(LocalEmulationClass)

i <- 2
patientname = DataBaseMaster$AFPatinetsNames[[i]]
if(DP_CheckDistributionSummariesExists(path , patientname) == FALSE){next}
MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017  , patientname)
RRData <- DP_LoadRpeaksfile(path , patientname) 
t <- RRData$RRCombined$t
# Some logic to avoid dodgy data

DistributionSummaries <- DP_LoadDistributionSummaries(path , patientname)
AdjustedBeliefs <- BC_EstimateGlobalandLocalParametersDisSum( DistributionSummaries )
Probabilities <- BC_CalulateDenistiesGMM( W = DP_RemoveNaRows(AdjustedBeliefs$W) , LocalDistributionStruct = LocalDistributionStruct )[ , 2:3]

Prevision <- matrix(0 , dim(Probabilities)[1]+1 , 2)
Prevision[1 , 1] =  DataSetPriorProbabilities$B
Prevision[1 , 2] = 0.01
referenceprobabilties <- matrix(0 , dim(Probabilities)[1]+1 , 1)
referenceprobabilties[1 , 1] =  DataSetPriorProbabilities$B

for(ii in 2:dim(Probabilities)[1]){
  
  UpdatedPrevision <- BC_ReifiedBeliefUpdateWithUncertainDensities(f_i = Probabilities[ii - 1, 1:2] , E_d = Prevision[ii -1 , 1] , V_d = Prevision[ii -1 , 2] ,LocalEmulationClass =  LocalEmulationClass , numbersamples = 25 )
  
  Prevision[ii , 1] <-  0.995*Prevision[ii-1 , 1]  +  0.005*UpdatedPrevision$E_Pt
  Prevision[ii , 2] <-  0.995*Prevision[ii-1 , 2]  +  0.005*UpdatedPrevision$V_Pt
  referenceprobabilties[ii,1] <- ((referenceprobabilties[ii-1,1])*Probabilities[ii - 1, 2]) / ((referenceprobabilties[ii-1,1])*Probabilities[ii - 1, 2] + (1-referenceprobabilties[ii-1,1])*Probabilities[ii - 1, 1] )
  if(mod(ii , 100) == 0){DP_WaitBar(ii/dim(Probabilities)[1])}

}


p1 <- BC_PlotCreateRRTimesPlots( RRData , MetaData = MetaData )
p2 <- BC_PlotPrevision( Prevision ) 

x11()
grid.arrange( p1 , p2 , nrow = 2 , ncol = 1)
