
DataSetPriorProbabilities <- BC_ExtractValidationPriors(Priorprobabilities , DataBaseMaster , DataBase)

GlobalUpdateDiagnostics <- matrix(0 , dim(DataBase[[3]])[1] + dim(DataBase[[4]])[1] , 1 )
counter <-1
rm(DataBase)
rm(DataBaseMaster)

# Load or create master data.
if( BCOptions[[1]] == 'DistributionSummaries' ){
    print('Pre computed database master not found. Processing.')
    DataBaseMaster <- BC_LoadInAllDistributionSummaries(path , listAllPatients , PatIndex2017)
}

if( BCOptions[[1]] == 'CDFs' ){
    print('Pre computed database master not found. Processing.')
    DataBaseMaster <- BC_LoadInAllCDFS(path , listAllPatients , PatIndex2017)
}


# Load or create database.
if(BCOptions[[2]] == 'AFClassifier' ){
    print('Pre computed database master not found. Processing.')
    DataBase <- BC_CreateAFandNOAFDataStructure(DataBaseMaster =  DataBaseMaster,PatIndex2017 = PatIndex2017)
  #for(i in 1:length(DataBase)){DataBase[[i]] <- DataBase[[i]][ , 2:31] }
  #source('BC_ValidationPlots.R')  
}

if(BCOptions$DensityEstimationGlobal == 'MVN'){
  for(i in 1:dim(DataBase[[3]])[1]){
    GlobalUpdateDiagnostics[counter,] <- BC_GlobalBayesianBeliefUpdateMVN(M = DataBase[[3]][i,] , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = DataSetPriorProbabilities)$A
    counter <- counter +1
  }
  for(i in 1:dim(DataBase[[4]])[1]){
    GlobalUpdateDiagnostics[counter,] <- BC_GlobalBayesianBeliefUpdateMVN(M = DataBase[[4]][i,] , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = DataSetPriorProbabilities)$A
    counter <- counter +1
  }
}
if(BCOptions$DensityEstimationGlobal == 'GMM'){
  for(i in 1:dim(DataBase[[3]])[1]){
    if(sum(is.na(t(as.matrix(DataBase[[3]][i,]))))>0){next}
    GlobalUpdateDiagnostics[counter,] <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix(DataBase[[3]][i,])) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = DataSetPriorProbabilities)$A
    counter <- counter +1
  }
  for(i in 1:dim(DataBase[[4]])[1]){
    if(sum(is.na(t(as.matrix(DataBase[[4]][i,]))))>0){next}
    GlobalUpdateDiagnostics[counter,] <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix(DataBase[[4]][i,])) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = DataSetPriorProbabilities)$A
    counter <- counter +1
  }
  GlobalUpdateDiagnostics[1:dim(CrossValidatedProbability)[1] , ] <- CrossValidatedProbability
}  


GlobalProbCalibrationStruct <- cbind(GlobalUpdateDiagnostics , 0*GlobalUpdateDiagnostics)
GlobalProbCalibrationStruct[1:dim(DataBaseMaster$AFPatientsDatabase)[1] , 2] <- 1 
GlobalProbCalibrationStruct <- DP_RemoveNaRows(GlobalProbCalibrationStruct)


#plot(GlobalProbCalibrationStruct[,1] , GlobalProbCalibrationStruct[,2]  , pch = 16 , col = rgb(0,0,1,alpha = 0.1))
GlobalProbCalibrationStruct <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(GlobalProbCalibrationStruct , BinWidth = 0.1))

x11(20,14)
print(BC_PlotCreateProbabilityCalibrationPlot(GlobalProbCalibrationStruct))

#bound = 0.5
#sensitivity <- sum( GlobalUpdateDiagnostics[ c(1:dim(DataBaseMaster$AFPatientsDatabase)[1]), ] > bound )/dim(DataBaseMaster$AFPatientsDatabase)[1]
#specificity <- sum( GlobalUpdateDiagnostics[ -c(1:dim(DataBaseMaster$AFPatientsDatabase)[1]), ] < bound)/dim(DataBaseMaster$NAFPatientsDatabase)[1]
#PPV <- ( 0.2*sensitivity )/( (0.2*sensitivity) + ( (1 - specificity )*0.8 ))
#NPV <- ( 0.8*specificity )/( (0.2*(1-sensitivity)) + ( ( specificity )*0.8 ))


GlobalEmulationClass = BE_CreateDefaultEmulationClass()
GlobalEmulationClass$X = GlobalProbCalibrationStruct$x
GlobalEmulationClass$CorrelationLength <- function(X , n){
  return(0.2)
}
GlobalEmulationClass$Y = GlobalProbCalibrationStruct$y - GlobalProbCalibrationStruct$x
GlobalEmulationClass$w <- function(X){
  return( diag(GlobalProbCalibrationStruct$sd^2) )
}
xstar = as.matrix(seq(0,1,0.001))

emulatoroutput <- BE_BayesLinearEmulatorLSEstimates(xstar =xstar , EmulatorSettings = GlobalEmulationClass)

BE_PlotOneDOutput(emulatoroutput , GlobalEmulationClass , xstar)

tmp <- BE_BayesLinearEmulatorLSEstimates(xstar = as.matrix(GlobalUpdateDiagnostics) , EmulatorSettings = GlobalEmulationClass)

E_Pt <- GlobalUpdateDiagnostics + tmp$E_D_fX
V_Pt <- as.matrix(diag(tmp$V_D_fX))

GlobalProbCalibrationStruct <- cbind( E_Pt , 0*GlobalUpdateDiagnostics)
GlobalProbCalibrationStruct[1:dim(DataBaseMaster$AFPatientsDatabase)[1] , 2] <- 1 
GlobalProbCalibrationStruct <- DP_RemoveNaRows(GlobalProbCalibrationStruct)


#plot(GlobalProbCalibrationStruct[,1] , GlobalProbCalibrationStruct[,2]  , pch = 16 , col = rgb(0,0,1,alpha = 0.1))
GlobalProbCalibrationStruct <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(GlobalProbCalibrationStruct , BinWidth = 0.05))

x11(20,14)
print(BC_PlotCreateProbabilityCalibrationPlot(GlobalProbCalibrationStruct))

