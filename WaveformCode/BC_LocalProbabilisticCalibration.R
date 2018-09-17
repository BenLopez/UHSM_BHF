
DataSetPriorProbabilities <- BC_ExtractValidationPriors(Priorprobabilities , DataBaseMaster , DataBase)

LocalUpdateDiagnostics <- list()
counter <-1

for(i in 1:dim(DataBaseMaster$AFPatientsDatabase)[1]){
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017  , DataBaseMaster$AFPatinetsNames[[i]])
  t <- DP_LoadRpeaksfile(path ,  DataBaseMaster$AFPatinetsNames[[i]])$RRCombined$t
  t <- t[ t %in% DataBaseMaster$AFPatientsDatabase[i , , 12]]
  if(length(t)< 5000){
    next
  }
  Z <- DataBaseMaster$AFPatientsDatabase[i , , 1:11]
  Z[is.na(Z)] <- 100
  Z[Z[,8] == 1 , 9] <- 0
  Z[Z[,8] == 1,10] <- 0
  AdjustedBeliefs <- BC_EstimateGlobalandLocalParameters(Z)
  
  if(BCOptions$GlobalUpdate == 'Yes'){
    if(BCOptions$DensityEstimationGlobal == 'MVN'){
      GlobalUpdatedBeliefs <- BC_GlobalBayesianBeliefUpdateMVN(M = AdjustedBeliefs$M , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = DataSetPriorProbabilities)
    }
    if( BCOptions$DensityEstimationGlobal == 'GMM' ){
      if( sum(is.na(t(as.matrix( AdjustedBeliefs$M))))>0 ){
        next}
      GlobalUpdatedBeliefs <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix( AdjustedBeliefs$M)) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = DataSetPriorProbabilities)
    }
  }else{
    GlobalUpdatedBeliefs <- DataSetPriorProbabilities 
  }
  
  LocalUpdateDiagnostics[[counter]] <- matrix(0 , length(t) , 2 )
  Probabilities <- BC_BayesianBeliefUpdateGMM( W = AdjustedBeliefs$W , LocalDistributionStruct =  LocalDistributionStruct , Probabilities = GlobalUpdatedBeliefs , n = BCParameters$TS_Likelihood_clique )[1:length(t),1]
  Implausability <- BC_CalulateCulmulativeImplausability( Implausability = BC_TestingCalulateRegularisedImplausabilty(AdjustedBeliefs$W , SecondOrderStruct = LocalSecondOrderStruct , ImSecondOrderStruct) , n = BCParameters$TS_Likelihood_clique )[1:length(t),]
  Probabilities[Implausability[,1]>3] <- 0.00000000000001
  LocalUpdateDiagnostics[[counter]][1:length(t) , 1] <-  Probabilities
  LocalUpdateDiagnostics[[counter]][1:length(t) , 2] <-  BC_CreateAFAnnotationFomMetaData(t , MetaData)
  Implausability[is.na(Implausability)] <- 100
  if(sum(((Implausability[,1]>3)*(Implausability[,2]>3)) == 1) >0 ){
  LocalUpdateDiagnostics[[counter]] <- LocalUpdateDiagnostics[[counter]][-which(((Implausability[,1]>3)*(Implausability[,2]>3)) == 1) , ]
  }
  DP_WaitBar(counter / (dim(DataBaseMaster$AFPatientsDatabase)[1] + dim(DataBaseMaster$NAFPatientsDatabase)[1]))  
  counter <- counter +1
}

for(i in 1:dim(DataBaseMaster$NAFPatientsDatabase)[1]){
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017  , DataBaseMaster$NAFPatinetsNames[[i]])
  if(nrow(MetaData) == 0){next}
  t <- DP_LoadRpeaksfile(path ,  DataBaseMaster$NAFPatinetsNames[[i]])$RRCombined$t
  t <- t[ t %in% DataBaseMaster$NAFPatientsDatabase[i , , 12]]
  if(length(t)< 5000){
    next
  }
  Z <- DataBaseMaster$NAFPatientsDatabase[i , , 1:11]
  Z[is.na(Z)] <- 100
  Z[Z[,8] == 1 , 9] <- 0
  Z[Z[,8] == 1,10] <- 0
  AdjustedBeliefs <- BC_EstimateGlobalandLocalParameters(Z)
  
  if(BCOptions$GlobalUpdate == 'Yes'){
  if(BCOptions$DensityEstimationGlobal == 'MVN'){
  GlobalUpdatedBeliefs <- BC_GlobalBayesianBeliefUpdateMVN(M = AdjustedBeliefs$M , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = DataSetPriorProbabilities)
  }
  if( BCOptions$DensityEstimationGlobal == 'GMM' ){
  if( sum(is.na(t(as.matrix( AdjustedBeliefs$M))))>0 ){
      next}
    GlobalUpdatedBeliefs <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix( AdjustedBeliefs$M)) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = DataSetPriorProbabilities)
  }
  }else{
    GlobalUpdatedBeliefs <- DataSetPriorProbabilities 
  }
  
  LocalUpdateDiagnostics[[counter]] <- matrix(0 , length(t) , 2 )
  Probabilities <- BC_BayesianBeliefUpdateGMM( W = AdjustedBeliefs$W , LocalDistributionStruct =  LocalDistributionStruct , Probabilities = GlobalUpdatedBeliefs , n = BCParameters$TS_Likelihood_clique )[1:length(t),1]
  Implausability <- BC_CalulateCulmulativeImplausability( Implausability = BC_TestingCalulateRegularisedImplausabilty(AdjustedBeliefs$W , SecondOrderStruct = LocalSecondOrderStruct , ImSecondOrderStruct) , n = BCParameters$TS_Likelihood_clique )[1:length(t),]
  Probabilities[Implausability[,1]>3] <- 0.00000000000001
  LocalUpdateDiagnostics[[counter]][1:length(t) , 1] <-  Probabilities
  LocalUpdateDiagnostics[[counter]][1:length(t) , 2] <-  0
  Implausability[is.na(Implausability)] <- 100
  if(sum(((Implausability[,1]>3)*(Implausability[,2]>3)) == 1) >0 ){
    LocalUpdateDiagnostics[[counter]] <- LocalUpdateDiagnostics[[counter]][-which(((Implausability[,1]>3)*(Implausability[,2]>3)) == 1) , ]
  }
  DP_WaitBar(counter / (dim(DataBaseMaster$AFPatientsDatabase)[1] + dim(DataBaseMaster$NAFPatientsDatabase)[1]))  
  counter <- counter +1
}

LocalProbCalibrationStruct <- matrix(0 , length(unlist(LocalUpdateDiagnostics))/2 , 2 )
StartIndex <- 1
for(i in 1:length(LocalUpdateDiagnostics )){
  LocalProbCalibrationStruct[StartIndex:(StartIndex + dim(LocalUpdateDiagnostics[[i]])[1] - 1) , ] <- LocalUpdateDiagnostics[[i]]
  StartIndex <- (StartIndex + dim(LocalUpdateDiagnostics[[i]])[1])
}


numberofsamples = 1000000
ProbabiliticCalibrationOutput <- BC_CreateCalibrationStructure(DP_RemoveNaRows(LocalProbCalibrationStruct[sample(1:dim(LocalProbCalibrationStruct)[1] , numberofsamples),]) , BinWidth = 0.05)
BC_PlotCreateProbabilityCalibrationPlot(ProbabiliticCalibrationOutput) + ggtitle('Local Calibration Probabilities')

#EmulatorOurput <- BE_BayesLinearEmulatorGLSEstimates(y = as.matrix(ProbabiliticCalibrationOutput$x)  , 
#                                  x =  as.matrix(ProbabiliticCalibrationOutput$y) , 
#                                 xstar = as.matrix(ProbabiliticCalibrationOutput$y) , 
#                                w = diag(ProbabiliticCalibrationOutput$sd^2)/0.1,
#                               l = 0.075 , p = 2 , h = function(X){cbind(matrix(1 , dim(X)[1]) , as.matrix(X)  )})



#plot(  as.matrix(ProbabiliticCalibrationOutput$x) , EmulatorOurput$E_D_fX )
#abline(0,1)
#EmulatorOurput <- BE_BayesLinearEmulatorGLSEstimates(y = as.matrix(ProbabiliticCalibrationOutput$x)  , 
#                                                    x =  as.matrix(ProbabiliticCalibrationOutput$y) , 
#                                                   xstar = as.matrix(seq(0,1,0.001)) , 
#                                                  w = diag(ProbabiliticCalibrationOutput$sd^2)/0.1,
#                                                 l = 0.075 , p = 2 , h = function(X){cbind(matrix(1 , dim(X)[1]) , as.matrix(X)  )})

# lines(EmulatorOurput$E_D_fX , as.matrix(seq(0,1,0.001)) )
# points( as.matrix(ProbabiliticCalibrationOutput$x) , as.matrix(ProbabiliticCalibrationOutput$y) , col ='blue')
