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


DataStructure <- array(0 , c(length(listAllPatients) , 13 , 400))


for(i in 1:length(listAllPatients) ){
  PatientID <- listAllPatients[[i]]
  if(DP_CheckFileExists(path , PatientID , Name = paste0(PatientID , '_RRDistributionSummariesForeCasting'))){
    tmp <-  DP_LoadFile(path , PatientID , Name = paste0(PatientID , '_RRDistributionSummariesForeCasting'))
    tmp <- as.matrix(tmp[,tmp[1,]!=0  ])
    DataStructure[i , , (dim(DataStructure)[3] - dim(tmp)[2]+1):dim(DataStructure)[3] ] <- tmp
  }
}

# Remove rows with all zero rows
PatientsInStudy <- listAllPatients[-which(apply(DataStructure[, 1 , ] , 1 , function(X){sum(X!=0)} ) < 30 )]
AFLogical <- apply(as.matrix(PatientsInStudy) , 1 , function(X){DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017  , X))})
alpha = c(sum(AFLogical)/length(AFLogical) , (1-(sum(AFLogical)/length(AFLogical))))
DataStructure <- DataStructure[-which(apply(DataStructure[, 1 , ] , 1 , function(X){sum(X!=0)} ) < 30 ) ,  , ]
rownames(DataStructure) <- PatientsInStudy

DataStructure[ , 3 , ]   <-  log(DataStructure[, 3, ])
DataStructure[ , 4 , ]  <-  sqrt(DataStructure[ , 4, ])
DataStructure[ , 7 , ]  <-  log(DataStructure[ , 7, ])
DataStructure[ , 10 , ] <-  sqrt(DataStructure[ , 10, ])
DataStructure[ , 11 , ] <-  sqrt(DataStructure[ , 11, ])
DataStructure[ , 12 , ] <-  sqrt(DataStructure[ , 12, ])


# Train Pre Op model
{
  MasterPreOpData <- POM_CreateDataStructure(PatIndex2017 , DetIndex2017 , BioChemIndex2017)
  ReducedMasterPreOpData <- MasterPreOpData[MasterPreOpData$Pre_OperativeHeartRhythm != 'Atrial fibrillation/flutter' , ]
  #ReducedMasterPreOpData <- ReducedMasterPreOpData[ReducedMasterPreOpData$PseudoId %in% PatientsInStudy , ]
  BLBModelStruct <- POM_CreateBLBModellingStructure(MasterPreOpData =  ReducedMasterPreOpData  , listofcovariates = POM_CreateDefaultSetofCov() )
  PosteriorProbability <- POM_CalculateBLBProbability(MasterPreOpData = ReducedMasterPreOpData , BLBModelStruct, listofcovariates = POM_CreateDefaultSetofCov() , ReducedMasterPreOpData$PseudoId )
  BLBCalibrationStruct <- POM_CreateBLBCalibrationStructure(PosteriorProbability , ReducedMasterPreOpData$AFLogical )
  PreOpProbs <- POM_CalculatePreOpProbability(MasterPreOpData = ReducedMasterPreOpData , BLBModelStruct = BLBModelStruct , BLBCalibrationStruct = BLBCalibrationStruct , patientID = ReducedMasterPreOpData$PseudoId)
  PreOpProbs <- POM_CalculatePreOpProbability(MasterPreOpData = ReducedMasterPreOpData , BLBModelStruct = BLBModelStruct , BLBCalibrationStruct = BLBCalibrationStruct , patientID = PatientsInStudy , PriorProbability = sum(AFLogical)/length(AFLogical))
  PreOpProbs[is.na(PreOpProbs)] <- sum(AFLogical)/length(AFLogical)
  PreOpProbs$E_D_P <-POM_CalculateBLBProbability(MasterPreOpData = ReducedMasterPreOpData , BLBModelStruct, listofcovariates = POM_CreateDefaultSetofCov() , PatientsInStudy  , PriorProbability = sum(AFLogical)/length(AFLogical))
}


#indextoview <- ((dim(DataStructure)[3] - 29) : (dim(DataStructure)[3]))
indextoview <- (dim(DataStructure)[3] - 29  ) 
{
  Im <- BLBL_CrossValidateImplausabilitiesALL(TrainingData = DataStructure[, -c( 13), indextoview] , AFLogical)
  ImAF <- Im[,1]
  ImNAF <- Im[,2]
 
  ImAF <- (ImAF - mean(ImAF[AFLogical] , na.rm = T))/sqrt(var(ImAF[AFLogical]))
  ImNAF <- (ImNAF - mean(ImNAF[AFLogical ==0] , na.rm = T))/sqrt(var(ImNAF[AFLogical==0], na.rm = T))
  
   
  x11()
  par(mfrow =c(1,2))
  print(BC_PlotCompareSingleHists((ImAF[AFLogical ==F]) , (ImAF[AFLogical ==T])  ,breaks = 15 ,  main = 'Implausabilty of Going into AF', xlab = 'Im'))
  print(BC_PlotCompareSingleHists((ImNAF[AFLogical ==F]) , (ImNAF[AFLogical ==T]) ,breaks = 15 , main = 'Implausabilty of Not Going into AF' , xlab = 'Im'))
  
  x11()
  plot(ImAF , ImNAF, col = 'blue' , pch = 16)
  points(ImAF[AFLogical == T] , ImNAF[AFLogical == T], col = 'red', pch = 16)
  abline(h = 3)
  abline(v = 3)
  
  PosteriorProbability <- BLBF_CrossValidateKDEall( ImAF , ImNAF , AFLogical  )
  PosteriorProbability[is.na(PosteriorProbability)] = c(sum(AFLogical)/length(AFLogical))

  x11(20,14)
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , AFLogical == T), BinWidth = 0.1))
  print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Forecasting Emperical Probabilities'))
  PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability , AFLogical == T))
  PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability , AFLogical == T))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('Forecasting ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Forecasting PPV vs NPV Curves')
  x11(20,14)
  grid.arrange(ROCplot , NPVPPVPlot , nrow= 2)
}


StartIndex <- 382
NumberofPoints <- 1
VariablesToView <- c(1:12)
test <- BLBF_DoubleCrossValidateInference(DataStructure ,AFLogical, StartIndex , NumberofPoints , VariablesToView  )
 

test[,1]= (test[,1] - mean(test[AFLogical,1]))/sqrt(var(test[AFLogical,1]))
test[,2]= (test[,2] - mean(test[AFLogical ==0,2]))/sqrt(var(test[AFLogical==0,2]))


x11(20,14)
par(mfrow = c(1 , 2))
print(BC_PlotCompareSingleHists((test[AFLogical ==F , 1]) , (test[AFLogical ==T , 1]) ))
print(BC_PlotCompareSingleHists((test[AFLogical ==F , 2]) , (test[AFLogical ==T , 2]) ))

AFModel <- kde( cbind( test[AFLogical ==T , 1] , test[AFLogical ==T , 2] )  )
NAFModel <- kde( cbind( test[AFLogical ==F , 1] , test[AFLogical ==F , 2] )  )

x11(20,14)
plot(NAFModel  , col ='blue' , xlab =c('Implausability AF') , ylab =c('Implausability NAF'))
plot(AFModel, col = 'red' , add = T)

points(test[AFLogical == F , 1] , test[AFLogical == F ,2] , col = 'blue' , pch = 16)
points(test[AFLogical == T , 1] , test[AFLogical == T , 2], col = 'red', pch = 16)
abline(3 , 0)
abline(v = 3)

PosteriorProbability2 <- 0*ImAF
PosteriorProbability2[((ImAF > 3)*(ImAF>3))==1] <- sum(AFLogical)/length(AFLogical)
PosteriorProbability2[((ImAF > 3)*(ImAF>3))==0] <- BLBF_CrossValidateKDEall( test[((ImAF > 3)*(ImAF>3))==0,1] , test[((ImAF > 3)*(ImAF>3))==0,2] , AFLogical[((ImAF > 3)*(ImAF>3))==0]  , PriorProbability = PosteriorProbability[((ImAF > 3)*(ImAF>3))==0 , ])


x11(20,14)
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability2 , AFLogical == T), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Forecasting Emperical Probabilities'))
PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability2 , AFLogical == T))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('Forecasting ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Forecasting PPV vs NPV Curves')

x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot , nrow= 2)

IndividualProbabilityUpdateStruct <- matrix(0 , length(AFLogical) , 30)
for(i in c(0:29)){
  IndividualProbabilityUpdateStruct[,i+1] <- BLBF_CalculatePosteriorProbabilityNBLA(DataStructure = DataStructure , AFLogical = AFLogical   , indextoview = 371 + i )
}

IndependentProbabilityUpdateStruct <- matrix(0 , length(AFLogical) , 30)
for(i in c(0:29)){
  if(i == 0){
    IndependentProbabilityUpdateStruct[,i+1] <- BLBF_CalculatePosteriorProbabilityNBLA(DataStructure = DataStructure , AFLogical = AFLogical  , indextoview = 371 + i )}
  if(i > 0){
    IndependentProbabilityUpdateStruct[,i+1] <- BLBF_CalculatePosteriorProbabilityNBLA(DataStructure = DataStructure , AFLogical = AFLogical  , indextoview = 371 + i , PriorProbability = IndependentProbabilityUpdateStruct[,i])}
}

BLAProbabilityUpdateStruct <- matrix(0 , length(AFLogical) , 30)
BLAProbabilityUpdateStruct[ , 1] <-  IndividualProbabilityUpdateStruct[ , 1]
for(i in c(1:29)){
  BLAProbabilityUpdateStruct[,i+1] <- BLBF_CalculatePosteriorProbabilityBLA(DataStructure = DataStructure , AFLogical = AFLogical  , StartIndex = max(371 , 371 + (i-2)) , NumberofPoints = min(2 , i) , PriorProbability = BLAProbabilityUpdateStruct[,i] , VariablesToView = c(1:12))
}

BLAProbabilityUpdateStruct <- matrix(0 , length(AFLogical) , 30)
BLAProbabilityUpdateStruct[ , 1] <-  IndividualProbabilityUpdateStruct[ , 1]
for(i in c(1)){
  BLAProbabilityUpdateStruct[,i+1] <- BLBF_CalculatePosteriorProbabilityBLA(DataStructure = DataStructure , AFLogical = AFLogical  , StartIndex = max(371 , 371 + (i-1)) , NumberofPoints =  1, PriorProbability = BLAProbabilityUpdateStruct[,i] , VariablesToView = c(1:12))
}


BLBF_PlotProbabilityTrajectories( IndividualProbabilityUpdateStruct  , AFLogical )
title('RR-Times Probability Trajectory')
BLBF_PlotSenSpecPerformanceSweep( IndividualProbabilityUpdateStruct , AFLogical )
title('ROC Curves Over Time')
BLBF_PlotProbabilityTrajectories( t(apply(IndividualProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical )
title('Mean RR-Times Probability Trajectory')
BLBF_PlotSenSpecPerformanceSweep( t(apply(IndividualProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical )
title('ROC Curves Over Time')

BLBF_PlotProbabilityTrajectories(IndependentProbabilityUpdateStruct  , AFLogical )
BLBF_PlotProbabilityTrajectories(t(apply(IndependentProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical )
BLBF_PlotSenSpecPerformanceSweep(IndependentProbabilityUpdateStruct , AFLogical )
BLBF_PlotSenSpecPerformanceSweep(t(apply(IndependentProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical)

BLBF_PlotProbabilityTrajectories(BLAProbabilityUpdateStruct  , AFLogical )
BLBF_PlotSenSpecPerformanceSweep(BLAProbabilityUpdateStruct , AFLogical )
BLBF_PlotProbabilityTrajectories(t(apply(BLAProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical )
BLBF_PlotSenSpecPerformanceSweep(t(apply(BLAProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical)


