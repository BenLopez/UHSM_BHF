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
  #listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

# Train model
{
  MasterPreOpData <- POM_CreateDataStructure(PatIndex2017 , DetIndex2017 , BioChemIndex2017)
  ReducedMasterPreOpData <- MasterPreOpData[MasterPreOpData$Pre_OperativeHeartRhythm != 'Atrial fibrillation/flutter' , ]
  BLBModelStruct <- POM_CreateBLBModellingStructure(MasterPreOpData =  ReducedMasterPreOpData  , listofcovariates = POM_CreateDefaultSetofCov() )
  PosteriorProbability <- POM_CalculateBLBProbability(MasterPreOpData = ReducedMasterPreOpData , BLBModelStruct, listofcovariates = POM_CreateDefaultSetofCov() , ReducedMasterPreOpData$PseudoId )
  BLBCalibrationStruct <- POM_CreateBLBCalibrationStructure(PosteriorProbability , ReducedMasterPreOpData$AFLogical )

  }



MasterPreOpData <- POM_CreateDataStructure(PatIndex2017 , DetIndex2017 , BioChemIndex2017)

# Remove patients with preoperative flutter
ReducedMasterPreOpData <- MasterPreOpData[MasterPreOpData$Pre_OperativeHeartRhythm != 'Atrial fibrillation/flutter' , ]
ReducedMasterPreOpData <- data.frame(AFLogical = !is.na(ReducedMasterPreOpData$FirstNewAF) , ReducedMasterPreOpData)

for(ii in 1:dim(ReducedMasterPreOpData)[2]){
if(is.numeric(ReducedMasterPreOpData[, ii]) ){
     ReducedMasterPreOpData[is.na(ReducedMasterPreOpData[, ii]), ii] <- mean(ReducedMasterPreOpData[!is.na(ReducedMasterPreOpData[, ii]), ii])
  }
}

summary(glm(formula = AFLogical ~ Age  + Gender + Weight + CPB + AdditiveEUROScore +LogisticEUROScore + SCTSLogisticEuroSCORE+ PreOpNa+PreOpK+PreopUrea+PreopCreat+PreOpCRP+PreOpAlb+PreopBili+PreopMg,family=binomial(link='logit') , data=ReducedMasterPreOpData))

model <- (glm(formula = AFLogical ~ Age  + CPB + AdditiveEUROScore  + PreOpNa+PreopUrea+PreOpCRP+PreOpAlb,family=binomial(link='logit') , data=ReducedMasterPreOpData))
summary(model)

AFLogical <- ReducedMasterPreOpData$AFLogical
BC_PlotCompareSingleHists(model$fitted.values[AFLogical == F] , model$fitted.values[AFLogical == T])

listofcovariates <- c('Age'  , 'CPB' ,  'AdditiveEUROScore' , 'PreOpNa', 'PreopUrea' ,'PreOpCRP', 'PreOpAlb' , 'PreopBili' , 'PreopCreat')
indexesofcovariates <-  which(names(ReducedMasterPreOpData) %in%listofcovariates )

DataMatrix <- ReducedMasterPreOpData[ , indexesofcovariates]

BC_PlotCompareTwoHists(DataMatrix[AFLogical == 1 ,] ,DataMatrix[AFLogical == 0 , ]  )
BC_PlotPairsFromTwoVariables(DataMatrix[AFLogical == 1 ,] ,DataMatrix[which(AFLogical == 0)[sample(1:sum(AFLogical == 0) , sum(AFLogical == 1))] ,]  , alpha = 0.1)

PosteriorProbability <- BLBC_FitBayesLinearBayesClassifier(Data =  DataMatrix , Labels = AFLogical )


{
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , Labels == 1), BinWidth = 0.1))

  CovariateEmulationClass = BE_CreateDefaultEmulationClass()
  CovariateEmulationClass$X = as.matrix(EmpericalProbabilityStructure$x)
  CovariateEmulationClass$Y = as.matrix(EmpericalProbabilityStructure$y - EmpericalProbabilityStructure$x)
  CovariateEmulationClass$MeanFunction = function(X){
    X <- as.matrix(X)
    H = cbind( as.matrix(1 + 0*X[,1]) , X  )
    return(H)
  }
  CovariateEmulationClass$CorrelationLength <- function(X , n){
    return(c(0.1 , 0.1))
  }
  CovariateEmulationClass$w <- function(X){
    return( 0.1*diag(EmpericalProbabilityStructure$sd^2) )
  }
  CovariateEmulationClass <- BE_PerformPreCalulationForLSEEmulator( CovariateEmulationClass )
  
  xstar = seq(0,1,0.01)
  EmOutput <- BE_BayesLinearEmulatorLSEstimates(EmulatorSettings =  CovariateEmulationClass , xstar = xstar )
  BE_PlotOneDOutput(EmOutput , CovariateEmulationClass , xstar)
  
  
  xstar = as.matrix(PosteriorProbability)
  EmOutput <- BE_BayesLinearEmulatorLSEstimates(EmulatorSettings =  CovariateEmulationClass , xstar = xstar )
  
  x11(20,14)
  E_D_Prevision <- as.matrix(PosteriorProbability + EmOutput$E_D_fX)
  E_D_Prevision[E_D_Prevision > 1] <- 0.99999999
  E_D_Prevision[E_D_Prevision < 0] <- 0.00000001
  
  CalibratedProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(E_D_Prevision , Labels == 1), BinWidth = 0.1))
  BC_PlotCreateProbabilityCalibrationPlot(CalibratedProbabilityStructure) + ggtitle('Calibrated Covariate Emperical Probabilities')
}

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(E_D_Prevision , Labels == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept =( 1-PriorProbability)) + geom_hline(yintercept = PriorProbability)+  ggtitle('PPV vs NPV Curves')
x11()


