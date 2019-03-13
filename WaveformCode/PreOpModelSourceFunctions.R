POM_CreateDefaultSetofCovForPatIndex <- function(){
  return( c( 'PseudoId',
             'NewPseudoId',
             'Age' ,
             'Gender' ,
             'Weight' , 
             'ProcDetails' ,
             'CPB',
             'ConfirmedFirstNewAF',
             'FirstNewAF',
             'Pre_OperativeHeartRhythm'))
}
POM_CreateDefaultSetofCovForDet <- function(){
  return(c("LogisticEUROScore" ,
           "AdditiveEUROScore" ,
           "SCTSLogisticEuroSCORE"))
}
POM_CreateDefaultSetofCov <- function(){
  listofcovariates <- c('Age'  , 'CPB' ,  'AdditiveEUROScore' , 'PreOpNa', 'PreopUrea' ,'PreOpCRP', 'PreOpAlb' , 'PreopBili' , 'PreopCreat')
  return(listofcovariates)
}
POM_CreateDataStructure <- function(PatIndex2017 , DetIndex2017 , BioChemIndex2017, SetofCovariatesPatientIndex = POM_CreateDefaultSetofCovForPatIndex(), SetofCovariatesDetIndex = POM_CreateDefaultSetofCovForDet()){
  
  #Patient index data
  ReducedPatientIndex <-  PatIndex2017[ , which(names(PatIndex2017) %in% SetofCovariatesPatientIndex)] 
  ReducedPatientIndex <-  ReducedPatientIndex[ReducedPatientIndex$NewPseudoId %in% DetIndex2017$NewPseudoId , ]
  ReducedPatientIndex <-  ReducedPatientIndex[order(ReducedPatientIndex$NewPseudoId) , ]
  ReducedPatientIndex$CPB[is.na(ReducedPatientIndex$CPB)] <- "00:00"
  ReducedPatientIndex$CPB <- as.numeric(difftime(strptime(x = ReducedPatientIndex$CPB ,  format = '%R') ,  strptime(x = ReducedPatientIndex$CPB ,  format = '%R')[2],units = c('mins') ))
  
  ReducedDetIndex <- DetIndex2017[order(DetIndex2017$NewPseudoId) , which(names(DetIndex2017) %in% SetofCovariatesDetIndex)]
  
  #dendtite data
  MasterPreOpData <- cbind(ReducedPatientIndex , ReducedDetIndex)
  
  #Bio chemistry data
  PreOpBioChemData <- matrix(0 , length(BioChemIndex2017) , 8)
  for(i in 1:length(BioChemIndex2017)){
    PreOpBioChemData[i,1:8] <- as.numeric(BioChemIndex2017[[i]]$MetaData[ ,  6:13])
  }
  
  colnames(PreOpBioChemData) <- c(names(BioChemIndex2017[[i]]$MetaData[ ,  6:13]) )
  rownames(PreOpBioChemData) <- c(names(BioChemIndex2017) )
  PreOpBioChemData <- data.frame( NewPseudoId = rownames(PreOpBioChemData) , PreOpBioChemData)
  PreOpBioChemData$NewPseudoId <- as.character(PreOpBioChemData$NewPseudoId)
  
  MasterPreOpData <- MasterPreOpData[MasterPreOpData$NewPseudoId  %in% PreOpBioChemData$NewPseudoId , ]
  
  ReducedPreOpBioChemData <-  PreOpBioChemData[ PreOpBioChemData$NewPseudoId  %in%  MasterPreOpData$NewPseudoId , ] 
  ReducedPreOpBioChemData <-  ReducedPreOpBioChemData[order(ReducedPreOpBioChemData$NewPseudoId) , 2:9] 
  
  MasterPreOpData <- cbind(MasterPreOpData , ReducedPreOpBioChemData)
  MasterPreOpData$Gender <- as.factor(MasterPreOpData$Gender)
  
  for(ii in 1:dim(MasterPreOpData)[2]){
    if(is.numeric(MasterPreOpData[, ii]) ){
      MasterPreOpData[is.na(MasterPreOpData[, ii]), ii] <- mean(MasterPreOpData[!is.na(MasterPreOpData[, ii]), ii])
    }
  }
  
  MasterPreOpData <- MasterPreOpData[duplicated(MasterPreOpData$PseudoId) == 0 , ]
  
  MasterPreOpData$ConfirmedFirstNewAF[is.na(MasterPreOpData$ConfirmedFirstNewAF)] <- 'NA'
  MasterPreOpData <- data.frame(AFLogical = ((!is.na(MasterPreOpData$FirstNewAF))*(MasterPreOpData$ConfirmedFirstNewAF != 'CNAF')) == 1 , MasterPreOpData)
  
  
  return(MasterPreOpData)
}
POM_CreateBLBModellingStructure <- function( MasterPreOpData  , listofcovariates = POM_CreateDefaultSetofCov() ){
  indexesofcovariates <-  which(names(MasterPreOpData) %in%listofcovariates )
  
  AFLogical <- MasterPreOpData$AFLogical
  Data <- MasterPreOpData[ , indexesofcovariates]
  
  AFData <- Data[AFLogical == T, ]
  NAFData <- Data[AFLogical == F, ]
  
  mAF <- apply(AFData , 2 , mean) 
  vAF <- cov(AFData)
  mNAF <- apply(NAFData , 2 , mean) 
  vNAF <- cov(NAFData)
  
  ImAF  <- sqrt(apply(Data,   1 , function(X){(X - mAF)%*%(solve(vAF)%*%(X - mAF))}))
  ImNAF <- sqrt(apply(Data,   1 , function(X){(X - mNAF)%*%(solve(vNAF)%*%(X - mNAF))}))
  
  AFModel <- kde(cbind(ImAF[AFLogical == T ] , ImNAF[AFLogical == T ]) )
  NAFModel <- kde(cbind(ImAF[AFLogical == F ] , ImNAF[AFLogical == F ])) 
  
  output <- list(1,1,1,1,1,1)
  output[[1]] <- mAF
  output[[2]] <- vAF
  output[[3]] <- mNAF
  output[[4]] <- vNAF
  output[[5]] <- AFModel
  output[[6]] <- NAFModel
  
  return(setNames(output , c('mAF' , 'vAF' , 'mNAF' , 'vNAF' , 'AFModel' , 'NAFModel')))
  
}
POM_CalculateBLBProbability <- function(MasterPreOpData ,BLBModelStruct, listofcovariates = POM_CreateDefaultSetofCov() , patientID , PriorProbability = sum(MasterPreOpData$AFLogical)/length(MasterPreOpData$AFLogical)){
  indexesofcovariates <-  which(names(MasterPreOpData) %in%listofcovariates )
  
  Data <- MasterPreOpData[ MasterPreOpData$PseudoId %in% patientID, indexesofcovariates]
  
  ImAF  <- sqrt(apply(Data,   1 , function(X){(X - BLBModelStruct$mAF)%*%(solve(BLBModelStruct$vAF)%*%(X - BLBModelStruct$mAF))}))
  ImNAF <- sqrt(apply(Data,   1 , function(X){(X - BLBModelStruct$mNAF)%*%(solve(BLBModelStruct$vNAF)%*%(X - BLBModelStruct$mNAF))}))
  
  return(BLBF_CalculatePosteriorProbabilities(AFModel =BLBModelStruct$AFModel, NAFModel =BLBModelStruct$NAFModel , ImMatrix = cbind(ImAF , ImNAF) , alpha = PriorProbability ))
  
}
POM_CreateBLBCalibrationStructure <- function(PosteriorProbability , AFLogical , BinWidth  = 0.1){
  
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , AFLogical ), BinWidth = BinWidth ))
  
  CovariateEmulationClass = BE_CreateDefaultEmulationClass()
  CovariateEmulationClass$X = as.matrix(EmpericalProbabilityStructure$x)
  CovariateEmulationClass$Y = as.matrix(EmpericalProbabilityStructure$y - EmpericalProbabilityStructure$x)
  CovariateEmulationClass$MeanFunction = function(X){
    X <- as.matrix(X)
    H = cbind( matrix(1 ,dim(X)[1] , 1) , X  )
    return(H)
  }
  CovariateEmulationClass$CorrelationLength <- function(X , n){
    return(c(BinWidth ))
  }
  CovariateEmulationClass$w <- function(X){
    return( 0.1*diag(EmpericalProbabilityStructure$sd^2) )
  }
  CovariateEmulationClass <- BE_PerformPreCalulationForLSEEmulator( CovariateEmulationClass )
  
  xstar = seq(0,1,0.001)
  EmOutput <- BE_BayesLinearEmulatorLSEstimates(EmulatorSettings =  CovariateEmulationClass , xstar = xstar )
  
  output <- list(1 , 1, 1)
  output[[1]] <- xstar
  output[[2]] <- EmOutput$E_D_fX
  output[[3]] <- diag(EmOutput$V_D_fX)
  output <- setNames(output , c('Xstar' , 'E_D_fX' , 'V_D_fX'))
  return(output)
}
POM_CalculateCalibratedProbability <- function(PosteriorProbability  , BLBCalibrationStruct){
  
  PosteriorProbability <- as.matrix(PosteriorProbability)
  
  minindexes <- apply(PosteriorProbability , 1 , function(X){which.min( abs(X - BLBCalibrationStruct$Xstar)) } )
  
  E_D_P <- PosteriorProbability + BLBCalibrationStruct$E_D_fX[minindexes]
  V_D_P <-  BLBCalibrationStruct$V_D_fX[minindexes]
  
  return( setNames(list( as.matrix(E_D_P) , as.matrix(V_D_P) )  , c('E_D_P' , 'V_D_P')) )
  
}
POM_CalculatePreOpProbability <- function(MasterPreOpData , BLBModelStruct,BLBCalibrationStruct,listofcovariates = POM_CreateDefaultSetofCov(),patientID,PriorProbability = sum(MasterPreOpData$AFLogical)/length(MasterPreOpData$AFLogical)){
  PosteriorProbability <- POM_CalculateBLBProbability(MasterPreOpData = MasterPreOpData , BLBModelStruct = BLBModelStruct, listofcovariates = listofcovariates , patientID = patientID  , PriorProbability=PriorProbability)
  outut <- POM_CalculateCalibratedProbability(PosteriorProbability = PosteriorProbability , BLBCalibrationStruct  = BLBCalibrationStruct )
  outut$E_D_P[outut$E_D_P < 0] <- 0.001
  outut$E_D_P[outut$E_D_P > 1] <- 0.999
  
  return(outut)
}
