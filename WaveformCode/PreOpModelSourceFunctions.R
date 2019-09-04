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
             'Pre_OperativeHeartRhythm',
             'PreopRRT',
             'Height'))
}
POM_CreateDefaultSetofCovForDet <- function(){
  return(c("LogisticEUROScore" ,
           "AdditiveEUROScore" ,
           "SCTSLogisticEuroSCORE",
           "EjectionFractionCategory",
           "NYHAGrade",
           "AnginaGrade",
           "Urgency"))
}
POM_CreateDefaultSetofCovForFLow <- function(){
  return(c("HR" ,
           "ReliableART.S",
           "ReliableART.M",
           "ReliableART.D",
           "CVP",
           "SpO2",
           "Lac",
           "FiO26num",
           "ArtPO2",
           "IABP2",
           "GCSnum"))
}
POM_CreateDefaultSetofCov <- function(){
  listofcovariates <- c('Age'  ,
                        'CPB' ,
                        'AdditiveEUROScore' ,
                        'PreOpNa',
                        'PreopUrea' ,
                        'PreOpCRP',
                        'PreOpAlb' ,
                        'PreopBili' ,
                        'PreopCreat' )
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
  
  #for(ii in 1:dim(MasterPreOpData)[2]){
  #  if(is.numeric(MasterPreOpData[, ii]) ){
  #    MasterPreOpData[is.na(MasterPreOpData[, ii]), ii] <- mean(MasterPreOpData[!is.na(MasterPreOpData[, ii]), ii])
  #  }
  #}
  
  MasterPreOpData <- MasterPreOpData[duplicated(MasterPreOpData$PseudoId) == 0 , ]
  
  MasterPreOpData$ConfirmedFirstNewAF[is.na(MasterPreOpData$ConfirmedFirstNewAF)] <- 'NA'
  
  MasterPreOpData <- data.frame(AFLogical = ((!is.na(MasterPreOpData$FirstNewAF))*(MasterPreOpData$ConfirmedFirstNewAF != 'CNAF')) , MasterPreOpData)
  
  
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
POM_ExtractPreOpFromHaem <- function(HaemIndex2017){
  HaemDataMatrix <- matrix(0 , length(HaemIndex2017) , 5)
  PatientNames <- matrix(0 , length(HaemIndex2017) , 1)
  counter <- 1
  for( i in 1:dim(HaemDataMatrix)[1] ){
    PatientNames[counter] = HaemIndex2017[[i]]$MetaData$NewPseudoId
    HaemDataMatrix[counter,] = as.numeric(HaemIndex2017[[i]]$MetaData[ , c(4,5,6,8,11 ) ])
    counter <- counter + 1
  }
  
  rownames(HaemDataMatrix) = PatientNames
  colnames(HaemDataMatrix) = names(HaemIndex2017[[1]]$MetaData[ , c(4,5,6,8,11 ) ])
  
  return(HaemDataMatrix)
}
#POM_ExtractAFMedication <- function(FluidsIndex2017){
#  DidoxinLogical <- matrix(0 , length(FluidsIndex2017) , 1)
#  AmiodaroneLogical <- matrix(0 , length(FluidsIndex2017) , 1)   
#  for(i in 1:length(FluidsIndex2017)){
#    
#    DidoxinLogical[i,] <-   sum(!is.na(FluidsIndex2017[[i]]$TimeSeriesData$Digoxin..)) > 0
#    AmiodaroneLogical[i,] <- sum(!is.na(FluidsIndex2017[[i]]$TimeSeriesData[ , which(grepl('Amiodarone', names(FluidsIndex2017[[i]]$TimeSeriesData) ))]))>0
#  }
#  outputstruct <- data.frame(NewPseudoId = names(FluidsIndex2017) , DidoxinLogical = DidoxinLogical , AmiodaroneLogical = AmiodaroneLogical)
#  return(outputstruct)
#}
POM_ExtractAFMedication <- function(AllDataStructure){
  DidoxinLogical <- matrix(0 , length(AllDataStructure) , 1)
  AmiodaroneLogical <- matrix(0 , length(AllDataStructure) , 1)   
  for(i in 1:length(AllDataStructure)){
    
    DidoxinLogical[i,] <-   sum(!is.na(AllDataStructure[[i]]$timeseriesvariables$Digoxin..)) > 0
    AmiodaroneLogical[i,] <- sum(!is.na(AllDataStructure[[i]]$timeseriesvariables[ , which(grepl('Amiodarone', names(AllDataStructure[[i]]$timeseriesvariables) ))]))>0
  }
  outputstruct <- data.frame(NewPseudoId = names(AllDataStructure) , DidoxinLogical = DidoxinLogical , AmiodaroneLogical = AmiodaroneLogical)
  return(outputstruct)
}
POM_ExtractPostOpFromBioChem <- function(BioChemIndex2017){
  
  PostOpBioChem <- matrix(0 , length(BioChemIndex2017) , 8)
  
  for(i in 1:length(BioChemIndex2017)){
    PostOpBioChem[i , 1:8] <- as.numeric(BioChemIndex2017[[i]]$TimeSeriesData[1 , 2:9])
  }
  
  colnames(PostOpBioChem) <- c(names(BioChemIndex2017[[i]]$TimeSeriesData[1 , 2:9]) )
  rownames(PostOpBioChem) <- names(BioChemIndex2017)
  return(PostOpBioChem)
}

POM_ExtractPostOpFromHaem <- function(HaemIndex2017){
  {
  PostopHaemDataMatrix <- matrix(0 , length(HaemIndex2017) , 7)
  for( i in 1:dim(PostopHaemDataMatrix)[1] ){
    if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[1 , 2])){
      PostopHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(5,6,8,9) ])
      PostopHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[1 , c(2,3 , 4) ])
    }
    if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[2 , 2])){
      PostopHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(2,3 , 4) ])
      PostopHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[1 , c(5,6,8,9) ])
    }
    if(is.na(HaemIndex2017[[i]]$TimeSeriesData[2 , 2]) & is.na(HaemIndex2017[[i]]$TimeSeriesData[1 , 2])){
      PostopHaemDataMatrix[i, ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(2,3 , 4 , 5,6,8,9) ])
    }
  }
  rownames(PostopHaemDataMatrix) = names(HaemIndex2017)
  colnames(PostopHaemDataMatrix) = names(HaemIndex2017[[1]]$TimeSeriesData[1 , c(2,3 , 4, 5,6,8,9) ])
  }
  return(PostopHaemDataMatrix)
}
POM_ExtractPostOpFromFlow<- function(FlowIndex2017 , PatIndex2017){
  PostopFlowDataMatrix <- matrix(0 , length(FlowIndex2017) , length(POM_CreateDefaultSetofCovForFLow()))
  
  for(i in 1:length(FlowIndex2017)){
    PatientRecord <- PatIndex2017[ PatIndex2017$NewPseudoId == names(FlowIndex2017)[i], ]
    dayvector <- abs(difftime(DP_StripTime(PatientRecord$FirstITUEntry) ,FlowIndex2017[[i]]$TimeSeriesData$time , units= 'hours')) < 24
      
    FlowIndex2017[[i]]$TimeSeriesData[ , dayvector]

    FlowIndex2017[[i]]$TimeSeriesData$CVP[FlowIndex2017[[i]]$TimeSeriesData$CVP == '---'] = NA
    FlowIndex2017[[i]]$TimeSeriesData$CVP[FlowIndex2017[[i]]$TimeSeriesData$CVP == ''] = NA
    FlowIndex2017[[i]]$TimeSeriesData$CVP <- as.numeric(FlowIndex2017[[i]]$TimeSeriesData$CVP)
    FlowIndex2017[[i]]$TimeSeriesData$IABP2 <- as.numeric(FlowIndex2017[[i]]$TimeSeriesData$IABP2 == 'Yes')
    PostopFlowDataMatrix[i ,] <- as.numeric(FlowIndex2017[[i]]$TimeSeriesData[1 , names(FlowIndex2017[[i]]$TimeSeriesData)%in%POM_CreateDefaultSetofCovForFLow() ])
  }
  rownames(PostopFlowDataMatrix) <- names(FlowIndex2017)
  colnames(PostopFlowDataMatrix) <- names(FlowIndex2017[[i]]$TimeSeriesData)[names(FlowIndex2017[[i]]$TimeSeriesData)%in%POM_CreateDefaultSetofCovForFLow()]
  return(PostopFlowDataMatrix)
}

POM_ExtractPostOpFromFluids<- function(FluidsIndex2017){
  PostopFluidsDataMatrix <- matrix(0 , length(FluidsIndex2017) , 2)
  for(i in 1:length(FluidsIndex2017)){
    
      diff <- FluidsIndex2017[[i]]$TimeSeriesData$Hourly.In - FluidsIndex2017[[i]]$TimeSeriesData$Hourly.Out
      timevec <- FluidsIndex2017[[i]]$TimeSeriesData$time
      timevec <- timevec[!is.na(diff)]
      diff <- diff[!is.na(diff)]
      diff <- diff[!is.na(timevec)]
      timevec <- timevec[!is.na(timevec)] 
      
      PostopFluidsDataMatrix[i ,1] <- c(diff[which(!is.na(diff))[1] ])
      PostopFluidsDataMatrix[i ,2] <- FluidsIndex2017[[i]]$TimeSeriesData$Filter[1]
  }
  rownames(PostopFluidsDataMatrix) <- names(FluidsIndex2017)
  colnames(PostopFluidsDataMatrix) <- c('FluidBalance', 'CVVHDialysis')
return(PostopFluidsDataMatrix)
}
POM_MeanImputation <- function(MasterPreOpData){
  for(ii in 1:dim(MasterPreOpData)[2]){
    if(is.numeric(MasterPreOpData[, ii]) ){
     MasterPreOpData[is.na(MasterPreOpData[, ii]), ii] <- mean(MasterPreOpData[!is.na(MasterPreOpData[, ii]), ii])
    }
  }
  return(MasterPreOpData)
}
POM_GaussianImputation <- function(MasterPreOpData){
for(ii in 1:dim(MasterPreOpData)[2]){
    if(is.numeric(MasterPreOpData[, ii]) ){
      if(sum(is.na(MasterPreOpData[, ii]))> 0){
      MasterPreOpData[is.na(MasterPreOpData[, ii]), ii] <- rnorm(sum(is.na(MasterPreOpData[, ii])) , mean(MasterPreOpData[!is.na(MasterPreOpData[, ii]), ii]) , sqrt(var(MasterPreOpData[!is.na(MasterPreOpData[, ii]), ii])))
      }
      }
  }
  return(MasterPreOpData)
}
POM_SampledImputation <- function(MasterPreOpData){
  for(ii in 1:dim(MasterPreOpData)[2]){
      if( (sum(is.na(MasterPreOpData[, ii]))> 0) && (sum(is.na(MasterPreOpData[, ii]))< dim(MasterPreOpData)[1]) ){
      MasterPreOpData[is.na(MasterPreOpData[, ii]), ii] <- sample(x = MasterPreOpData[!is.na(MasterPreOpData[, ii]), ii] , size = sum(is.na(MasterPreOpData[, ii])) , replace = T )
      }
  }
  return(MasterPreOpData)
}
POM_RemoveOutliers <- function(MasterPreOpData , threshold = 5){
  for(ii in 1:dim(MasterPreOpData)[2]){
    if(is.numeric(MasterPreOpData[, ii]) ){
     logicalVector = abs(DP_NormaliseData(MasterPreOpData[,ii]) ) > threshold
     logicalVector[is.na(logicalVector)] <- FALSE
     MasterPreOpData[logicalVector , ii] <- NA
     }
  }
  return(MasterPreOpData)
}
POM_ExtraPostOpScores <- function(PostOpScoresIndex2017){
  output <- matrix(NA , length(PostOpScoresIndex2017) , 3)
  for(ii in 1:length(PostOpScoresIndex2017)){
    output[ii,] <- as.numeric(PostOpScoresIndex2017[[ii]][1 , c(2:4)])
  }
  colnames(output) <- names(PostOpScoresIndex2017[[1]][1 , c(2:4)]) 
  rownames(output) <- names(PostOpScoresIndex2017)
  return(output)
}
POM_TreatmentBeforeAFib <- function(FluidsIndex2017 , PatIndex2017){
  output <- matrix(NA , length(FluidsIndex2017) , 3)

  colnames(output) <- c('Noradrenaline' , 'Dopamine' , 'Adrenaline')
  rownames(output) <- names(FluidsIndex2017)
  
  NoradrenalineIndicies <- which(grepl('Noradrenaline' , names( FluidsIndex2017$`z1026-1-1-1`$TimeSeriesData)))
  DopanineIndicies <- which(grepl('Dopamine' , names( FluidsIndex2017$`z1026-1-1-1`$TimeSeriesData)))
  adremalineIndicies <- which(grepl('Adrenaline' , names( FluidsIndex2017$`z1026-1-1-1`$TimeSeriesData)))
  
  for(ii in 1:length(FluidsIndex2017)){
    output[ii,1] = sum(FluidsIndex2017[[ii]]$TimeSeriesData[1,NoradrenalineIndicies] , na.rm = T) > 0
    output[ii,2] = sum(FluidsIndex2017[[ii]]$TimeSeriesData[1,DopanineIndicies] , na.rm = T) > 0
    output[ii,3] = sum(FluidsIndex2017[[ii]]$TimeSeriesData[1,adremalineIndicies] , na.rm = T) > 0
    
  }
  
  return(output)
}
POM_FillEmptyRecords <- function(PatientData){
  PatientData$K[PatientData$K == ''] = NA
  for( i in 2:dim(PatientData)[1] ){
    PatientData[i,is.na(PatientData[i,])] <- PatientData[i-1,is.na(PatientData[i,])]
    
  }
  return(PatientData)
}
POM_FindFirstGoodRecord <- function(PatientData){
  dayOfRecord <- as.Date.POSIXct(PatientData$time)
  DayOne <- which((dayOfRecord - dayOfRecord[1]) ==1)
  
  #if(sum(!is.na(PatientData$Platelets[DayOne])*!is.na(PatientData$Bilirubin[DayOne])*!is.na(PatientData$Creatinine[DayOne])) == 0 ){
    output <- PatientData$time[DayOne[length(DayOne)]] #<- PatientData$time[DayOne[1]]
  #}
  #if(sum(!is.na(PatientData$Platelets[DayOne])*!is.na(PatientData$Bilirubin[DayOne])*!is.na(PatientData$Creatinine[DayOne])) >= 1 ){
  #  output <- PatientData$time[which(!is.na(PatientData$Platelets[DayOne])*!is.na(PatientData$Bilirubin[DayOne])*!is.na(PatientData$Creatinine[DayOne])>0)[1]]
  #}
  return(output)
}
POM_ExtractFirstGoodRecord <- function(PatientData){
  PatientData <- POM_FillEmptyRecords( PatientData )
  timestamp <- POM_FindFirstGoodRecord( PatientData )
  
  output <- PatientData[which.min(abs(PatientData$time - timestamp)) , ]
  return(output)
}
POM_ExtractFirstRecords <- function(AllDataStructure){
  output <- matrix(NA , length(AllDataStructure) ,  length(AllDataStructure[[1]]$DiscreteData) + length(names(AllDataStructure[[1]]$timeseriesvariables)) )
  colnames(output) <- c(names(AllDataStructure[[1]]$DiscreteData) , names(AllDataStructure[[1]]$timeseriesvariables))
  output <- data.frame(output)
  for(ii in 1:dim(output)[1]){
    output[ii,1:(length(AllDataStructure[[1]]$DiscreteData))] <- AllDataStructure[[ii]]$DiscreteData
    if(nrow(POM_ExtractFirstGoodRecord(PatientData = AllDataStructure[[ii]]$timeseriesvariables))> 0 ){
      output[ii,(length(AllDataStructure[[1]]$DiscreteData) +1):(length(AllDataStructure[[1]]$DiscreteData)+ length(names(AllDataStructure[[1]]$timeseriesvariables)))] <- POM_ExtractFirstGoodRecord(PatientData = AllDataStructure[[ii]]$timeseriesvariables)
    }  
    DP_WaitBar(ii/dim(output)[1])
  }
  return(output)
}
POM_CategoricalCreateFieldNames <- function(){
  return(c('VariableNames' , 'n' , 'nNAFib' , 'nAFib' , 'PerNAFib' , 'PerAFib' ,'PValueProp','LO' , 'L_CI' , 'U_CI' , 'PValue' ))
}
POM_CategoricalUnivariateAnalysis <- function(CategoricalData , AFLogical , MasterData , NameofVariable ){
  
  CategoricalData <- as.factor(CategoricalData)
  
  output <- data.frame(matrix(NA , length(levels(CategoricalData)) , 11 ))
  output[,1] <- levels(CategoricalData)
  names(output) <- POM_CategoricalCreateFieldNames()
  
  for(i in 1:length(levels(CategoricalData))){
    # basic analysis  
    tmp <- as.factor(CategoricalData == levels(CategoricalData)[i] )
    output[i , 2] <- as.numeric(summary(tmp)[2])
    output[i , 3] <- as.numeric(summary(tmp[AFLogical == 0])[2])
    output[i , 4] <- as.numeric(summary(tmp[AFLogical == 1])[2])
    output[i , 5] <- summary(tmp[AFLogical == 0])[2]/sum(AFLogical == 0)*100
    output[i , 6] <- summary(tmp[AFLogical == 1])[2]/sum(AFLogical == 1)*100
    testtmp <- prop.test( c(summary(tmp[AFLogical == 1])[1] ,summary(tmp[AFLogical == 0])[1]) , c(  sum(AFLogical == 1), sum(AFLogical == 0) ))
    output[i , 7] <- signif(testtmp$p.value , 3)
    # logistic regression analysis
  }
  
  model <- glm(FB_CreateFormula('AFLogical' , unique(c(which(names(MasterData) == NameofVariable), which(names(MasterData) == 'Valve'), which(names(MasterData) == 'CABG'), which(names(MasterData) == 'Aortic') ,  which(names(MasterData) == 'Complex') ,which(names(MasterData) == 'LogisticEUROScore'),which(names(MasterData) == 'CPB')  ) ) , MasterData ) ,
               family=binomial(link='logit'),
               data = MasterData)
  
  
  output[ , 8] <- exp(coef(model)[1:length(levels(CategoricalData))])
  output[ , 9:10] <- exp(confint(model, trace = FALSE) )[1:length(levels(CategoricalData)),]
  output[ , 11] <- summary(model)$coefficients[1:length(levels(CategoricalData)),4]  
  output[ , 2:11] <- signif(output[,2:11] , 4)
  
  output[1 , 8:11] <- NA
  return(output)
}
POM_ContinuousCreateFieldNames <- function(){
  return(c('VariableNames' , 'n' , 'M-NAF' , 'SD-NAF' , 'M-AF' , 'SD-AF' ,'P-M' , 'LO' , 'L_CI' , 'U_CI' , 'PValue' ))
}
POM_ContinuousUnivariateAnalysis <- function(ContinuousData , AFLogical , MasterData , NameofVariable){
  ContinuousData <- as.numeric(ContinuousData)
  output <- data.frame(matrix(NA , 1 , length(POM_ContinuousCreateFieldNames()) ))
  output[,1] <- NameofVariable
  names(output) <- POM_ContinuousCreateFieldNames()
  
  output[,2] <- sum(!is.na(ContinuousData))
  output[,3] <- mean(ContinuousData[AFLogical ==0] , na.rm = T) 
  output[,4] <- sqrt(var(ContinuousData[AFLogical ==0] , na.rm = T))
  output[,5] <- mean(ContinuousData[AFLogical ==1] , na.rm = T) 
  output[,6] <- sqrt(var(ContinuousData[AFLogical ==1] , na.rm = T))
  
  tmptest <- t.test(ContinuousData[AFLogical ==1] , ContinuousData[AFLogical ==0])
  output[,7] <- tmptest$p.value
  
  model <- glm(FB_CreateFormula('AFLogical' , unique(c(which(names(MasterData) == NameofVariable), which(names(MasterData) == 'ProcDetails') ,which(names(MasterData) == 'LogisticEUROScore'),which(names(MasterData) == 'CPB')  ) ) , MasterData ) ,
               family=binomial(link='logit'),
               data = MasterData)
  output[ , 8] <- exp(coef(model))[2]
  output[ , 9:10] <- exp(confint(model) )[2,]
  output[ , 11] <- summary(model)$coefficients[2,4]  
  output[ , 2:11] <- signif(output[,2:11] , 4)
  
  return(output)
}
POM_GroupComparison <- function(NamesofVariables , MasterData , ControlModel = 0){
  set.seed(1)
  IndiciesVariables <- which(names(MasterData) %in% NamesofVariables)
  formulaformodel <- FB_CreateFormula('AFLogical' , IndiciesVariables , MasterData)
  
  DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  ,  ]))
  model <- (glm(formula = formulaformodel,family=binomial(link='logit') , data=DataForLogistic))
  stepaicoutput <- stepAIC(model)
  
  if(is.na(all.vars(stepaicoutput$formula)[2])){
  AUC1 = 0
  }else{
  AUC1 <- FC_CalculateCrossValidatedROC( formulaformodel = stepaicoutput$formula , PreoperativeIndices = IndiciesVariables , MasterData )
  }
  StepOutputAUC <- FC_StepWiseForwardAUC(PreoperativeIndices =  IndiciesVariables , MasterData =MasterData,matrixforstep = diag(length(IndiciesVariables)) )
  AUC2 <- StepOutputAUC[[1]]
  #StepOutputAUC <- list(NA,NA)
  #AUC2 <- NA
  
  if(ControlModel!= 0){
  formulaformodel <- FB_CreateFormula('AFLogical' , which(names(MasterData) == ControlModel) , MasterData)
  AUC3 <- FC_CalculateCrossValidatedROC(formulaformodel ,which(names(MasterData) ==ControlModel), MasterData)
  }else{
    AUC3 <- NA
  }
  
  output <- list(AUC1 , AUC2 , AUC3 , stepaicoutput$formula ,  StepOutputAUC[[2]] )
  output <- setNames(output , c('AUCAIC' , 'AUCStepAUC' , 'AUCLogisticEuro' , 'ModelAIC' , 'ModelAUC'))
  return(output)
}
POM_CalulateProbabilitiesFromModel <- function(formula , DataForLogistic){
  model <- glm(formula = formula
               , family = binomial(link = "logit"),  data = DataForLogistic)
  return(LogisticProbility <- predict(model , DataForLogistic , type = c('response')))
  
}
