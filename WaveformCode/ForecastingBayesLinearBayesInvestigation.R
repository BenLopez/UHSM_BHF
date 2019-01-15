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

DataStructure[, 3, ]  <-  log(DataStructure[, 3, ])
DataStructure[ , 4, ] <-  sqrt(DataStructure[ , 4, ])
DataStructure[ , 7, ] <-  log(DataStructure[ , 7, ])
DataStructure[ , 10, ] <- sqrt(DataStructure[ , 10, ])
DataStructure[ , 11, ] <- sqrt(DataStructure[ , 11, ])
DataStructure[ , 12, ] <- sqrt(DataStructure[ , 12, ])


#indextoview <- (dim(DataStructure)[3] - 29) : (dim(DataStructure)[3])
indextoview <- (dim(DataStructure)[3] - 28) 

{
  x11(20,14)
  par(mfrow = c(2,6))
  BC_PlotCompareSingleHists(X = DataStructure[AFLogical ==F , 1, indextoview] , Y = DataStructure[AFLogical == T , 1, indextoview] , main = 'Mean')
  BC_PlotCompareSingleHists(X = DataStructure[AFLogical ==F , 2, indextoview] , Y = DataStructure[AFLogical == T , 2, indextoview], main = 'Median')
  BC_PlotCompareSingleHists(X = (DataStructure[AFLogical ==F , 3, indextoview]) , Y = (DataStructure[AFLogical == T , 3, indextoview]), main = 'Variance' )
  BC_PlotCompareSingleHists(X = (DataStructure[AFLogical ==F , 4, indextoview]) , Y = (DataStructure[AFLogical == T , 4, indextoview]), main = 'IQR' )
  BC_PlotCompareSingleHists(X = DataStructure[AFLogical ==F , 5, indextoview] , Y = DataStructure[AFLogical == T , 5, indextoview], main = 'Skewness')
  BC_PlotCompareSingleHists(X = DataStructure[AFLogical ==F , 6, indextoview] , Y = DataStructure[AFLogical == T , 6, indextoview], main = 'Kurtosis' )
  BC_PlotCompareSingleHists(X = (DataStructure[AFLogical ==F , 7, indextoview]) , Y = (DataStructure[AFLogical == T , 7, indextoview]), main = 'max Denisty' )
  BC_PlotCompareSingleHists(X = (DataStructure[AFLogical ==F , 8, indextoview]) , Y = (DataStructure[AFLogical == T , 8, indextoview]), main = 'Variance Denisty' )
  BC_PlotCompareSingleHists(X = DataStructure[AFLogical ==F , 9, indextoview] , Y = DataStructure[AFLogical == T , 9, indextoview] , main = 'Number of Modes')
  BC_PlotCompareSingleHists(X = (DataStructure[AFLogical ==F , 10, indextoview]) , Y = (DataStructure[AFLogical == T , 10, indextoview]), main = 'Variance in Mode Denisty' )
  BC_PlotCompareSingleHists(X = (DataStructure[AFLogical ==F , 11, indextoview]) , Y = (DataStructure[AFLogical == T , 11, indextoview]), main = 'Variance in Mode Location' )
  BC_PlotCompareSingleHists(X = (DataStructure[AFLogical ==F , 12, indextoview]) , Y = (DataStructure[AFLogical == T , 12, indextoview]), main = 'Approximate Entropy' )
  
  BC_PlotPairsFromTwoVariables(DataStructure[which(AFLogical ==F)[sample(1:sum(AFLogical  == 0) , sum(AFLogical == 1))] , -c( 13) , indextoview] ,DataStructure[AFLogical == T , -c( 13), indextoview]  , alpha = 0.1)
  
  
  mAF <- (apply(DataStructure[AFLogical == T , -c( 13), indextoview]  , 2 , mean))
  CAF <- cov(DataStructure[AFLogical == T , -c( 13), indextoview])
  
  mNAF <- (apply(DataStructure[AFLogical == F , -c( 13), indextoview]  , 2 , mean))
  CNAF <- cov(DataStructure[AFLogical == F , -c( 13), indextoview])
  
  ImAF <- BLBF_ForecastingCalculateImplausabilities(DataStructure[, -c( 13), indextoview] , mAF , CAF)
  ImNAF <- BLBF_ForecastingCalculateImplausabilities(DataStructure[, -c( 13), indextoview] , mNAF , CNAF)
  
  ImAF <- (ImAF - mean(ImAF[AFLogical ==T])) / sqrt(var(ImAF[AFLogical ==T])) 
  ImNAF <- (ImNAF - mean(ImNAF[AFLogical ==F])) / sqrt(var(ImNAF[AFLogical == F])) 
  
  x11()
  par(mfrow =c(1,2))
  print(BC_PlotCompareSingleHists((ImAF[AFLogical ==F]) , (ImAF[AFLogical ==T])  ,breaks = 15 ,  main = 'Implausabilty of Going into AF', xlab = 'Im'))
  print(BC_PlotCompareSingleHists((ImNAF[AFLogical ==F]) , (ImNAF[AFLogical ==T]) ,breaks = 15 , main = 'Implausabilty of Not Going into AF' , xlab = 'Im'))
  
  AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  x11(20,14)
  plot(NAFModel  , col ='blue', xlab = 'Implausabilty of Going into AF' , ylab = 'Implausabilty of Not Going into AF' , main = 'KDE Implausabilities')
  plot(AFModel, add = T, col ='red')
  
  points(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F] , col = 'blue' , pch = 16 )
  points(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T], col = 'red' ,pch = 16 )
  
  PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical)
  
  
  x11(20,14)
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , AFLogical == T), BinWidth = 0.1))
  print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Forecasting Emperical Probabilities'))
  
  PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability , AFLogical == T))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('Forecasting ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Forecasting PPV vs NPV Curves')
  x11(20,14)
  grid.arrange(ROCplot , NPVPPVPlot , nrow= 2)
  }

{
  StartIndex <- 371
  NumberofPoints <- 1
  VariablesToView <- c(1:12)[c(-1 , -3 )]
  LengthVector <- length(VariablesToView)
  
  tmpdata <- DataStructure[ , VariablesToView, c(StartIndex:(StartIndex+NumberofPoints) )]
  AFSOS <- BLBF_CalculateSecondOrderStruct( tmpdata[which(AFLogical)[sort(sample(1:sum(AFLogical) , sum(AFLogical) -1))] , , ])
  NAFSOS <- BLBF_CalculateSecondOrderStruct( tmpdata[which(AFLogical == F)[sort(sample(1:sum(AFLogical==F) , sum(AFLogical==F) -1))] , , ])
  
  AdjustedVersionsAF <- t(apply( tmpdata , 1 , function(X){BLBF_CalculateForecastAdjustedVersion(AFSOS , X[1:(LengthVector*(NumberofPoints))] , X[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints + 1) )])}))
  AdjustedVersionsNAF <- t(apply( tmpdata , 1 , function(X){BLBF_CalculateForecastAdjustedVersion(NAFSOS , X[1:(LengthVector*(NumberofPoints))] , X[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints +1) )])}))
  
  AdjustedVersionsAF[AFLogical == T] <- BLBF_CrossValidateAdjustedVersionAll(tmpdata[AFLogical == T , , ]  , LengthVector , NumberofPoints )
  AdjustedVersionsNAF[AFLogical == F] <- BLBF_CrossValidateAdjustedVersionAll(tmpdata[AFLogical == F , , ]  , LengthVector , NumberofPoints )
  
  # Could be over fitting.
  AdjustedVersionsSOSAF <- BLU_CreateSOS(mu = apply(AdjustedVersionsAF[AFLogical == T ,] , 2 , mean) , Sigma = cov(AdjustedVersionsAF[AFLogical == T ,]) )
  AdjustedVersionsSOSNAF <- BLU_CreateSOS(mu = apply(AdjustedVersionsNAF[AFLogical == F ,]  , 2 , mean) , Sigma = cov(AdjustedVersionsNAF[AFLogical == F , ]) )

  ImAF <- apply(AdjustedVersionsAF , 1 , function(X){BLBF_Discrepancy(AdjustedVersionsSOSAF , X )})
  ImNAF <- apply(AdjustedVersionsNAF , 1 , function(X){BLBF_Discrepancy(AdjustedVersionsSOSNAF, X)} ) 
  ImAF <- (ImAF - mean(ImAF[AFLogical ==T])) / sqrt(var(ImAF[AFLogical ==T])) 
  ImNAF <- (ImNAF - mean(ImNAF[AFLogical ==F])) / sqrt(var(ImNAF[AFLogical == F])) 
  
  x11()
  par(mfrow = c(1 , 2))
  print(BC_PlotCompareSingleHists((ImAF[AFLogical ==F]) , (ImAF[AFLogical ==T]) ))
  print(BC_PlotCompareSingleHists((ImNAF[AFLogical ==F]) , (ImNAF[AFLogical ==T]) ))
  
  AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  
  x11(20,14)
  plot(NAFModel  , col ='red' , xlab =c('Implausability AF') , ylab =c('Implausability NAF'))
  plot(AFModel, add = T)
  
  points(ImAF[AFLogical == F] , ImNAF[AFLogical == F] , col = 'red' , pch = 16)
  points(ImAF[AFLogical == T] , ImNAF[AFLogical == T], pch = 16)
  
  PosteriorProbability3 <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical , alpha = PosteriorProbability)
  
  PerformanceSweep2 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability3 , AFLogical == T))
  ROCplot2 <- BC_PlotsCreateROC(PerformanceSweep2) + ggtitle('Forecasting ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot2 <- BC_PlotsCreateNPVPPV(PerformanceSweep2) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Forecasting PPV vs NPV Curves')
  
  x11(20,14)
  grid.arrange(ROCplot2 , NPVPPVPlot2 , nrow= 2)
  
  x11(20,14)
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability3 , AFLogical == T), BinWidth = 0.1))
  print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Forecasting Emperical Probabilities'))
  x11(20,14)
  plot(1-PerformanceSweep[,2] , PerformanceSweep[,1] , col = 'red' , type ='l')
  lines(1-PerformanceSweep2[,2] , PerformanceSweep2[,1]  , col = 'blue' , type ='l')
}

PosteriorProbability1 <- BLBF_CalculatePosteriorProbabilityNBLA( DataStructure[ ,VariablesToView <- c(1:12)[c(-1 , -3 )] , ] , AFLogical , VariablesToView ,  indextoview )
PosteriorProbability2 <- BLBF_CalculatePosteriorProbabilityBLA(DataStructure , AFLogical , PriorProbability = PosteriorProbability , StartIndex = 371 , NumberofPoints = 1 , VariablesToView = c(1:12)[c(-1 , -3 )])
PosteriorProbability3 <- BLBF_CalculatePosteriorProbabilityBLA(DataStructure , AFLogical , PriorProbability = PosteriorProbability2 , StartIndex = 371 , NumberofPoints = 2 , VariablesToView = c(1:12)[c(-1 , -3 )])

PerformanceSweep1 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability1 , AFLogical == T))
PerformanceSweep2 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability2 , AFLogical == T))
PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability3 , AFLogical == T))

x11()
plot( 1-PerformanceSweep1[ , 2] , PerformanceSweep1[ , 1] ,  col = rgb(0,0,1 , alpha = 0.5) , type='l' , ylim = c(0,1)  , xlab = c('1 - Specificity') , ylab = c('Sensitivity'))
lines( 1-PerformanceSweep2[ , 2] , PerformanceSweep2[ , 1] , col = rgb(1,0,0 , alpha = 0.5) ,  type ='l')
title('Blue: ROC curve Using One Time Point. Red: ROC Curve Using Two Time Points ')

PosteriorProbabilityStruct <- matrix(0  , length(AFLogical) , length((dim(DataStructure)[3] - 29) : (dim(DataStructure)[3])))
PerformanceSweepStruct <- array(0 , c(29 , dim(PerformanceSweep1)) )

PosteriorProbabilityStruct[,1] <- BLBF_CalculatePosteriorProbabilityNBLA( DataStructure[ ,VariablesToView <- c(1:12)[c(-1 , -3 )] , ] , AFLogical , VariablesToView ,  indextoview )
for(i in 2:29){
  PosteriorProbabilityStruct[,i] <- BLBF_CalculatePosteriorProbabilityBLA(DataStructure , AFLogical , PriorProbability = PosteriorProbabilityStruct[,i -1] , StartIndex = 371 +(i-2) , NumberofPoints = 1 , VariablesToView = c(1:12)[c(-1 , -3 )])
  DP_WaitBar(i/29)
}
PerformanceSweepStruct[1 , ,] <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbabilityStruct[,1] , AFLogical == T))
for(i in 2:29){
  PerformanceSweepStruct[i , ,] <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbabilityStruct[,i] , AFLogical == T))
}

x11(20,14)
XX <- flipud(as.matrix(seq(1,29,1))*-500 )
plot(XX , PosteriorProbabilityStruct[ 3, 1:29] ,  col = rgb(0,0,1 , alpha = 0.5) , type='l' , ylim = c(0,1) , ylab = 'P' , xlab = c('Beats Before AF') )
for(i in 4:length(AFLogical)){
if(AFLogical[i] == T){
lines(XX , PosteriorProbabilityStruct[ i, 1:29] ,  col = rgb(0,0,1 , alpha = 0.5) , type='l')} 
if(AFLogical[i] == F){
    lines(XX , PosteriorProbabilityStruct[ i, 1:29] ,  col = rgb(1,0,0 , alpha = 0.1) , type='l')} 
}
title('Probability Trajectories')

x11(20,14)
plot(1-PerformanceSweepStruct[1 ,  , 2] , PerformanceSweepStruct[1 ,  , 1] , type ='l' , col = rgb(0,0,1, alpha = 0.5), xlab = '1 -Specificity' , ylab = 'Sensitvity')
for(i in 1:29){
lines(1-PerformanceSweepStruct[i ,  , 2] , PerformanceSweepStruct[i ,  , 1] , type ='l' , col = rgb(0,0,1, alpha = 0.5))
}
abline(0,1)  
title('ROC Curves over Time')

x11(20,14)
plot(1-PerformanceSweep1[,2] , PerformanceSweep1[,1] , col = 'red' , type ='l' , xlab = '1 -Specificity' , ylab = 'Sensitvity')
lines(1-PerformanceSweep2[,2] , PerformanceSweep2[,1]  , col = 'blue' , type ='l')
lines(1-PerformanceSweep3[,2] , PerformanceSweep3[,1]  , col = 'green' , type ='l')
abline(0,1)
title('ROC Curves over Time')



PosteriorProbabilityStruct <- matrix(0  , length(AFLogical) , length((dim(DataStructure)[3] - 29) : (dim(DataStructure)[3])))
counter <- 1
for(indextoview in (dim(DataStructure)[3] - 29) : (dim(DataStructure)[3])){
  
  mAF <- (apply(DataStructure[AFLogical == T , -c( 13), indextoview]  , 2 , mean))
  CAF <- cov(DataStructure[AFLogical == T , -c( 13), indextoview])
  mNAF <- (apply(DataStructure[AFLogical == F , -c( 13), indextoview]  , 2 , mean))
  CNAF <- cov(DataStructure[AFLogical == F , -c( 13), indextoview])
  
  ImAF <- BLBF_ForecastingCalculateImplausabilities(DataStructure[ , -c( 13), indextoview] , mAF , CAF)
  ImNAF <- BLBF_ForecastingCalculateImplausabilities(DataStructure[ , -c( 13), indextoview] , mNAF , CNAF)
  
  ImAF <- (ImAF - mean(ImAF[AFLogical ==T])) / sqrt(var(ImAF[AFLogical ==T])) 
  ImNAF <- (ImNAF - mean(ImNAF[AFLogical ==F])) / sqrt(var(ImNAF[AFLogical == F])) 
  
  AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical)
  PosteriorProbabilityStruct[ , counter ] <- PosteriorProbability
  
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , AFLogical == T), BinWidth = 0.1))
  #print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Forecasting Emperical Probabilities'))
  
  PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability , AFLogical == T))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('Forecasting ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Forecasting PPV vs NPV Curves')
  #x11(20,14)
  #grid.arrange(ROCplot , NPVPPVPlot , nrow= 2)
  counter <- counter + 1
  }

{plot(0, xlim = c(0,30) , ylim = c(0,1))
  for(i in 1:dim(PosteriorProbabilityStruct)[1]){
    if(AFLogical[i] == T){lines(DP_cummean(PosteriorProbabilityStruct[i,]) , col =  rgb(0,0,1,alpha = 0.5))}
    if(AFLogical[i] == F){lines(DP_cummean(PosteriorProbabilityStruct[i,])  , col = rgb(1,0,0,alpha = 0.05) )}  
  }
  
  
  tmp <- apply(PosteriorProbabilityStruct , 1 , function(X){DP_cummean(X)})
  lines(apply(tmp[ , AFLogical == F] , 1 , function(X){quantile(X , 0.90)}))
  
  requiredspecificity <- 0.9
  Threshold <- apply(tmp[ , AFLogical == F] , 1 , function(X){quantile(X , requiredspecificity)})
  
  senstruct <- matrix(0,30,1)
  Specstruct <- matrix(0,30,1)
  
  for(i in 1:30){
    senstruct[i] <- sum(tmp[i,AFLogical == T] > Threshold[i])/sum(AFLogical == T)
    Specstruct[i] <- sum(tmp[i,AFLogical == F] < Threshold[i])/sum(AFLogical == F)
  }
  
  plot(flipud(as.matrix(seq(1,30,1))*-500 ),senstruct , type = 'l' , ylim = c(0,1) , ylab=c('Sensitivity and Specificity') , xlab = c('Number of Beats Before AF') )
  lines(flipud(as.matrix(seq(1,30,1))*-500 ),Specstruct , type = 'l' , col ='red')
  title('Sensitivity and Specifictity Over Time')
}


