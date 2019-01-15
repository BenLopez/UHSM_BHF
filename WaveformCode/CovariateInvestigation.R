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


{DP_LoadPatientDendrite()
DetIndex2017 <- DP_ExtractRecordsFromDendrite(DetIndex2017, listAllPatients)
DetIndex2017 <- DetIndex2017[!duplicated(DetIndex2017$NewPseudoId) , ]
DetIndex2017 <- DetIndex2017[order(DetIndex2017$NewPseudoId) , ]

SetofCovariatesDet <- c("LogisticEUROScore" ,
                        "AdditiveEUROScore",
                        "SCTSLogisticEuroSCORE")

ReducedDetIndex2017 <- DetIndex2017[ , which(names(DetIndex2017) %in% SetofCovariatesDet)]}

{
ReducedPatientIndex <- DP_ExtractMetaDataForMultiplePatinets(PatIndex2017 = PatIndex2017 , listAllPatients = listAllPatients)
  SetofCovariates <- c( 'PseudoId',
                      'Age' ,
                      'Gender' ,
                      'Weight' , 
                      'AKIReason' , 
                      'ProcDetails' ,
                      'CPB',
                      'ConfirmedFirstNewAF',
                      'FirstNewAF')
  
ReducedPatientIndex <- ReducedPatientIndex[ , which(names(ReducedPatientIndex) %in% SetofCovariates)]
ReducedPatientIndex$Gender <- as.factor(ReducedPatientIndex$Gender)
ReducedPatientIndex$AKIReason[!is.na(ReducedPatientIndex$AKIReason)] <- 'Yes' 
ReducedPatientIndex$AKIReason[is.na(ReducedPatientIndex$AKIReason)] <- 'No'
ReducedPatientIndex$AKIReason <- as.factor(ReducedPatientIndex$AKIReason)
ReducedPatientIndex$ProcDetails <- as.factor(ReducedPatientIndex$ProcDetails)
ReducedPatientIndex$ConfirmedFirstNewAF[!is.na(ReducedPatientIndex$FirstNewAF)] <- 'Yes'
ReducedPatientIndex$ConfirmedFirstNewAF[is.na(ReducedPatientIndex$FirstNewAF)] <- 'No'
ReducedPatientIndex$ConfirmedFirstNewAF <- as.factor(ReducedPatientIndex$ConfirmedFirstNewAF)
ReducedPatientIndex$CPB[!is.na(ReducedPatientIndex$CPB)] <- 'Yes'
ReducedPatientIndex$CPB[is.na(ReducedPatientIndex$CPB)] <- 'No'
ReducedPatientIndex$CPB <- as.factor(ReducedPatientIndex$CPB)


ReducedPatientIndex <- ReducedPatientIndex[-which(ReducedPatientIndex$PseudoId == 'z973')[2] , ]
ReducedPatientIndex <- ReducedPatientIndex[!duplicated(ReducedPatientIndex$PseudoId) , ]
ReducedPatientIndex <- ReducedPatientIndex[order(ReducedPatientIndex$PseudoId) , ]
}

CompleteReduced <- cbind(ReducedPatientIndex , ReducedDetIndex2017)
CompleteReduced$AdditiveEUROScore[is.na(CompleteReduced$AdditiveEUROScore)] <- median(CompleteReduced$AdditiveEUROScore[!is.na(CompleteReduced$AdditiveEUROScore)])
CompleteReduced$SCTSLogisticEuroSCORE[is.na(CompleteReduced$SCTSLogisticEuroSCORE)] <- median(CompleteReduced$SCTSLogisticEuroSCORE[!is.na(CompleteReduced$SCTSLogisticEuroSCORE)])

summary(glm(formula = ConfirmedFirstNewAF ~ Age  + AKIReason + LogisticEUROScore + SCTSLogisticEuroSCORE + AdditiveEUROScore ,family=binomial(link='logit') , data=CompleteReduced))

model <- glm(formula = ConfirmedFirstNewAF ~ Age  + AKIReason + LogisticEUROScore + SCTSLogisticEuroSCORE+ AdditiveEUROScore, family=binomial(link='logit') , data=CompleteReduced)
summary(model)

tmp <- hist(model$fitted.values[CompleteReduced$ConfirmedFirstNewAF == 'No'] , freq = FALSE , col = rgb(0,0,1,alpha = 0.5) , breaks = seq(0,1,0.05)  , xlab = 'PAmp' , main = c('Logistic Regression Fit') , ylim =c(0,15))
hist(model$fitted.values[CompleteReduced$ConfirmedFirstNewAF == 'Yes'] , freq = FALSE , add = T , col = rgb(1,0,0,alpha = 0.5) , breaks = tmp$breaks)

x11()
par(mfrow = c(1 , 3))
plot(CompleteReduced$ConfirmedFirstNewAF , CompleteReduced$AKIReason  , xlab=c('Exhibit AF in Future') , ylab=c('AKI')  )
title('Acute Kidney Injury')
plot(CompleteReduced$ConfirmedFirstNewAF , CompleteReduced$CPB  , xlab=c('Exhibit AF in Future') , ylab=c('CPB')  )
title('Cornary-By-Pass')
plot(CompleteReduced$ConfirmedFirstNewAF , CompleteReduced$Gender  , xlab=c('Exhibit AF in Future') , ylab=c('Gender')  )
title('Gender')


{
x11()
par(mfrow = c(1 , 4))
BC_PlotCompareSingleHists(CompleteReduced$Age[CompleteReduced$ConfirmedFirstNewAF == 'No'] , CompleteReduced$Age[CompleteReduced$ConfirmedFirstNewAF == 'Yes'] , xlab = 'Age' , main = c('Age') )
BC_PlotCompareSingleHists(CompleteReduced$LogisticEUROScore[CompleteReduced$ConfirmedFirstNewAF == 'No'] , CompleteReduced$LogisticEUROScore[CompleteReduced$ConfirmedFirstNewAF == 'Yes'] , xlab = 'LogisticEUROScore' , main = c('LogisticEUROScore') )
BC_PlotCompareSingleHists(CompleteReduced$AdditiveEUROScore[CompleteReduced$ConfirmedFirstNewAF == 'No'] , CompleteReduced$AdditiveEUROScore[CompleteReduced$ConfirmedFirstNewAF == 'Yes'] , xlab = 'AdditiveEUROScore' , main = c('AdditiveEUROScore') )
BC_PlotCompareSingleHists(CompleteReduced$SCTSLogisticEuroSCORE[CompleteReduced$ConfirmedFirstNewAF == 'No'] , CompleteReduced$SCTSLogisticEuroSCORE[CompleteReduced$ConfirmedFirstNewAF == 'Yes'] , xlab = 'SCTSLogisticEuroSCORE' , main = c('SCTSLogisticEuroSCORE') )
tmp <- CompleteReduced[CompleteReduced$ConfirmedFirstNewAF == 'No' , c(1,10,11,12)]
BC_PlotPairsFromTwoVariables( CompleteReduced[CompleteReduced$ConfirmedFirstNewAF == 'Yes' , c(1,10,11,12)] ,tmp[sample(1:dim(tmp)[1] , sum(CompleteReduced$ConfirmedFirstNewAF == 'Yes')), ], alpha = 0.2 )
}

AFData <- CompleteReduced[CompleteReduced$ConfirmedFirstNewAF == 'Yes' , c(1,6,10,12)]
AFData$AKIReason <- as.numeric(AFData$AKIReason)

NAFData <- CompleteReduced[CompleteReduced$ConfirmedFirstNewAF == 'No' , c(1,6,10,12)]
NAFData$AKIReason <- as.numeric(NAFData$AKIReason)

BC_PlotCompareTwoHists(AFData , NAFData)


mAF <- apply(AFData , 2 , mean) 
vAF <- cov(AFData)

mNAF <- apply(NAFData , 2 , mean) 
vNAF <- cov(NAFData)

TotalData <- rbind(AFData , NAFData)

ImAF  <- sqrt(apply(TotalData,   1 , function(X){(X - mAF)%*%(solve(vAF)%*%(X - mAF))}))
ImNAF <- sqrt(apply(TotalData,   1 , function(X){(X - mNAF)%*%(solve(vNAF)%*%(X - mNAF))}))
ImAF <- (ImAF - mean(ImAF[1:dim(AFData)[1]]))/sqrt(var(ImAF - mean(ImAF[1:dim(AFData)[1]])))
ImNAF <- (ImNAF - mean(ImNAF[(dim(AFData)[1] + 1):(dim(TotalData)[1])]))/sqrt(var(ImNAF - mean(ImNAF[(dim(AFData)[1] + 1):(dim(TotalData)[1])])))

BC_PlotCompareTwoHists(cbind(ImAF[1:dim(AFData)[1]] , ImNAF[1:dim(AFData)[1]]) , cbind(ImAF[(dim(AFData)[1] + 1):(dim(TotalData)[1])] , ImNAF[(dim(AFData)[1] + 1):dim(TotalData)[1]]) )
BC_PlotPairsFromTwoVariables( cbind(ImAF[1:dim(AFData)[1]] , ImNAF[1:dim(AFData)[1]]) , cbind(ImAF[(dim(AFData)[1] + 1):(dim(TotalData)[1])] , ImNAF[(dim(AFData)[1] + 1):dim(TotalData)[1]])[sample(1:(dim(TotalData)[1] - dim(AFData)[1]) , dim(AFData)[1]) , ] , alpha = 0.2 )


AFModel <- kde(cbind(ImAF[1:dim(AFData)[1]] , ImNAF[1:dim(AFData)[1]]))
NAFModel <- kde(cbind(ImAF[(dim(AFData)[1] + 1):(dim(TotalData)[1])] , ImNAF[(dim(AFData)[1] + 1):dim(TotalData)[1]]))   

f_1 <- matrix(0 , dim(TotalData)[1] , 2)

f_1[,1] <- predict(AFModel , x = cbind(ImAF , ImNAF))
f_1[,2] <- predict(NAFModel , x = cbind(ImAF , ImNAF))

alpha = c(sum(CompleteReduced$ConfirmedFirstNewAF == 'Yes')/length(CompleteReduced$ConfirmedFirstNewAF) , (1-(sum(CompleteReduced$ConfirmedFirstNewAF == 'Yes')/length(CompleteReduced$ConfirmedFirstNewAF))))
PosteriorProbability <- alpha[1]*f_1[,1]/(alpha[1]*f_1[,1] + alpha[2]*f_1[,2])

plot(PosteriorProbability)
abline(0.2 , 0)
abline(0.5 , 0)

{
tmp <- BC_CreateProbCalibrationStruct(as.matrix(PosteriorProbability) , alpha , length(CompleteReduced$ConfirmedFirstNewAF))
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(tmp, BinWidth = 0.075))

CovariateEmulationClass = BE_CreateDefaultEmulationClass()
CovariateEmulationClass$X = as.matrix(EmpericalProbabilityStructure$x)
CovariateEmulationClass$Y = as.matrix(EmpericalProbabilityStructure$y - EmpericalProbabilityStructure$x)
CovariateEmulationClass$MeanFunction = function(X){
  X <- as.matrix(X)
  H = cbind( as.matrix(1 + 0*X[,1]) , X  )
  return(H)
}
CovariateEmulationClass$CorrelationLength <- function(X , n){
  return(c(0.2 , 0.2))
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
BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Covariate Emperical Probabilities')

x11(20,14)
E_D_Prevision <- as.matrix(PosteriorProbability + EmOutput$E_D_fX)
E_D_Prevision[E_D_Prevision > 1] <- 0.99999999
E_D_Prevision[E_D_Prevision < 0] <- 0.00000001

tmp <- BC_CreateProbCalibrationStruct(E_D_Prevision , alpha , length(CompleteReduced$ConfirmedFirstNewAF))
CalibratedProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(tmp, BinWidth = 0.075))
BC_PlotCreateProbabilityCalibrationPlot(CalibratedProbabilityStructure) + ggtitle('Calibrated Covariate Emperical Probabilities')
}

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = tmp)

ROCplot <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('Covariate ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Covariate PPV vs NPV Curves')

x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot , nrow= 2)



