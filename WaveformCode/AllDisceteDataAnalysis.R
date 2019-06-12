{
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
    set.seed( 1 )
  }
}


BC_CalculateAreaUnderCurve <- function(PerformanceSweep){
  return(sum(diff(PerformanceSweep[,2])*PerformanceSweep[2:dim(PerformanceSweep)[1],1]))
}


##### Preoperation Data #####

# PatIndex2017, DetIndex2017 , BioChemIndex2017, HaemIndex2017 

{
AFMedication <- POM_ExtractAFMedication(FluidsIndex2017)
MasterPreOpData <- POM_CreateDataStructure(PatIndex2017 , DetIndex2017 , BioChemIndex2017)
PreOpHaem <- POM_ExtractPreOpFromHaem(HaemIndex2017)

PreOpHaem <- PreOpHaem[ which(rownames(PreOpHaem) %in% MasterPreOpData$NewPseudoId) ,  ]
PreOpHaem <- PreOpHaem[order(rownames(PreOpHaem)) ,  ]

MasterPreOpData <- cbind(MasterPreOpData , PreOpHaem)
AFMedication <- AFMedication[ which(AFMedication$NewPseudoId %in% MasterPreOpData$NewPseudoId) ,  ]
AFMedication <- AFMedication[ order(AFMedication$NewPseudoId) , ]

# Remove pre operative atrial-fibrillation/flutter
AFMedication <- AFMedication[MasterPreOpData$Pre_OperativeHeartRhythm != "Atrial fibrillation/flutter" , ]
MasterPreOpData <- MasterPreOpData[MasterPreOpData$Pre_OperativeHeartRhythm != "Atrial fibrillation/flutter" , ]

# Remove patinets who do not exhibit AFib but are treated with either amioderone or didoxin
tmplogical <- ((MasterPreOpData$AFLogical ==0)*(apply(AFMedication[ , 2:3] , 1, sum ) > 0)) == 0
MasterPreOpData <- MasterPreOpData[tmplogical, ]
AFMedication <- AFMedication[tmplogical , ]
rm(tmplogical)
MasterPreOpData$AFLogical <- as.factor( MasterPreOpData$AFLogical )
MasterPreOpData$Gender <- as.factor( MasterPreOpData$Gender )

{
  PostOpBioChem <- POM_ExtractPostOpFromBioChem(BioChemIndex2017 = BioChemIndex2017)
  PostOpHaem <- POM_ExtractPostOpFromHaem(HaemIndex2017 = HaemIndex2017)
  PostOpFluids <- POM_ExtractPostOpFromFluids(FluidsIndex2017 = FluidsIndex2017)
  
  PostOpBioChem <- PostOpBioChem[ which(rownames(PostOpBioChem) %in% MasterPreOpData$NewPseudoId) ,  ]
  PostOpBioChem <- PostOpBioChem[order(rownames(PostOpBioChem)) ,  ]
  
  PostOpHaem <- PostOpHaem[ which(rownames(PostOpHaem) %in% MasterPreOpData$NewPseudoId) ,  ]
  PostOpHaem <- PostOpHaem[order(rownames(PostOpHaem)) ,  ]
  
  PostOpFluids <- as.matrix(PostOpFluids[ which(rownames(PostOpFluids) %in% MasterPreOpData$NewPseudoId) ,  ])
  PostOpFluids <- PostOpFluids[order(rownames(PostOpFluids) ),   ]
}

MasterData <- cbind(MasterPreOpData , PostOpBioChem , PostOpHaem , PostOpFluids)
MasterData$ProcDetails <- DP_AssignSurgeryLabels(MasterData$ProcDetails)
MasterData$LogisticEUROScore <- log(MasterData$LogisticEUROScore)

MasterData$BMI = MasterData$Weight/(MasterData$Height/100)

}


# Patient demographic infomortion
{
{
sum(MasterData$AFLogical == 1)
sum(MasterData$AFLogical == 0)

# Age
t.test(MasterData$Age[MasterData$AFLogical == 1] , MasterData$Age[MasterData$AFLogical == 0])
sqrt(var(MasterData$Age[MasterData$AFLogical == 1]))
sqrt(var(MasterData$Age[MasterData$AFLogical == 0]))

# Gender
(summary(MasterData$Gender[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
(summary(MasterData$Gender[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
prop.test( c(summary(MasterData$Gender[MasterData$AFLogical == 1])[1] ,summary(MasterData$Gender[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) ))

# Weight
t.test(MasterData$Weight[MasterData$AFLogical == 1] , MasterData$Weight[MasterData$AFLogical == 0])
sqrt(var(MasterData$Weight[MasterData$AFLogical == 1]))
sqrt(var(MasterData$Weight[MasterData$AFLogical == 0]))

# BMI

t.test(MasterData$BMI[MasterData$AFLogical == 1] , MasterData$BMI[MasterData$AFLogical == 0])
sqrt(var(MasterData$BMI[MasterData$AFLogical == 1] , na.rm =T))
sqrt(var(MasterData$BMI[MasterData$AFLogical == 0], na.rm =T))


# Pre operative heart rhythm
tmp <- as.factor(MasterData$Pre_OperativeHeartRhythm == 'Sinus Rhythm')
(summary(tmp[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
(summary(tmp[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
prop.test( c(summary(tmp[MasterData$AFLogical == 1])[1] ,summary(tmp[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) ))


# Pre operative RRT

tmp <- as.factor(MasterData$PreopRRT == 'Yes')
(summary(tmp[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
(summary(tmp[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
prop.test( c(summary(tmp[MasterData$AFLogical == 1])[1] ,summary(tmp[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) ))
}

# Surgery 
ProcDetails <- as.factor(MasterData$ProcDetails)

for(i in 1:length(levels(ProcDetails))){
print(levels(ProcDetails)[i])
tmp <- as.factor(MasterData$ProcDetails == levels(ProcDetails)[i] )
print(summary(tmp[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
print(summary(tmp[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
print(prop.test( c(summary(tmp[MasterData$AFLogical == 1])[1] ,summary(tmp[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) )))
}  
}
# Surgery type boxplots
for(i in c(1,2,3,4,13:28)){
  plot(data.frame(as.factor(MasterData$ProcDetails) , MasterData[ , i]) , ylab = names(MasterData)[i] , xlab = 'ProxDetails' )
}


rownames(tteststruct) = names(MasterData)[c(2,4,13:28 )]

MasterData <- POM_RemoveOutliers(MasterPreOpData = MasterData)

NamesofPreOperativePredictors <- names(MasterData)[c(1:4,14:30 , 46) ]
counter <- 1
tteststruct <- matrix(0 , length( c(2,4,14:29 , 46 )) , 1)
for(i in c(2,4,14:29 , 46)){
BC_PlotCompareSingleHists(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i] , main =  names(MasterData)[i])
print(names(MasterData)[i])
print(t.test(MasterData[MasterData$AFLogical == 1 , i] , MasterData[MasterData$AFLogical == 0 , i]))
print(sqrt(var(MasterData[MasterData$AFLogical == 1 , i] , na.rm = T)))
print(sqrt(var(MasterData[MasterData$AFLogical == 0 , i] , na.rm = T)))
tteststruct[counter] <- (t.test(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i])$p.value)
counter <- counter + 1
}

# Logistic regression analysis individual variables.
DataForLogistic <- data.frame(MasterData[  , c(1,2,4,10,13:28)])

for(i in 2:dim(DataForLogistic)[2]){
  model <- (glm(formula = as.formula(paste0("AFLogical ~" , names(DataForLogistic)[i]) )
                ,family=binomial(link='logit') , data=DataForLogistic))
# odds ratio
# print( names(DataForLogistic)[i] )
  print(exp(coef(model)[2]))
  print( exp(confint(model)[2,] ))  
  print( summary( model )$coefficients[2,4] )
}
MasterData$ProcDetails <- as.factor(MasterData$ProcDetails)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,2,3,4,5,10,14:29 , 46)]))
model <- (glm(formula = AFLogical ~ 
                ProcDetails  +
                Age +
                Gender +
                Weight +
               BMI +
                LogisticEUROScore +
                PreOpNa +
                PreOpK +
                PreopUrea +
                PreopCreat +
                PreOpCRP +
                PreOpAlb +
                PreopBili +
                PreopMg+
                PreopHb +
                PreopPLT +
                PreopWBC +
                PreopAPTT+
                
                ,family=binomial(link='logit') , data=DataForLogistic))
summary(model)$coefficients[summary(model)$coefficients[,4] < 0.1,4]
summary( model )
stepAIC( model )

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,2,3,4,5,11,14:29 , 46)]))
model <- (glm(formula = AFLogical ~
                (ProcDetails == 'CABG' ) +
                Age +
                Gender +
                AdditiveEUROScore +
                PreOpNa +
                PreopUrea +
                PreOpAlb  +
                PreopAPTT
              , family=binomial(link='logit') , data=DataForLogistic))
summary(model)


LogisticProbility <- matrix(0,length(MasterData$AFLogical),1)
for(i in 1:dim(LogisticProbility)[1]){
  DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,2,3,4,5,11,14:29 , 46)]))
  model <- (glm(formula = AFLogical ~
                  (ProcDetails == 'CABG' ) +
                  Age +
                  Gender +
                  AdditiveEUROScore +
                  PreOpNa +
                  PreopUrea +
                  PreOpAlb  +
                  PreopAPTT
                , family=binomial(link='logit') , data=DataForLogistic[-i,]))
LogisticProbility[i,] <- predict(object= model , newdata = DataForLogistic[i , ] , type = c( "response") )
}

x11()
BC_PlotCompareSingleHists(LogisticProbility[MasterData$AFLogical ==0] , LogisticProbility[MasterData$AFLogical ==1] , main ='Histogram of Logistic Model Output'  )
t.test(LogisticProbility[MasterData$AFLogical ==0] ,  LogisticProbility[MasterData$AFLogical ==1])

x11()
plot(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$AFLogical ==0] , seq(0,1,0.01) ) , col = 'blue' , type ='l' , xlab = 'Output Logistic Model' , ylab = 'Culmulative Prob')
lines(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$AFLogical ==1] , seq(0,1,0.01) ) , col = 'red' , type ='l')
abline(v = 0.2)
title('CDF Analysis')

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$AFLogical == 1))
ROCplot2 <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot2 <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept =( 1-0.2)) + geom_hline(yintercept = 0.2)+  ggtitle('PPV vs NPV Curves')

x11(20,14)
grid.arrange(ROCplot2 , NPVPPVPlot2 + xlim(( 1-0.2) , 1) , nrow= 2)

x11(20,14)
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$AFLogical == 1), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))

# BayesLinearBayes Model

{
DataForBayes <- POM_SampledImputation(MasterData)
DataForBayes <- DataForBayes[ , names(MasterData) %in% c('ProcDetails' , 'Age' , 'AdditiveEUROScore' , 'PreOpNa' , 'PreopUrea' , 'PreopCreat' , 'PreOpAlb' , 'PreopBili')]
DataForBayes$ProcDetails <-DataForBayes$ProcDetails == 'CABG'
AFLogical <- as.numeric(MasterPreOpData$AFLogical)==2

PosteriorProbability <- BLBC_FitBayesLinearBayesClassifier(Data = as.matrix(DP_RemoveNaRows(DataForBayes) ) , Labels = AFLogical[DP_FindNARows(DataForBayes)] )

x11()
BC_PlotCompareSingleHists(PosteriorProbability[AFLogical ==0] ,PosteriorProbability[AFLogical ==1] , main ='Histogram of Logistic Model Output'  , breaks = 20)
t.test(PosteriorProbability[AFLogical ==0] , PosteriorProbability[AFLogical ==1])

x11()
plot(seq(0,1,0.01) , FM_CalculateCDFS(PosteriorProbability[AFLogical==0] , seq(0,1,0.01) ) , col = 'blue' , type ='l' , xlab = 'Output Logistic Model' , ylab = 'Culmulative Prob')
lines(seq(0,1,0.01) , FM_CalculateCDFS(PosteriorProbability[AFLogical ==1] , seq(0,1,0.01) ) , col = 'red' , type ='l')
abline(v = 0.2)
title('CDF Analysis')

PerformanceSweep2 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability ,AFLogical[DP_FindNARows(DataForBayes)] == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep2) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep2) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')

x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot + xlim(( 1-0.2) , 1) , nrow= 2)

x11(20,14)
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , MasterData$AFLogical == 1), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))

}

plot(1- PerformanceSweep2[,2] , PerformanceSweep2[,1] , type ='l' , col ='blue')
lines(1- PerformanceSweep[,2] , PerformanceSweep[,1] , type ='l' , col = 'red')
abline(v = 0.25)

plot(LogisticProbility[MasterData$AFLogical ==0] , PosteriorProbability[MasterData$AFLogical ==0] , col ='blue' , pch = 16 , xlim = c(0,1) , ylim = c(0,1))
points(LogisticProbility[MasterData$AFLogical ==1] , PosteriorProbability[MasterData$AFLogical ==1] , col ='red' , pch = 16 )

###### PostOp Analysis ######
PostOpIndexes <- c(1, 9 , 29:44)
PostOpNames <- names(MasterData)[PostOpIndexes]

NamesofPreOperativePredictors <- names(MasterData)[PostOpIndexes ]
counter <- 1
tteststruct <- matrix(0 , length( PostOpIndexes[2:length(PostOpIndexes)]) , 1)
for(i in PostOpIndexes[2:length(PostOpIndexes)]){
  BC_PlotCompareSingleHists(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i] , main =  names(MasterData)[i])
  print(names(MasterData)[i])
  print(t.test(MasterData[MasterData$AFLogical == 1 , i] , MasterData[MasterData$AFLogical == 0 , i]))
  print(sqrt(var(MasterData[MasterData$AFLogical == 1 , i] , na.rm = T)))
  print(sqrt(var(MasterData[MasterData$AFLogical == 0 , i] , na.rm = T)))
  tteststruct[counter] <- (t.test(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i])$p.value)
  counter <- counter + 1
}

DataForLogistic <- data.frame(MasterData[  , PostOpIndexes])

for(i in 2:dim(DataForLogistic)[2]){
  model <- (glm(formula = as.formula(paste0("AFLogical ~" , names(DataForLogistic)[i]) )
                ,family=binomial(link='logit') , data=DataForLogistic))
  # odds ratio
  # print( names(DataForLogistic)[i] )
  print(exp(coef(model)[2]))
  print( exp(confint(model)[2,] ))  
  print( summary( model )$coefficients[2,4] )
}

DataForLogistic <- POM_SampledImputation(data.frame(MasterData[  ,  ]))
model <- (glm(formula = AFLogical ~ 
                ProcDetails +
                Age +
                Gender +
                Weight +
                AdditiveEUROScore +
                LogisticEUROScore +
                SCTSLogisticEuroSCORE +
                PreOpNa +
                PreOpK +
                PreopUrea +
                PreopCreat +
                PreOpCRP +
                PreOpAlb +
                PreopBili +
                PreopMg+
                PreopHb +
                PreopPLT +
                PreopWBC +
                PreopAPTT +
                CPB +
                Na +
                K +
                Urea +
                Creatinine +
                CRP +
                Albumin +
                Bilirubin +
                Mg +
                WBC +
                Hb +
                Platelets +
                PT +
                PT+
                APTT +
                INR +
                Fibrinogen +
                PostOpFluids
               ,family=binomial(link='logit') , data=DataForLogistic))
summary(model)$coefficients[summary(model)$coefficients[,4] < 0.1,4]
summary( model )


model <- glm(formula = AFLogical ~ (ProcDetails == 'CABG') + Age + Gender + AdditiveEUROScore + 
                 SCTSLogisticEuroSCORE + PreOpNa + PreopUrea + PreopCreat + 
                 PreOpAlb + PreopPLT + PreopAPTT + CPB + CRP + Mg + WBC + 
                 Hb + Platelets, family = binomial(link = "logit"), data = DataForLogistic)
summary( model )

LogisticProbility <- predict(object= model , newdata = DataForLogistic , type = c( "response") )

for(i in 1:dim(LogisticProbility)[1]){
  DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
  model <- (glm(formula = AFLogical ~ 
                  (ProcDetails == 'CABG' ) + Age + Gender + AdditiveEUROScore + 
                  SCTSLogisticEuroSCORE + PreOpNa + PreopUrea + PreopCreat + 
                  PreOpAlb + PreopPLT + PreopAPTT + CPB + CRP + Mg + WBC + 
                  Hb + Platelets
                 , family=binomial(link='logit') , data=DataForLogistic[-i , ]))
  LogisticProbility[i,] <- predict(object= model , newdata = DataForLogistic[i , ] , type = c( "response") )
}

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$AFLogical == 1))
BC_CalculateAreaUnderCurve(PerformanceSweep)

x11()
BC_PlotCompareSingleHists(LogisticProbility[MasterData$AFLogical ==0] , LogisticProbility[MasterData$AFLogical ==1] , main ='Histogram of Logistic Model Output'  )
t.test(LogisticProbility[MasterData$AFLogical ==0] ,  LogisticProbility[MasterData$AFLogical ==1])

x11()
plot(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$AFLogical ==0] , seq(0,1,0.01) ) , col = 'blue' , type ='l' , xlab = 'Output Logistic Model' , ylab = 'Culmulative Prob')
lines(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$AFLogical ==1] , seq(0,1,0.01) ) , col = 'red' , type ='l')
abline(v = 0.2)
title('CDF Analysis')

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$AFLogical == 1))
ROCplot2 <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot2 <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept =( 1-0.2)) + geom_hline(yintercept = 0.2)+  ggtitle('PPV vs NPV Curves')

x11(20,14)
grid.arrange(ROCplot2 , NPVPPVPlot2 + xlim(( 1-0.2) , 1) , nrow= 2)

x11(20,14)
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$AFLogical == 1), BinWidth = 0.05))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))

#### Bayes Linear Bayes ####

{ 
DataForBayes <- POM_SampledImputation(MasterData)
DataForBayes <- DataForBayes[ , names(MasterData) %in% c('ProcDetails' , 
                                                         'Age' , 
                                                         'Gender',
                                                         'AdditiveEUROScore' , 
                                                         'SCTSLogisticEuroSCORE',
                                                         'PreOpNa' , 
                                                         'PreopUrea' , 
                                                         'PreopCreat' , 
                                                         'PreOpAlb' , 
                                                         'PreopPLT' ,
                                                         'PreopAPTT',
                                                         'CPB' ,
                                                         'CRP' ,
                                                         'Mg',
                                                         'Platelets',
                                                         'WBC',
                                                         'Hb')]
DataForBayes$ProcDetails <-DataForBayes$ProcDetails == 'CABG'
DataForBayes$Gender <- DataForBayes$Gender  == 'Female'
AFLogical <- as.numeric(MasterPreOpData$AFLogical)==2

}


PosteriorProbability <- BLBC_FitBayesLinearBayesClassifier(Data = as.matrix(DP_RemoveNaRows(DataForBayes) ) , Labels = AFLogical[DP_FindNARows(DataForBayes)] )

x11()
BC_PlotCompareSingleHists(PosteriorProbability[AFLogical ==0] ,PosteriorProbability[AFLogical ==1] , main ='Histogram of Logistic Model Output'  , breaks = 20)
t.test(PosteriorProbability[AFLogical ==0] , PosteriorProbability[AFLogical ==1])

x11()
plot(seq(0,1,0.01) , FM_CalculateCDFS(PosteriorProbability[AFLogical==0] , seq(0,1,0.01) ) , col = 'blue' , type ='l' , xlab = 'Output Logistic Model' , ylab = 'Culmulative Prob')
lines(seq(0,1,0.01) , FM_CalculateCDFS(PosteriorProbability[AFLogical ==1] , seq(0,1,0.01) ) , col = 'red' , type ='l')
abline(v = 0.2)
title('CDF Analysis')

PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability ,AFLogical[DP_FindNARows(DataForBayes)] == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')

BC_CalculateAreaUnderCurve(PerformanceSweep3)

x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot + xlim(( 1-0.2) , 1) , nrow= 2)

x11(20,14)
EmpericalProbabilityStructure3 <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , MasterData$AFLogical == 1), BinWidth = 0.05))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure3) + ggtitle('Emperical Probabilities'))

plot(1- PerformanceSweep2[,2] , PerformanceSweep2[,1] , type ='l' , col ='blue')
lines(1- PerformanceSweep[,2] , PerformanceSweep[,1] , type ='l' , col = 'red')
lines(1- PerformanceSweep3[,2] , PerformanceSweep3[,1] , type ='l' , col = 'green')
abline(0,1)
abline(v = 0.3)
abline(h = 0.6)


#### Gaussian Process Classfier #####

library('kernlab')

GaussianLogistic <- LogisticProbility
for(i in 1:dim(DataForBayes)[1]){  
  
  DataForBayes <- POM_SampledImputation(MasterData)
  DataForBayes <- DataForBayes[ , names(MasterData) %in% c('ProcDetails' , 
                                                           'Age' , 
                                                           'Gender',
                                                           'AdditiveEUROScore' , 
                                                           'SCTSLogisticEuroSCORE',
                                                           'PreOpNa' , 
                                                           'PreopUrea' , 
                                                           'PreopCreat' , 
                                                           'PreOpAlb' , 
                                                           'PreopPLT' ,
                                                           'PreopAPTT',
                                                           'CPB' ,
                                                           'CRP' ,
                                                           'Mg',
                                                           'Platelets',
                                                           'WBC',
                                                           'Hb')]
  DataForBayes$ProcDetails <-DataForBayes$ProcDetails == 'CABG'
  DataForBayes$Gender <- DataForBayes$Gender  == 'Female'
  AFLogical <- as.numeric(MasterPreOpData$AFLogical)==2
  
foo <- gausspr(DataForBayes[-i,], as.factor(AFLogical[-i]) )
GaussianLogistic[i] <- predict(foo , DataForBayes[i,], type = c( "probabilities"))[,2]
DP_WaitBar(i/dim(DataForBayes)[1])
}

foo <- gausspr(DataForBayes, as.factor(AFLogical) )
GaussianLogistic2 <- predict(foo , DataForBayes, type = c( "probabilities"))[,2]


x11()
BC_PlotCompareSingleHists(GaussianLogistic2[MasterData$AFLogical ==0] , GaussianLogistic2[MasterData$AFLogical ==1] , main ='Histogram of Logistic Model Output'  )
t.test(GaussianLogistic2[MasterData$AFLogical ==0] ,  GaussianLogistic2[MasterData$AFLogical ==1])

PerformanceSweep4 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(GaussianLogistic2 ,AFLogical[DP_FindNARows(DataForBayes)] == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep4) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep4) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')

BC_CalculateAreaUnderCurve(PerformanceSweep4)


x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot + xlim(( 1-0.2) , 1) , nrow= 2)



x11(20,14)
EmpericalProbabilityStructure3 <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(GaussianLogistic , MasterData$AFLogical == 1), BinWidth = 0.05))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure3) + ggtitle('Emperical Probabilities'))

plot(1- PerformanceSweep2[,2] , PerformanceSweep2[,1] , type ='l' , col ='blue')
lines(1- PerformanceSweep[,2] , PerformanceSweep[,1] , type ='l' , col = 'red')
lines(1- PerformanceSweep3[,2] , PerformanceSweep3[,1] , type ='l' , col = 'green')
lines(1- PerformanceSweep4[,2] , PerformanceSweep4[,1] , type ='l' , col = 'purple')

abline(0,1)
abline(v = 0.5)
