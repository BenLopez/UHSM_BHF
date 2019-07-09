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
  
  PostOpFlow <- POM_ExtractPostOpFromFlow(FlowIndex2017)
  PostOpFlow <- PostOpFlow[ which(rownames(PostOpFlow) %in% MasterPreOpData$NewPseudoId) ,  ]
  PostOpFlow <- PostOpFlow[order(rownames(PostOpFlow)) , ]

  PostOpRiskScore <- POM_ExtraPostOpScores(PostOpScoresIndex2017 = PostOpScoresIndex2017)
  PostOpRiskScore <- PostOpRiskScore[ which(rownames(PostOpRiskScore) %in% MasterPreOpData$PseudoId ) ,  ]
  PostOpRiskScore <- PostOpRiskScore[ order( rownames( PostOpRiskScore ) ) , ]
  
  TreatmentsLogcal <- POM_TreatmentBeforeAFib(FluidsIndex2017 )
  TreatmentsLogcal <- TreatmentsLogcal[ which(rownames(TreatmentsLogcal) %in% MasterPreOpData$NewPseudoId ) ,  ]
  TreatmentsLogcal <- TreatmentsLogcal[ order( rownames( TreatmentsLogcal ) ) , ]
  
  }

MasterData <- cbind(MasterPreOpData , PostOpBioChem , PostOpHaem , PostOpFluids , PostOpFlow , TreatmentsLogcal )
MasterData <- MasterData[which(MasterData$PseudoId %in% rownames( PostOpRiskScore )) , ]
MasterData <- MasterData[order(MasterData$PseudoId) , ]
MasterData <- cbind(MasterData , PostOpRiskScore)

MasterData$ProcDetails <- DP_AssignSurgeryLabels(MasterData$ProcDetails)

MasterData$BMI = MasterData$Weight/(MasterData$Height/100)

MasterData$NYHAGrade[MasterData$NYHAGrade == "Not Known"] = NA 
MasterData$NYHAGrade[MasterData$NYHAGrade == ""] = NA 

MasterData$AnginaGrade[MasterData$AnginaGrade == ""] = NA 
MasterData$AnginaGrade[MasterData$AnginaGrade == "Unknown"] = NA 

MasterData$EjectionFractionCategory[MasterData$EjectionFractionCategory == ""] = NA 
MasterData$EjectionFractionCategory[MasterData$EjectionFractionCategory == "Not measured"] = NA 

# Transformations
MasterData$LogisticEUROScore <- log(MasterData$LogisticEUROScore)
MasterData$SCTSLogisticEuroSCORE <- MasterData$SCTSLogisticEuroSCORE

MasterData$CRP <- log(MasterData$CRP)
MasterData$PreOpCRP <- log(MasterData$PreOpCRP)

MasterData$RACE <- log(MasterData$RACE)
MasterData$logCASUS <- log(MasterData$logCASUS)

MasterData$PreopBili <- log(MasterData$PreopBili)
MasterData$Bilirubin <- log(MasterData$Bilirubin)

{ #biochem
  MasterData$dMg <- MasterData$Mg - MasterData$PreopMg
  MasterData$dAlb <- MasterData$Albumin - MasterData$PreOpAlb
  MasterData$dUrea <- MasterData$Urea - MasterData$PreopUrea
  MasterData$dNa <- MasterData$Na - MasterData$PreOpNa
  MasterData$dBili <- MasterData$Bilirubin - MasterData$PreopBili
  MasterData$dCreat <- MasterData$Creatinine - MasterData$PreopCreat
  MasterData$dK <-   MasterData$K - MasterData$PreOpK
  MasterData$dPT <-  MasterData$PT - MasterData$PreopPT
  #bloods
  MasterData$dCRP <- MasterData$CRP - MasterData$PreOpCRP 
  MasterData$dPLT <- MasterData$Platelets - MasterData$PreopPLT
  MasterData$dHb <- MasterData$Hb - MasterData$PreopHb 
  MasterData$dAPTT <- MasterData$APTT - MasterData$PreopAPTT
  MasterData$dWBC <- MasterData$WBC - MasterData$PreopWBC
}

MasterData <- POM_RemoveOutliers(MasterPreOpData = MasterData)

MasterData$ProcDetails <- as.factor(MasterData$ProcDetails)

PerentageMissing <- (apply(MasterData , 2 , function(X){sum(is.na (X)) } )/apply(MasterData , 2 , function(X){length(X) } ))*100

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

# Post op meds
MasterData$Adrenaline <- as.factor(MasterData$Adrenaline)
(summary(MasterData$Adrenaline[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
(summary(MasterData$Adrenaline[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
prop.test( c(summary(MasterData$Adrenaline[MasterData$AFLogical == 1])[1] ,summary(MasterData$Adrenaline[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) ))

MasterData$Noradrenaline <- as.factor(MasterData$Noradrenaline)
(summary(MasterData$Noradrenaline[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
(summary(MasterData$Noradrenaline[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
prop.test( c(summary(MasterData$Noradrenaline[MasterData$AFLogical == 1])[1] ,summary(MasterData$Noradrenaline[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) ))

MasterData$Dopamine <- as.factor(MasterData$Dopamine)
(summary(MasterData$Dopamine[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
(summary(MasterData$Dopamine[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
prop.test( c(summary(MasterData$Dopamine[MasterData$AFLogical == 1])[1] ,summary(MasterData$Dopamine[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) ))

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

ProcDetails <- as.factor(MasterData$EjectionFractionCategory)

for(i in 1:length(levels(ProcDetails))){
  print(levels(ProcDetails)[i])
  tmp <- as.factor(MasterData$EjectionFractionCategory == levels(ProcDetails)[i] )
  print(summary(tmp[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
  print(summary(tmp[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
  print(prop.test( c(summary(tmp[MasterData$AFLogical == 1])[1] ,summary(tmp[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) )))
}
# NYHA Grade

ProcDetails <- as.factor(MasterData$NYHAGrade)

for(i in 1:length(levels(ProcDetails))){
  print(levels(ProcDetails)[i])
  print(sum(MasterData$NYHAGrade == levels(ProcDetails)[i] , na.rm = T))
  print(sum(MasterData$NYHAGrade[MasterData$AFLogical ==1] == levels(ProcDetails)[i] , na.rm = T)/sum(MasterData$NYHAGrade == levels(ProcDetails)[i] , na.rm = T))
  tmp <- as.factor(MasterData$NYHAGrade == levels(ProcDetails)[i] )
  print(summary(tmp[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
  print(summary(tmp[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
  print(prop.test( c(summary(tmp[MasterData$AFLogical == 1])[1] ,summary(tmp[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) )))
}

# Angina Grade

ProcDetails <- as.factor(MasterData$AnginaGrade)

for(i in 1:length(levels(ProcDetails))){
  print(levels(ProcDetails)[i])
  print(sum(MasterData$AnginaGrade == levels(ProcDetails)[i] , na.rm = T))
  print(sum(MasterData$AnginaGrade[MasterData$AFLogical ==1] == levels(ProcDetails)[i] , na.rm = T)/sum(MasterData$AnginaGrade == levels(ProcDetails)[i] , na.rm = T))
  tmp <- as.factor(MasterData$AnginaGrade == levels(ProcDetails)[i] )
  print(summary(tmp[MasterData$AFLogical == 1])/sum(MasterData$AFLogical == 1))*100
  print(summary(tmp[MasterData$AFLogical == 0])/sum(MasterData$AFLogical == 0))*100
  print(prop.test( c(summary(tmp[MasterData$AFLogical == 1])[1] ,summary(tmp[MasterData$AFLogical == 0])[1]) , c(  sum(MasterData$AFLogical == 1), sum(MasterData$AFLogical == 0) )))
}

# Urgency 
ProcDetails <- as.factor(MasterData$Urgency)

for(i in 1:length(levels(ProcDetails))){
  print(levels(ProcDetails)[i])
  print(sum(MasterData$Urgency == levels(ProcDetails)[i] , na.rm = T))
  print(sum(MasterData$Urgency[MasterData$AFLogical ==1] == levels(ProcDetails)[i] , na.rm = T)/sum(MasterData$Urgency == levels(ProcDetails)[i] , na.rm = T))
  tmp <- as.factor(MasterData$Urgency == levels(ProcDetails)[i] )
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

NamesofPreOperativePredictors <- c('Age' ,
                                   'Gender' ,
                                   'Weight' , 
                                   'ProcDetails' ,
                                   'BMI',
                                   "LogisticEUROScore" ,
                                   "EjectionFractionCategory",
                                   "NYHAGrade",
                                   "AnginaGrade",
                                   "Urgency",
                                   'PreOpNa',
                                   'PreOpK',
                                   'PreopMg',
                                   'PreopUrea' ,
                                   'PreOpCRP',
                                   'PreOpAlb' ,
                                   'PreopBili' ,
                                   'PreopCreat',
                                   "PreopHb" ,
                                   "PreopPLT" ,
                                   "PreopWBC",
                                   "PreopPT",
                                   "PreopAPTT")
PreoperativeIndices = which( names(MasterData) %in% NamesofPreOperativePredictors)


counter <- 1
tteststruct <- matrix(0 , length( PreoperativeIndices ) , 1)
for(i in PreoperativeIndices){
if(!is.numeric(MasterData[MasterData$AFLogical == 0 , i])){
  next
}  
BC_PlotCompareSingleHists(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i] , main =  names(MasterData)[i])
print(names(MasterData)[i])
print(t.test(MasterData[MasterData$AFLogical == 1 , i] , MasterData[MasterData$AFLogical == 0 , i]))
print(sqrt(var(MasterData[MasterData$AFLogical == 1 , i] , na.rm = T)))
print(sqrt(var(MasterData[MasterData$AFLogical == 0 , i] , na.rm = T)))
tteststruct[counter] <- (t.test(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i])$p.value)
counter <- counter + 1
}

# Logistic regression analysis individual variables.

DataForLogistic <- data.frame(MasterData[  , c(1,PreoperativeIndices)])
for(i in 2:dim(DataForLogistic)[2]){
  model <- (glm(formula = as.formula(paste0("AFLogical ~" , names(DataForLogistic)[i] , '+ LogisticEUROScore + ProcDetails') )
                ,family=binomial(link='logit') , data=DataForLogistic))
# odds ratio
# print( names(DataForLogistic)[i] )
  print(exp(coef(model)[2]))
  print( exp(confint(model)[2,] ))  
  print( summary( model )$coefficients )
  #print(auc(model$y, predict(model , DataForLogistic ,  type = c( "response")) ) )
}

MasterData$ProcDetails <- as.factor(MasterData$ProcDetails)


StepResults <- matrix(0 , 100, length(NamesofPreOperativePredictors) )
colnames(StepResults) <- NamesofPreOperativePredictors

for(ii in 1:100){
#imputed_Data <- mice(MasterData[ , c(1,PreoperativeIndices)], m=1)
#DataForLogistic <- MasterData[ , c(1,PreoperativeIndices)]
DataForLogistic <- POM_SampledImputation(MasterData[ , c(1,PreoperativeIndices)] )
#  POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,PreoperativeIndices)]))
#DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,PreoperativeIndices)]))
model <- (glm(formula = AFLogical ~ 
                ProcDetails  +
                Age +
                Gender +
                Weight +
                BMI +
                Urgency +
                AnginaGrade +
                NYHAGrade +
                EjectionFractionCategory+
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
                PreopAPTT
                ,family=binomial(link='logit') , data=DataForLogistic))
summary(model)$coefficients[summary(model)$coefficients[,4] < 0.1,4]
stepoutput <- stepAIC(model)

StepResults[ii,] <- apply(as.matrix(NamesofPreOperativePredictors) , 1 , function(X){grepl(X ,as.character(stepoutput$formula)[3])} )
}

indexesfromstep <- which(names(MasterData) %in% colnames(StepResults[ , apply(StepResults , 2 , function(X){ (sum(X)/length(X))*100 })>50]))

model <- glm(formula = FB_CreateFormula('AFLogical' , indexesfromstep , MasterData), 
             family = binomial(link = "logit"), data = DataForLogistic)
summary(model)
LogisticProbility <- predict(model , DataForLogistic , type = c( "response"))

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$AFLogical == 1))
BC_CalculateAreaUnderCurve(PerformanceSweep)

LogisticProbility <- matrix(0,length(MasterData$AFLogical),1)
for(i in 1:dim(LogisticProbility)[1]){
  DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,PreoperativeIndices)]))
  model <- glm(formula = FB_CreateFormula('AFLogical' , indexesfromstep , MasterData) , family = binomial(link = "logit"), data = DataForLogistic[-i,])
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

BC_CalculateAreaUnderCurve(PerformanceSweep)

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

NamesofPostOperativePredictors <- c('CPB',
                                   'Na', 
                                   'K',
                                   'Urea',
                                   'Creatinine',
                                   'CRP',
                                   'Albumin',
                                   'Bilirubin',
                                   'Mg',
                                   'WBC',
                                   'Hb',
                                   'Platelets',
                                   'PT',
                                   'APTT',
                                   'INR',
                                   'Fibrinogen',
                                   'PostOpFluids',
                                   'HR',
                                   "ReliableART.S",
                                   "ReliableART.M",
                                   "ReliableART.D",
                                   "CVP",
                                   "SpO2",
                                   'dNa',
                                   'dK',
                                   'dMg',
                                   'dUrea' ,
                                   'dCRP',
                                   'dAlb' ,
                                   'dBili' ,
                                   'dCreat',
                                   "dHb" ,
                                   "dPLT" ,
                                   "dWBC",
                                   "dPT",
                                   "dAPTT",
                                   "SOFA",
                                   "RACE",
                                   'logCASUS',
                                   'Noradrenaline',
                                   'Dopamine',
                                   'Adrenaline')
PostOpIndexes <- which(names(MasterData) %in% NamesofPostOperativePredictors)

NamesofPreOperativePredictors <- names(MasterData)[PostOpIndexes ]
counter <- 1
tteststruct <- matrix(0 , length( PostOpIndexes[2:length(PostOpIndexes)]) , 1)
for(i in PostOpIndexes[2:length(PostOpIndexes)]){
  if(!is.numeric(MasterData[MasterData$AFLogical == 0 , i])){next}
  BC_PlotCompareSingleHists(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i] , main =  names(MasterData)[i])
  
  print(names(MasterData)[i])
  print(t.test(MasterData[MasterData$AFLogical == 1 , i] , MasterData[MasterData$AFLogical == 0 , i]))
  print(sqrt(var(MasterData[MasterData$AFLogical == 1 , i] , na.rm = T)))
  print(sqrt(var(MasterData[MasterData$AFLogical == 0 , i] , na.rm = T)))
  tteststruct[counter] <- (t.test(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i])$p.value)
  counter <- counter + 1
}

# difference transforms
{
MasterData$dUrea <- abs(DP_NormaliseData(MasterData$dUrea) )
MasterData$dK <-abs(DP_NormaliseData(MasterData$dK))
MasterData$dCreat <- abs(DP_NormaliseData(MasterData$dCreat) )
MasterData$dNa <- abs(DP_NormaliseData(MasterData$dNa) )
MasterData$dCRP <- abs(DP_NormaliseData(abs(DP_NormaliseData(MasterData$dCRP) )))
MasterData$dWBC <- abs(DP_NormaliseData(MasterData$dWBC))
MasterData$dAlb <- abs(DP_NormaliseData(MasterData$dAlb) )
}

NamesofPreOperativePredictors <- names(MasterData)[PostOpIndexes ]
counter <- 1
tteststruct <- matrix(0 , length( PostOpIndexes[2:length(PostOpIndexes)]) , 1)
for(i in PostOpIndexes[2:length(PostOpIndexes)]){
  if(!is.numeric(MasterData[MasterData$AFLogical == 0 , i])){next}
  BC_PlotCompareSingleHists(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i] , main =  names(MasterData)[i])
  
    print(names(MasterData)[i])
  print(t.test(MasterData[MasterData$AFLogical == 1 , i] , MasterData[MasterData$AFLogical == 0 , i]))
  print(sqrt(var(MasterData[MasterData$AFLogical == 1 , i] , na.rm = T)))
  print(sqrt(var(MasterData[MasterData$AFLogical == 0 , i] , na.rm = T)))
  tteststruct[counter] <- (t.test(MasterData[MasterData$AFLogical == 0 , i] , MasterData[MasterData$AFLogical == 1 , i])$p.value)
  counter <- counter + 1
}

DataForLogistic <- data.frame(MasterData[  ,c(1, PostOpIndexes) ])

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
model <- (glm(formula = FB_CreateFormula('AFLogical' , c(PreoperativeIndices,PostOpIndexes) , MasterData)
               ,family=binomial(link='logit') , data=DataForLogistic))
summary( model )

StepResults2 <- matrix(0 , 100, length(c(PreoperativeIndices,PostOpIndexes)) )
colnames(StepResults2) <- names(MasterData[c(PreoperativeIndices,PostOpIndexes)])

for(ii in 1:100){
DataForLogistic <- POM_SampledImputation(data.frame(MasterData[  ,  ]))
model <- (glm(formula = FB_CreateFormula('AFLogical' , c(PreoperativeIndices,PostOpIndexes) , MasterData)
              ,family=binomial(link='logit') , data=DataForLogistic))
stepoutput <- stepAIC(model)
StepResults2[ii,] <- apply(as.matrix(names(MasterData[c(PreoperativeIndices,PostOpIndexes)])) , 1 , function(X){grepl(X ,as.character(stepoutput$formula)[3])} )
}

indexesfromstep <- which(names(MasterData) %in% colnames(StepResults2[ , apply(StepResults2 , 2 , function(X){ (sum(X)/length(X))*100 })>50]))

tmp <- (cov2cor(cov(StepResults2))[apply(StepResults2 , 2 , function(X){ (sum(X)/length(X))*100 })>50 , ])

LogisticProbility <- matrix(0 , dim(DataForLogistic)[1],1)
for(i in 1:dim(LogisticProbility)[1]){
  DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
  model <- glm(formula = AFLogical ~ 
                 Age + 
                 Weight +
                 Urgency + 
                 LogisticEUROScore +
                 PreopUrea +
                 PreopCreat+
                 PreOpAlb +
                 PreopPLT +
                 PreopAPTT+
                 CPB +
                 Na +
                 Urea +
                 CRP +
                 Mg +
                 Hb +
                 Platelets +
                 PT +
                 APTT +
                 ReliableART.D +
                 SOFA +
                 dBili, family = binomial(link = "logit"),  data = DataForLogistic[-i,])
  LogisticProbility[i,] <- predict(object= model , newdata = DataForLogistic[i , ] , type = c( "response") )
}

summary(model)

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
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$AFLogical == 1), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))

#### Bayes Linear Bayes ####

{
DataForBayes <- POM_SampledImputation(MasterData)
DataForBayes <- DataForBayes[ , names(MasterData) %in% c('Age' , 
                                                         'Weight',
                                                         'Urgency',
                                                         'LogisticEUROScore',
                                                         'PreopUrea' , 
                                                         'PreopCreat' , 
                                                         'PreOpAlb' , 
                                                         'PreopPLT' ,
                                                         'PreopAPTT',
                                                         'CPB' ,
                                                         'CRP',
                                                         'Na',
                                                         'Urea',
                                                         'CRP' ,
                                                         'Mg',
                                                         'Hb',
                                                         'Platelets',
                                                         'PT',
                                                         'APTT',
                                                         'ReliableART.D',
                                                         'SOFA',
                                                         'dBili')]
DataForBayes$Urgency <-DataForBayes$Urgency == 'Urgent'
AFLogical <- as.numeric(MasterData$AFLogical)==2
}


PosteriorProbability <- BLBC_FitBayesLinearBayesClassifier(Data = as.matrix(DataForBayes) , Labels = MasterData$AFLogical ==1 )

x11()
BC_PlotCompareSingleHists(PosteriorProbability[AFLogical ==0] ,PosteriorProbability[AFLogical ==1] , main ='Histogram of Logistic Model Output'  , breaks = 20)
t.test(PosteriorProbability[AFLogical ==0] , PosteriorProbability[AFLogical ==1])

x11()
plot(seq(0,1,0.01) , FM_CalculateCDFS(PosteriorProbability[AFLogical==0] , seq(0,1,0.01) ) , col = 'blue' , type ='l' , xlab = 'Output Logistic Model' , ylab = 'Culmulative Prob')
lines(seq(0,1,0.01) , FM_CalculateCDFS(PosteriorProbability[AFLogical ==1] , seq(0,1,0.01) ) , col = 'red' , type ='l')
abline(v = 0.2)
title('CDF Analysis')

PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability ,AFLogical[DataForBayes] == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')

BC_CalculateAreaUnderCurve(PerformanceSweep3)

x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot + xlim(( 1-0.2) , 1) , nrow= 2)

x11(20,14)
EmpericalProbabilityStructure3 <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , MasterData$AFLogical == 1), BinWidth = 0.05))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure3) + ggtitle('Emperical Probabilities'))

###### AUC Selection Criteria #####

FC_StepWiseForwardAUC( c(PreoperativeIndices) , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
model <- glm(formula = AFLogical ~ 
               LogisticEUROScore+
               PreopUrea +
               PreOpAlb +
               (Urgency == 'Urgent') +
               PreopAPTT +
               PreopCreat +
               Weight +
               PreopBili +
               PreopHb +
               PreOpK +
               PreOpCRP
               , family = binomial(link = "logit"),  data = DataForLogistic)
summary(model)

LogisticProbility <- predict(model , DataForLogistic , type = c('response'))

PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,DataForLogistic$AFLogical == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')

BC_CalculateAreaUnderCurve(PerformanceSweep3)

FC_StepWiseForwardAUC( c(PostOpIndexes) , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
model <- glm(formula = AFLogical ~ 
               logCASUS +
               Urea +
               CPB +
               CRP +
               Mg +
               dPLT +
               Noradrenaline +
               Creatinine +
               Hb +
               dPT +
               ReliableART.M
             , family = binomial(link = "logit"),  data = DataForLogistic)
summary(model)

LogisticProbility <- predict(model , DataForLogistic , type = c('response'))

PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,DataForLogistic$AFLogical == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')

BC_CalculateAreaUnderCurve(PerformanceSweep3)

x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot + xlim(( 1-0.2) , 1) , nrow= 2)

FC_StepWiseForwardAUC( c(PreoperativeIndices , PostOpIndexes ) , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
model <- glm(formula = AFLogical ~ 
               LogisticEUROScore+
               logCASUS+
               Urea +
               dK +
               Albumin +
               PreopCreat +
               PreopUrea+
               Weight +
               PreopBili +
               PreopHb +
               PreOpK +
               PreOpCRP +
               Platelets +
               Mg +
               dNa +
               ReliableART.M +
               Gender +
               Urgency +
               Hb +
               PreOpAlb +
               dBili +
               APTT +
               CPB +
               PreOpCRP +
               WBC
             , family = binomial(link = "logit"),  data = DataForLogistic)
summary(model)

LogisticProbility <- predict(model , DataForLogistic , type = c('response'))

PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,DataForLogistic$AFLogical == 1))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')

BC_CalculateAreaUnderCurve(PerformanceSweep3)

x11(20,14)
grid.arrange(ROCplot , NPVPPVPlot + xlim(( 1-0.2) , 1) , nrow= 2)

