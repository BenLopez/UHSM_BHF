
# Patient demographic infomortion
{
  {
    sum(MasterData$MasterData$AFlogical == 1)
    sum(MasterData$MasterData$AFlogical == 0)
    
    # Age
    t.test(MasterData$Age[MasterData$MasterData$AFlogical == 1] , MasterData$Age[MasterData$MasterData$AFlogical == 0])
    sqrt(var(MasterData$Age[MasterData$MasterData$AFlogical == 1]))
    sqrt(var(MasterData$Age[MasterData$MasterData$AFlogical == 0]))
    
    # Gender
    (summary(MasterData$Gender[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    (summary(MasterData$Gender[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    prop.test( c(summary(MasterData$Gender[MasterData$MasterData$AFlogical == 1])[1] ,summary(MasterData$Gender[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))
    
    # Post op meds
    MasterData$Adrenaline <- as.factor(MasterData$Adrenaline)
    (summary(MasterData$Adrenaline[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    (summary(MasterData$Adrenaline[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    prop.test( c(summary(MasterData$Adrenaline[MasterData$MasterData$AFlogical == 1])[1] ,summary(MasterData$Adrenaline[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))
    
    MasterData$Noradrenaline <- as.factor(MasterData$Noradrenaline)
    (summary(MasterData$Noradrenaline[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    (summary(MasterData$Noradrenaline[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    prop.test( c(summary(MasterData$Noradrenaline[MasterData$MasterData$AFlogical == 1])[1] ,summary(MasterData$Noradrenaline[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))
    
    MasterData$Dopamine <- as.factor(MasterData$Dopamine)
    (summary(MasterData$Dopamine[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    (summary(MasterData$Dopamine[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    prop.test( c(summary(MasterData$Dopamine[MasterData$MasterData$AFlogical == 1])[1] ,summary(MasterData$Dopamine[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))
    
    # Weight
    t.test(MasterData$Weight[MasterData$MasterData$AFlogical == 1] , MasterData$Weight[MasterData$MasterData$AFlogical == 0])
    sqrt(var(MasterData$Weight[MasterData$MasterData$AFlogical == 1]))
    sqrt(var(MasterData$Weight[MasterData$MasterData$AFlogical == 0]))
    
    # BMI
    
    t.test(MasterData$BMI[MasterData$MasterData$AFlogical == 1] , MasterData$BMI[MasterData$MasterData$AFlogical == 0])
    sqrt(var(MasterData$BMI[MasterData$MasterData$AFlogical == 1] , na.rm =T))
    sqrt(var(MasterData$BMI[MasterData$MasterData$AFlogical == 0], na.rm =T))
    
    
    # Pre operative heart rhythm
    tmp <- as.factor(MasterData$Pre_OperativeHeartRhythm == 'Sinus Rhythm')
    (summary(tmp[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    (summary(tmp[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    prop.test( c(summary(tmp[MasterData$MasterData$AFlogical == 1])[1] ,summary(tmp[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))
    
    
    # Pre operative RRT
    
    tmp <- as.factor(MasterData$PreopRRT == 'Yes')
    (summary(tmp[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    (summary(tmp[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    prop.test( c(summary(tmp[MasterData$MasterData$AFlogical == 1])[1] ,summary(tmp[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))
  }
  
  # Surgery 
  ProcDetails <- as.factor(MasterData$ProcDetails)
  
  for(i in 1:length(levels(ProcDetails))){
    print(levels(ProcDetails)[i])
    tmp <- as.factor(MasterData$ProcDetails == levels(ProcDetails)[i] )
    print(summary(tmp[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    print(summary(tmp[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    print(prop.test( c(summary(tmp[MasterData$MasterData$AFlogical == 1])[1] ,summary(tmp[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) )))
  }
  
  ProcDetails <- as.factor(MasterData$EjectionFractionCategory)
  
  for(i in 1:length(levels(ProcDetails))){
    print(levels(ProcDetails)[i])
    tmp <- as.factor(MasterData$EjectionFractionCategory == levels(ProcDetails)[i] )
    print(summary(tmp[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    print(summary(tmp[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    print(prop.test( c(summary(tmp[MasterData$MasterData$AFlogical == 1])[1] ,summary(tmp[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) )))
  }
  # NYHA Grade
  
  ProcDetails <- as.factor(MasterData$NYHAGrade)
  
  for(i in 1:length(levels(ProcDetails))){
    print(levels(ProcDetails)[i])
    print(sum(MasterData$NYHAGrade == levels(ProcDetails)[i] , na.rm = T))
    print(sum(MasterData$NYHAGrade[MasterData$MasterData$AFlogical ==1] == levels(ProcDetails)[i] , na.rm = T)/sum(MasterData$NYHAGrade == levels(ProcDetails)[i] , na.rm = T))
    tmp <- as.factor(MasterData$NYHAGrade == levels(ProcDetails)[i] )
    print(summary(tmp[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    print(summary(tmp[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    print(prop.test( c(summary(tmp[MasterData$MasterData$AFlogical == 1])[1] ,summary(tmp[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) )))
  }
  
  # Angina Grade
  
  ProcDetails <- as.factor(MasterData$AnginaGrade)
  
  for(i in 1:length(levels(ProcDetails))){
    print(levels(ProcDetails)[i])
    print(sum(MasterData$AnginaGrade == levels(ProcDetails)[i] , na.rm = T))
    print(sum(MasterData$AnginaGrade[MasterData$MasterData$AFlogical ==1] == levels(ProcDetails)[i] , na.rm = T)/sum(MasterData$AnginaGrade == levels(ProcDetails)[i] , na.rm = T))
    tmp <- as.factor(MasterData$AnginaGrade == levels(ProcDetails)[i] )
    print(summary(tmp[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    print(summary(tmp[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    print(prop.test( c(summary(tmp[MasterData$MasterData$AFlogical == 1])[1] ,summary(tmp[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) )))
  }
  
  # Urgency 
  ProcDetails <- as.factor(MasterData$Urgency)
  
  for(i in 1:length(levels(ProcDetails))){
    print(levels(ProcDetails)[i])
    print(sum(MasterData$Urgency == levels(ProcDetails)[i] , na.rm = T))
    print(sum(MasterData$Urgency[MasterData$MasterData$AFlogical ==1] == levels(ProcDetails)[i] , na.rm = T)/sum(MasterData$Urgency == levels(ProcDetails)[i] , na.rm = T))
    tmp <- as.factor(MasterData$Urgency == levels(ProcDetails)[i] )
    print(summary(tmp[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
    print(summary(tmp[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
    print(prop.test( c(summary(tmp[MasterData$MasterData$AFlogical == 1])[1] ,summary(tmp[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) )))
  }
  
}


# Surgery type boxplots
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
  if(!is.numeric(MasterData[MasterData$MasterData$AFlogical == 0 , i])){
    next
  }  
  BC_PlotCompareSingleHists(MasterData[MasterData$MasterData$AFlogical == 0 , i] , MasterData[MasterData$MasterData$AFlogical == 1 , i] , main =  names(MasterData)[i])
  print(names(MasterData)[i])
  print(t.test(MasterData[MasterData$MasterData$AFlogical == 1 , i] , MasterData[MasterData$MasterData$AFlogical == 0 , i]))
  print(sqrt(var(MasterData[MasterData$MasterData$AFlogical == 1 , i] , na.rm = T)))
  print(sqrt(var(MasterData[MasterData$MasterData$AFlogical == 0 , i] , na.rm = T)))
  tteststruct[counter] <- (t.test(MasterData[MasterData$MasterData$AFlogical == 0 , i] , MasterData[MasterData$MasterData$AFlogical == 1 , i])$p.value)
  counter <- counter + 1
}

# Logistic regression analysis individual variables.

DataForLogistic <- data.frame(MasterData[  , c(1,PreoperativeIndices)])
for(i in 2:dim(DataForLogistic)[2]){
  model <- (glm(formula = as.formula(paste0("MasterData$AFlogical ~" , names(DataForLogistic)[i] , '+ LogisticEUROScore + ProcDetails') )
                ,family=binomial(link='logit') , data=DataForLogistic))
  # odds ratio
  # print( names(DataForLogistic)[i] )
  print(exp(coef(model)[2]))
  print( exp(confint(model)[2,] ))  
  print( summary( model )$coefficients )
  #print(auc(model$y, predict(model , DataForLogistic ,  type = c( "response")) ) )
}

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + CRP
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)
#

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + WBC + Hb + Platelets
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + Urea
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)


DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + PreopUrea
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + Creatinine
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + PreopCreat
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + PreOpK + Creatinine + Urea  + PreopCreat + PreopUrea + weight
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)


DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + PreopHb
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)


DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + Hb
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)


DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + CPB + Mg
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + PreOpAlb
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)

DataForLogistic <- POM_SampledImputation(MasterData)
model <- glm(formula =MasterData$AFLogical ~ ProcDetails + LogisticEUROScore + Platelets
             ,family=binomial(link='logit') , data=DataForLogistic)
summary(model)



DataForLogistic <- data.frame(MasterData[  , c(1,PreoperativeIndices)])


StepResults <- matrix(0 , 100, length(NamesofPreOperativePredictors) )
colnames(StepResults) <- NamesofPreOperativePredictors

for(ii in 1:100){
  #imputed_Data <- mice(MasterData[ , c(1,PreoperativeIndices)], m=1)
  #DataForLogistic <- MasterData[ , c(1,PreoperativeIndices)]
  DataForLogistic <- POM_SampledImputation(MasterData[ , c(1,PreoperativeIndices)] )
  #  POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,PreoperativeIndices)]))
  #DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,PreoperativeIndices)]))
  model <- (glm(formula =MasterData$AFLogical ~ 
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

model <- glm(formula = FB_CreateFormula('MasterData$AFlogical' , indexesfromstep , MasterData), 
             family = binomial(link = "logit"), data = DataForLogistic)
summary(model)
LogisticProbility <- predict(model , DataForLogistic , type = c( "response"))

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$MasterData$AFlogical == 1))
BC_CalculateAreaUnderCurve(PerformanceSweep)

LogisticProbility <- matrix(0,length(MasterData$MasterData$AFlogical),1)
for(i in 1:dim(LogisticProbility)[1]){
  DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  , c(1,PreoperativeIndices)]))
  model <- glm(formula = FB_CreateFormula('MasterData$AFlogical' , indexesfromstep , MasterData) , family = binomial(link = "logit"), data = DataForLogistic[-i,])
  LogisticProbility[i,] <- predict(object= model , newdata = DataForLogistic[i , ] , type = c( "response") )
}

x11()
BC_PlotCompareSingleHists(LogisticProbility[MasterData$MasterData$AFlogical ==0] , LogisticProbility[MasterData$MasterData$AFlogical ==1] , main ='Histogram of Logistic Model Output'  )
t.test(LogisticProbility[MasterData$MasterData$AFlogical ==0] ,  LogisticProbility[MasterData$MasterData$AFlogical ==1])

x11()
plot(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$MasterData$AFlogical ==0] , seq(0,1,0.01) ) , col = 'blue' , type ='l' , xlab = 'Output Logistic Model' , ylab = 'Culmulative Prob')
lines(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$MasterData$AFlogical ==1] , seq(0,1,0.01) ) , col = 'red' , type ='l')
abline(v = 0.2)
title('CDF Analysis')

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$MasterData$AFlogical == 1))
ROCplot2 <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot2 <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept =( 1-0.2)) + geom_hline(yintercept = 0.2)+  ggtitle('PPV vs NPV Curves')

BC_CalculateAreaUnderCurve(PerformanceSweep)

x11(20,14)
grid.arrange(ROCplot2 , NPVPPVPlot2 + xlim(( 1-0.2) , 1) , nrow= 2)

x11(20,14)
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$MasterData$AFlogical == 1), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))


###### PostOp Analysis ######

{NamesofPostOperativePredictors <- c('CPB',
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
                                     'Adrenaline',
                                     'ArtPO2',
                                     'Lac',
                                     'FiO26num',
                                     'IABP2',
                                     'PAR',
                                     'VentilatedInOpLogical',
                                     'GCSnum',
                                     'PaO2OverFiO2')

MasterData$PaO2OverFiO2 <- MasterData$ArtPO2/MasterData$FiO26num
MasterData$VentilatedInOpLogical <- as.factor(MasterData$VentilatedInOpLogical)
MasterData$IABP2 <- as.factor(MasterData$VentilatedInOpLogical)
MasterData$CVVHDialysis <- as.factor(MasterData$CVVHDialysis)
PostOpIndexes <- which(names(MasterData) %in% NamesofPostOperativePredictors)
}

NamesofPreOperativePredictors <- names(MasterData)[PostOpIndexes ]
counter <- 1
tteststruct <- matrix(0 , length( PostOpIndexes[2:length(PostOpIndexes)]) , 1)
for(i in PostOpIndexes[2:length(PostOpIndexes)]){
  if(!is.numeric(MasterData[MasterData$MasterData$AFlogical == 0 , i])){next}
  BC_PlotCompareSingleHists(MasterData[MasterData$MasterData$AFlogical == 0 , i] , MasterData[MasterData$MasterData$AFlogical == 1 , i] , main =  names(MasterData)[i])
  
  print(names(MasterData)[i])
  print(t.test(MasterData[MasterData$MasterData$AFlogical == 1 , i] , MasterData[MasterData$MasterData$AFlogical == 0 , i]))
  print(sqrt(var(MasterData[MasterData$MasterData$AFlogical == 1 , i] , na.rm = T)))
  print(sqrt(var(MasterData[MasterData$MasterData$AFlogical == 0 , i] , na.rm = T)))
  tteststruct[counter] <- (t.test(MasterData[MasterData$MasterData$AFlogical == 0 , i] , MasterData[MasterData$MasterData$AFlogical == 1 , i])$p.value)
  counter <- counter + 1
}

(summary(MasterData$VentilatedInOpLogical[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
(summary(MasterData$VentilatedInOpLogical[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
prop.test( c(summary(MasterData$VentilatedInOpLogical[MasterData$MasterData$AFlogical == 1])[1] ,summary(MasterData$VentilatedInOpLogical[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))

(summary(MasterData$CVVHDialysis[MasterData$MasterData$AFlogical == 1])/sum(MasterData$MasterData$AFlogical == 1))*100
(summary(MasterData$CVVHDialysis[MasterData$MasterData$AFlogical == 0])/sum(MasterData$MasterData$AFlogical == 0))*100
prop.test( c(summary(MasterData$CVVHDialysis[MasterData$MasterData$AFlogical == 1])[1] ,summary(MasterData$CVVHDialysis[MasterData$MasterData$AFlogical == 0])[1]) , c(  sum(MasterData$MasterData$AFlogical == 1), sum(MasterData$MasterData$AFlogical == 0) ))


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
  if(!is.numeric(MasterData[MasterData$MasterData$AFlogical == 0 , i])){next}
  BC_PlotCompareSingleHists(MasterData[MasterData$MasterData$AFlogical == 0 , i] , MasterData[MasterData$MasterData$AFlogical == 1 , i] , main =  names(MasterData)[i])
  
  print(names(MasterData)[i])
  print(t.test(MasterData[MasterData$MasterData$AFlogical == 1 , i] , MasterData[MasterData$MasterData$AFlogical == 0 , i]))
  print(sqrt(var(MasterData[MasterData$MasterData$AFlogical == 1 , i] , na.rm = T)))
  print(sqrt(var(MasterData[MasterData$MasterData$AFlogical == 0 , i] , na.rm = T)))
  tteststruct[counter] <- (t.test(MasterData[MasterData$MasterData$AFlogical == 0 , i] , MasterData[MasterData$MasterData$AFlogical == 1 , i])$p.value)
  counter <- counter + 1
}

DataForLogistic <- data.frame(MasterData[  ,c(1, PostOpIndexes) ])

for(i in 2:dim(DataForLogistic)[2]){
  model <- (glm(formula = as.formula(paste0("MasterData$AFlogical ~" , names(DataForLogistic)[i]) )
                ,family=binomial(link='logit') , data=DataForLogistic))
  # odds ratio
  # print( names(DataForLogistic)[i] )
  print(exp(coef(model)[2]))
  print( exp(confint(model)[2,] ))  
  print( summary( model )$coefficients[2,4] )
}

##### SOFA Comparison #####

NamesofSofaScoreVariables <- c( 'FiO26num',
                                'ArtPO2',
                                'VentilatedInOpLogical',
                                'GCSnum',
                                'Platelets',
                                'Creatinine',
                                'Bilirubin',
                                'Dopamine',
                                'Adrenaline',
                                'ReliableART.M',
                                'PaO2OverFiO2'
)

IndiciesofSofaScoreVariables <- which(names(MasterData) %in% NamesofSofaScoreVariables)

FC_StepWiseForwardAUC(IndiciesofSofaScoreVariables , MasterData)

formulaformodel <- FB_CreateFormula('MasterData$AFlogical' , IndiciesofSofaScoreVariables , MasterData)
AUCSOFA1 <- FC_CalculateCrossValidatedROC(formulaformodel ,IndiciesofSofaScoreVariables, MasterData)

model <- (glm(formula = formulaformodel,family=binomial(link='logit') , data=DataForLogistic))
summary( model )

formulaformodel <- FB_CreateFormula('MasterData$AFlogical' , which(names(MasterData) == 'SOFA') , MasterData)
AUCSOFA2 <- FC_CalculateCrossValidatedROC(formulaformodel ,which(names(MasterData) == 'SOFA'), MasterData)

model <- (glm(formula =MasterData$AFLogical ~ SOFA
              ,family=binomial(link='logit') , data=DataForLogistic))
summary( model )


##### #####

DataForLogistic <- POM_SampledImputation(data.frame(MasterData[  ,  ]))
model <- (glm(formula = FB_CreateFormula('MasterData$AFlogical' , c(PreoperativeIndices,PostOpIndexes) , MasterData)
              ,family=binomial(link='logit') , data=DataForLogistic))
summary( model )

StepResults2 <- matrix(0 , 100, length(c(PreoperativeIndices,PostOpIndexes)) )
colnames(StepResults2) <- names(MasterData[c(PreoperativeIndices,PostOpIndexes)])

for(ii in 1:100){
  DataForLogistic <- POM_SampledImputation(data.frame(MasterData[  ,  ]))
  model <- (glm(formula = FB_CreateFormula('MasterData$AFlogical' , c(PreoperativeIndices,PostOpIndexes) , MasterData)
                ,family=binomial(link='logit') , data=DataForLogistic))
  stepoutput <- stepAIC(model)
  StepResults2[ii,] <- apply(as.matrix(names(MasterData[c(PreoperativeIndices,PostOpIndexes)])) , 1 , function(X){grepl(X ,as.character(stepoutput$formula)[3])} )
}

indexesfromstep <- which(names(MasterData) %in% colnames(StepResults2[ , apply(StepResults2 , 2 , function(X){ (sum(X)/length(X))*100 })>50]))

tmp <- (cov2cor(cov(StepResults2))[apply(StepResults2 , 2 , function(X){ (sum(X)/length(X))*100 })>50 , ])

LogisticProbility <- matrix(0 , dim(DataForLogistic)[1],1)
for(i in 1:dim(LogisticProbility)[1]){
  DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
  model <- glm(formula =MasterData$AFLogical ~ 
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

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$MasterData$AFlogical == 1))
BC_CalculateAreaUnderCurve(PerformanceSweep)

x11()
BC_PlotCompareSingleHists(LogisticProbility[MasterData$MasterData$AFlogical ==0] , LogisticProbility[MasterData$MasterData$AFlogical ==1] , main ='Histogram of Logistic Model Output'  )
t.test(LogisticProbility[MasterData$MasterData$AFlogical ==0] ,  LogisticProbility[MasterData$MasterData$AFlogical ==1])

x11()
plot(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$MasterData$AFlogical ==0] , seq(0,1,0.01) ) , col = 'blue' , type ='l' , xlab = 'Output Logistic Model' , ylab = 'Culmulative Prob')
lines(seq(0,1,0.01) , FM_CalculateCDFS(LogisticProbility[MasterData$MasterData$AFlogical ==1] , seq(0,1,0.01) ) , col = 'red' , type ='l')
abline(v = 0.2)
title('CDF Analysis')

PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$MasterData$AFlogical == 1))
ROCplot2 <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot2 <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept =( 1-0.2)) + geom_hline(yintercept = 0.2)+  ggtitle('PPV vs NPV Curves')

x11(20,14)
grid.arrange(ROCplot2 , NPVPPVPlot2 + xlim(( 1-0.2) , 1) , nrow= 2)

x11(20,14)
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$MasterData$AFlogical == 1), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))

###### AUC Selection Criteria #####

FC_StepWiseForwardAUC( c(PreoperativeIndices) , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
model <- glm(formula =MasterData$AFLogical ~ 
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

{LogisticProbility <- predict(model , DataForLogistic , type = c('response'))
  
  PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,DataForLogistic$MasterData$AFlogical == 1))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')
  
  BC_CalculateAreaUnderCurve(PerformanceSweep3)
  
  phist <- ggplot(data.frame( p =LogisticProbility[MasterData$MasterData$AFlogical ==0]) , aes(p)) + geom_density(col = rgb(0,0,1,alpha = 0.5)) +
    geom_density(data = data.frame( p =LogisticProbility[MasterData$MasterData$AFlogical ==1]) , aes(p), col = rgb(1,0,0,alpha = 0.5)) + 
    ggtitle('Histogram of Model Output')
  
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$MasterData$AFlogical == 1), BinWidth = 0.1))
  EmpericalProbabilityPlot <- BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities')
  
  
  
  x11(20,20)
  grid.arrange(phist , EmpericalProbabilityPlot , ROCplot , NPVPPVPlot , ncol =2 , nrow = 2)
}
BC_CalculateAreaUnderCurve(PerformanceSweep3)

FC_StepWiseForwardAUC( c(PostOpIndexes) , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
model <- glm(formula =MasterData$AFLogical ~ 
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

{LogisticProbility <- predict(model , DataForLogistic , type = c('response'))
  
  PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,DataForLogistic$MasterData$AFlogical == 1))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')
  
  BC_CalculateAreaUnderCurve(PerformanceSweep3)
  
  phist <- ggplot(data.frame( p =LogisticProbility[MasterData$MasterData$AFlogical ==0]) , aes(p)) + geom_density(col = rgb(0,0,1,alpha = 0.5)) +
    geom_density(data = data.frame( p =LogisticProbility[MasterData$MasterData$AFlogical ==1]) , aes(p), col = rgb(1,0,0,alpha = 0.5)) + 
    ggtitle('Histogram of Model Output')
  
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$MasterData$AFlogical == 1), BinWidth = 0.1))
  EmpericalProbabilityPlot <- BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities')
  
  
  
  x11(20,20)
  grid.arrange(phist , EmpericalProbabilityPlot , ROCplot , NPVPPVPlot , ncol =2 , nrow = 2)
}
FC_StepWiseForwardAUC( c(PreoperativeIndices , PostOpIndexes ) , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData))
model <- glm(formula =MasterData$AFLogical ~ 
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
               WBC
             , family = binomial(link = "logit"),  data = DataForLogistic)
summary(model)

{LogisticProbility <- predict(model , DataForLogistic , type = c('response'))
  
  PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,DataForLogistic$MasterData$AFlogical == 1))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')
  
  BC_CalculateAreaUnderCurve(PerformanceSweep3)
  
  phist <- ggplot(data.frame( p =LogisticProbility[MasterData$MasterData$AFlogical ==0]) , aes(p)) + geom_density(col = rgb(0,0,1,alpha = 0.5)) +
    geom_density(data = data.frame( p =LogisticProbility[MasterData$MasterData$AFlogical ==1]) , aes(p), col = rgb(1,0,0,alpha = 0.5)) + 
    ggtitle('Histogram of Model Output')
  
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$MasterData$AFlogical == 1), BinWidth = 0.1))
  EmpericalProbabilityPlot <- BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities')
  
  
  
  x11(20,20)
  grid.arrange(phist , EmpericalProbabilityPlot , ROCplot , NPVPPVPlot , ncol =2 , nrow = 2)
}

#### Missing Data ####
pecentagemissingplot <- ggplot(data.frame(x = 1:length(PerentageMissing[c(PreoperativeIndices)]) , y = PerentageMissing[c(PreoperativeIndices)]) , aes(x , y)) + 
  geom_point() +
  ggtitle('Percentage of Missing Values: Pre-operation Variables') + 
  xlab('Variable Name') +
  ylab('% missing') +
  #scale_x_continuous(breaks = c(1:length(PerentageMissing[c(PreoperativeIndices)])) )  
  scale_x_continuous(breaks = c(1:length(PerentageMissing[c(PreoperativeIndices)])) , labels=c(names(PerentageMissing[c(PreoperativeIndices)]) ))  + theme(axis.text.x = element_text(angle = 45, hjust = 1))


pecentagemissingplot2 <- ggplot(data.frame(x = 1:length(PerentageMissing[c(PostOpIndexes)]) , y = PerentageMissing[c(PostOpIndexes)]) , aes(x , y)) + 
  geom_point() +
  ggtitle('Percentage of Missing Values: Post-operation Variables') + 
  xlab('Variable Name') +
  ylab('% missing') +
  #scale_x_continuous(breaks = c(1:length(PerentageMissing[c(PostOpIndexes)])) )  
  scale_x_continuous(breaks = c(1:length(PerentageMissing[c(PostOpIndexes)])) , labels=c(names(PerentageMissing[c(PostOpIndexes)]) ))  + theme(axis.text.x = element_text(angle = 45, hjust = 1))

x11()
grid.arrange(pecentagemissingplot , pecentagemissingplot2 , nrow =1 , ncol = 2)
