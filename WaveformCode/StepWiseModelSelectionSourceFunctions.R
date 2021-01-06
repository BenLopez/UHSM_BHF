FB_CreateFormula <- function(y , x , DataForLogistic , intercept =1){
 
  output <- paste0(y, "~" , names(DataForLogistic)[x[1]])
    if(length(x) > 1){
      for(ii in 2:length(x)){
      output <- paste0(output , '+' ,  names(DataForLogistic)[x[ii]])
      }
    }
  if(intercept == 1){
  return(as.formula(output))}
  if(intercept == 0){
    return(as.formula(paste0(output , '-1')))}
}
FB_CalulateCrossValidatedProbilities<- function( formulaformodel , PreoperativeIndices , MasterData){
  LogisticProbility <- matrix(0,length(MasterData$AFLogical),1)
  set.seed(1)
  for(i in 1:dim(LogisticProbility)[1]){
    DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[ , names(MasterData) %in% all.vars(formulaformodel) ]) )
    model <- glm(formula = formulaformodel , family = binomial(link = "logit"), data = DataForLogistic[-i,])
    LogisticProbility[i,] <- predict(object= model , newdata = DataForLogistic[i , ] , type = c( "response") )
  }
  return(LogisticProbility)
}
FC_CalculateCrossValidatedROC <- function(formulaformodel ,PreoperativeIndices, MasterData){
  LogisticProbility <- FB_CalulateCrossValidatedProbilities(formulaformodel,PreoperativeIndices, MasterData)
  PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(LogisticProbility ,MasterData$AFLogical == 1))
  roc_obj <- BC_CalculateAreaUnderCurve(PerformanceSweep)
  return(roc_obj)
}  
FC_StepWiseForwardAUC <- function(PreoperativeIndices , MasterData , matrixforstep = diag(length(PreoperativeIndices))){
  
InitialAUC <- matrix(0 , dim(matrixforstep)[1], dim(matrixforstep)[1])

for(ii in 1:dim(InitialAUC)[1]){
formulaformodel <- FB_CreateFormula('AFLogical' , PreoperativeIndices[which(matrixforstep[ii,] ==1 )] , MasterData)
print(formulaformodel)
InitialAUC[ii,1] <- FC_CalculateCrossValidatedROC(formulaformodel = formulaformodel ,PreoperativeIndices = PreoperativeIndices,
                                                  MasterData = MasterData )
print(InitialAUC[ii,1])
}

maxindex <- which.max(InitialAUC[,1])
VariabesInModel <- maxindex
VariabesOutModel <- c(1:dim(matrixforstep)[1])[-maxindex]

AUCTop <- max(InitialAUC[,1])
for(j in 2:dim(matrixforstep)[1]){
  
for(ii in 1:length(VariabesOutModel)){
  if(j == 2){
  formulaformodel <- FB_CreateFormula('AFLogical' , unique(PreoperativeIndices[c(which(matrixforstep[VariabesInModel,] == 1) , which(matrixforstep[VariabesOutModel[ii],] ==1) )]) , MasterData)
  }else{
  formulaformodel <- FB_CreateFormula('AFLogical' , unique(PreoperativeIndices[c(which(apply(matrixforstep[VariabesInModel,] , 2 , sum) >= 1) , which(matrixforstep[VariabesOutModel[ii],] ==1) )]) , MasterData)
  }
  print(formulaformodel)
  InitialAUC[ii,j] <- FC_CalculateCrossValidatedROC(formulaformodel ,PreoperativeIndices, MasterData)
  print(paste0(round(AUCTop,3),' ', round(InitialAUC[ii,j],3) ) )
}

if(  round(max(InitialAUC[,j]) , 3) <= round(AUCTop , 3) ){
  break
}
  
if(  round(max(InitialAUC[,j]) , 3) > round(AUCTop , 3) ){
    maxindex <- which.max(InitialAUC[,j])
    VariabesInModel <- c(VariabesInModel, VariabesOutModel[maxindex])
    VariabesOutModel <- VariabesOutModel[-which(VariabesOutModel == VariabesOutModel[maxindex]) ]
    AUCTop <- max(InitialAUC[,j])
    next
}
}
if(length(VariabesInModel)>1){
return(list(AUCTop , FB_CreateFormula('AFLogical' , PreoperativeIndices[which(apply(matrixforstep[VariabesInModel,] , 2 , sum ) >= 1)] , MasterData) , InitialAUC))
}else{
return(list(AUCTop , FB_CreateFormula('AFLogical' , PreoperativeIndices[which(matrixforstep[VariabesInModel,] ==1)] , MasterData) , InitialAUC))
}
  }
FC_ExtractIndiciesFromFormula <- function(NamesofVariables , ModelFormula){
  output <- matrix(0 , 1 , length(NamesofVariables))
  output[, which(NamesofVariables %in% all.vars(ModelFormula)[-1]) ] <- 1 
  return(output)
}
FC_BackwardStep <- function(PreoperativeIndices , MasterData){
  
  InitialAUC <- matrix(0 , length(PreoperativeIndices), 1)
  
  formulaformodel <- FB_CreateFormula('AFLogical' , PreoperativeIndices , MasterData)
  print(formulaformodel)
  AUCTop <- FC_CalculateCrossValidatedROC(formulaformodel = formulaformodel ,PreoperativeIndices = PreoperativeIndices,
                                          MasterData = MasterData )
  for(ii in 1:length(PreoperativeIndices)){
    formulaformodel <- FB_CreateFormula('AFLogical' , PreoperativeIndices[-ii] , MasterData)
    print(formulaformodel)
    print(names(MasterData)[PreoperativeIndices[ii]] )
    InitialAUC[ii,1] <- FC_CalculateCrossValidatedROC(formulaformodel = formulaformodel ,PreoperativeIndices = unique(PreoperativeIndices),MasterData = MasterData )
    print(paste0(round(AUCTop,3),' ', round(InitialAUC[ii,1],3) ) )
  }
  
  if(sum((InitialAUC - AUCTop) >0) ==0 ){
    return(PreoperativeIndices)
  }
  if(sum((InitialAUC - AUCTop) >0) > 0 ){
    return(PreoperativeIndices[-which.max(InitialAUC[,1])])
  }
}
FC_StepWiseBackwardAUC <- function(PreoperativeIndices , MasterData){
Inmodel <- c(1)
StartLength <- length(PreoperativeIndices)
EndLength <- 0
Inmodel <- PreoperativeIndices
while( StartLength != EndLength ){
  StartLength <- length(Inmodel)
  Inmodel <-  FC_BackwardStep( Inmodel , MasterData )  
  EndLength <- length(Inmodel)  
}
  return(Inmodel)
}
FC_ForwardStep <- function(VariablesInModel ,VariablesOutModel, PreoperativeIndices , MasterData , matrixforstep = diag(length(PreoperativeIndices))){
  
  InitialAUC <- matrix(0 , length(VariablesOutModel), 1)
  formulaformodel <- FB_CreateFormula('AFLogical' , VariablesInModel , MasterData)
  AUCTop <- FC_CalculateCrossValidatedROC(formulaformodel = formulaformodel ,PreoperativeIndices = VariablesInModel ,
                                          MasterData = MasterData )
  
  for(ii in 1:length(VariablesOutModel)){
    formulaformodel <- FB_CreateFormula('AFLogical' , unique( c(VariablesInModel ,PreoperativeIndices[  which(matrixforstep[VariablesOutModel[ii],] == 1) ]) ) , MasterData)
    print(formulaformodel)
    InitialAUC[ii,1] <- FC_CalculateCrossValidatedROC(formulaformodel = formulaformodel ,PreoperativeIndices = PreoperativeIndices ,
                                                      MasterData = MasterData )
    print(paste0(AUCTop , ' ', InitialAUC[ii,1]))
  }
  
  if(round(AUCTop,3) > round(max(InitialAUC[,1]),3)) {
    return(list(VariablesInModel,VariablesOutModel) )
  }else{
    VariablesInModel <-  unique( c(VariablesInModel , PreoperativeIndices[ which(matrixforstep[VariablesOutModel[which.max(InitialAUC[,1])],] ==1)] ))
    VariablesOutModel <- VariablesOutModel[-which.max(InitialAUC[,1])] 
    return(list(VariablesInModel,VariablesOutModel))
  }
  
}
FC_StepWiseForwardBackwardAUC <- function(PreoperativeIndices , MasterData , matrixforstep = diag(length(PreoperativeIndices))){
  
  
  InitialAUC <- matrix(0 , dim(matrixforstep)[1], 1)
  
  for(ii in 1:dim(InitialAUC)[1]){
    if(sum(matrixforstep[ii,]) ==0){next}
    formulaformodel <- FB_CreateFormula('AFLogical' , PreoperativeIndices[which(matrixforstep[ii,] ==1 )] , MasterData)
    print(formulaformodel)
    InitialAUC[ii,1] <- FC_CalculateCrossValidatedROC(formulaformodel = formulaformodel ,PreoperativeIndices = PreoperativeIndices,
                                                      MasterData = MasterData )
    print(InitialAUC[ii,1])
  }
  
  maxindex <- which.max(InitialAUC[,1])
  VariablesInModel <- PreoperativeIndices[which(matrixforstep[maxindex,] ==1 )]
  VariablesOutModel <- c(1:dim(matrixforstep)[1])[-maxindex]
  numberinbefore <- length(VariablesInModel)
  numberinafter <- 0
  
while(numberinbefore != numberinafter){
    # backward selection
if( exists('forwardoutput') ){
if( length(forwardoutput[[1]])>3){
  disp('Backward step')
      VariablesInModel <- FC_StepWiseBackwardAUC( PreoperativeIndices = VariablesInModel , MasterData)
}
}
numberinbefore <- length(VariablesInModel)
# Forward selection
disp('Forward step')    
forwardoutput <- FC_ForwardStep(VariablesInModel,VariablesOutModel,PreoperativeIndices,MasterData , matrixforstep)
    disp(paste0('Number VariablesOutModel ',length(forwardoutput[[2]])))
    if(length(forwardoutput[[2]]) == 0 ){break}
    VariablesInModel <- forwardoutput[[1]]
    VariablesOutModel <- forwardoutput[[2]]
    numberinafter <- length(VariablesInModel)
}
  return(list(VariablesInModel, FB_CreateFormula('AFLogical' , VariablesInModel , MasterData) ))
}
FC_CalculateLeaveOneOutStats <- function( Formula , MasterData){
  
  indiciesinmodel <- which(FC_ExtractIndiciesFromFormula( names(MasterData) , Formula ) ==1)
  formulaformodel <- FB_CreateFormula('AFLogical' , indiciesinmodel[-i] , MasterData)
  AUCTop <- FC_CalculateCrossValidatedROC(formulaformodel ,PreoperativeIndices = indiciesinmodel , MasterData)
  
  if(length(indiciesinmodel)>1){
  output <- matrix(0 , length(indiciesinmodel) , 1)
  for(i in 1:dim(output)[1]){
    formulaformodel <- FB_CreateFormula('AFLogical' , indiciesinmodel[-i] , MasterData)
    output[i,] <- FC_CalculateCrossValidatedROC(formulaformodel ,PreoperativeIndices = indiciesinmodel , MasterData) - AUCTop
  }
  rownames(output ) <- names(MasterData)[indiciesinmodel]
  return(output)}
  if(length(indiciesinmodel)<=1){
  output <- matrix(NA , 1, 1)
  return(NA)
  }
}
FC_AlterModelSummary <- function(model , MasterData){
  
  output <- xtable(model)
  output[1:dim(output)[1],3] <- 0
  colnames(output)[3] <- 'Left Out CV AUC'
  
  leaveoneoutauc <- FC_CalculateLeaveOneOutStats(Formula = model$formula , MasterData )
  
  if(!is.null(dim(leaveoneoutauc)) ){
  tmpindiciesvector <- unlist(apply(apply(as.matrix(rownames(leaveoneoutauc)) , 1, function(X){grepl(X , rownames(output))} ),1,function(X){which(X)[1]}))
  output[1:dim(output)[1] , 3] <- leaveoneoutauc[tmpindiciesvector]
  output[1 , 3] <- 'NA'}else{
  output[1:dim(output)[1] , 3] <- 'NA'
  }
  return(output)
}
FB_CorrectModelTableNames <- function(modeltable){
  
  {
  tmp <- rownames(modeltable)
  tmp <- gsub(pattern = 'Preop' , replacement = 'Pre-operative ' , x = tmp) 
  tmp <- gsub(pattern = 'PreOp' , replacement = 'Pre-operative ' , x = tmp) 
  tmp <- gsub(pattern = 'PostOp' , replacement = 'Post-operative ' , x = tmp) 
  tmp <- gsub(pattern = 'Postop' , replacement = 'Post-operative ' , x = tmp) 
  tmp <- gsub(pattern = 'EjectionFractionCategoryLVEF' , replacement = 'Ejection Fraction Category(LVEF)' , x = tmp) 
  tmp <- gsub(pattern = ' Bili' , replacement = ' Bilirubin' , x = tmp)
  tmp <- gsub(pattern = 'Alb' , replacement = 'Albumin' , x = tmp)
  tmp <- gsub(pattern = 'PLT' , replacement = 'Platelets' , x = tmp) 
  tmp <- gsub(pattern = '_sq' , replacement = ' Squared' , x = tmp)
  tmp <- gsub(pattern = 'ExtracardiacArteriopathyExtracardiacArteriopathy' , replacement = 'Extracardiac Arteriopathy' , x = tmp)
  tmp <- gsub(pattern = 'Reliable' , replacement = '' , x = tmp)
  tmp <- gsub(pattern = 'ART.M' , replacement = 'ART Mean' , x = tmp)
  tmp <- gsub(pattern = 'ART.S' , replacement = 'ART Systolic' , x = tmp)
  tmp <- gsub(pattern = 'AKIUODayOneTRUE' , replacement = 'AKI UO Day One(Yes)' , x = tmp)
  tmp <- gsub(pattern = 'DopamineStandardised' , replacement = 'Dopamine Standardised' , x = tmp)
  tmp <- gsub(pattern = 'AdrenalineStandardised' , replacement = 'Adrenaline Standardised' , x = tmp)
  tmp <- gsub(pattern = 'MilrinoneStandardised' , replacement = 'Milrinone Standardised' , x = tmp)
  tmp <- gsub(pattern = 'DobutamineStandardised' , replacement = 'Dobutamine Standardised' , x = tmp)
  tmp <- gsub(pattern = 'NoradrenalineStandardised' , replacement = 'Noradrenaline Standardised' , x = tmp)
  tmp <- gsub(pattern = 'CardioVascularModel' , replacement = 'Cardiovascular Model' , x = tmp)
  tmp <- gsub(pattern = 'RespiratoryModel' , replacement = 'Respiratory Model' , x = tmp)
  tmp <- gsub(pattern = 'OperativeModel' , replacement = 'Operative Model' , x = tmp)
  tmp <- gsub(pattern = 'ProcDetails' , replacement = '' , x = tmp)
  tmp <- gsub(pattern = 'Planned.Valve.Surgery' , replacement = 'Planned Valve Surgery' , x = tmp)
  tmp <- gsub(pattern = 'Mitral Valve' , replacement = '(Mitral Valve)' , x = tmp)
  tmp <- gsub(pattern = 'Multiple' , replacement = '(Multiple)' , x = tmp)
  tmp <- gsub(pattern = 'None' , replacement = '(None)' , x = tmp)
  tmp <- gsub(pattern = 'Tricuspid Valve' , replacement = '(Tricuspid Valve)' , x = tmp)
  tmp <- gsub(pattern = 'Male' , replacement = '(Male)' , x = tmp)
  tmp <- gsub(pattern = 'Thoracic.AortaThoracic.Aorta' , replacement = 'Thoracic Aorta' , x = tmp)
  tmp <- gsub(pattern = ' Creat' , replacement = ' Creatinine' , x = tmp)
  tmp <- gsub(pattern = 'HistoryOfNeurologicalDysfunctionHistoryOfNeurologicalDysfunction' , replacement = 'History Of Neurological Dysfunction' , x = tmp)
  tmp <- gsub(pattern = 'Bilirubinrubin' , replacement = 'Bilirubin' , x = tmp)
  tmp <- gsub(pattern = 'FiO26num' , replacement = 'FiO_2' , x = tmp)
  tmp <- gsub(pattern = 'GCSnum' , replacement = 'GCS Number' , x = tmp)
  tmp <- gsub(pattern = 'ArtPO2' , replacement = 'PaO_2' , x = tmp)
  tmp <- gsub(pattern = 'NYHAGradeNYHA' , replacement = 'NYHA Grade' , x = tmp)
  tmp <- gsub(pattern = 'DiabetesDiabetes' , replacement = 'Diabetes' , x = tmp)
  tmp <- gsub(pattern = 'IntravenousInotropesPriorToAnaesthesiaIntravenousInotropesPriorToAnaesthesia' ,
              replacement = 'Intravenous Inotropes Prior To Anaesthesia' , x = tmp)
  tmp <- gsub(pattern = 'UrgencyUrgent' ,replacement = 'Urgency(Urgent)' , x = tmp)
  tmp <- gsub(pattern = 'UrgencyEmergency' ,replacement = 'Urgency(Emergency)' , x = tmp)
  tmp <- gsub(pattern = 'UrgencySalvage' ,replacement = 'Urgency(Salvage)' , x = tmp)
  tmp <- gsub(pattern = 'ReducedFALSE' ,replacement = 'Procedure Details CABG(Yes)' , x = tmp)
  tmp <- gsub(pattern = 'dBili' ,replacement = 'dBilirubin' , x = tmp)
  tmp <- gsub(pattern = 'dAlb' ,replacement = 'dAlbumin' , x = tmp)
  tmp <- gsub(pattern = 'dPLT' ,replacement = 'dPlatelets' , x = tmp)
  tmp <- gsub(pattern = 'dCreat' ,replacement = 'dCreatinine' , x = tmp)
  tmp <- gsub(pattern = 'SpO2' ,replacement = 'SpO_2' , x = tmp)
  tmp <- gsub(pattern = 'PaO2OverFiO2' ,replacement = 'PaO_2 Over FiO_2' , x = tmp)
  tmp <- gsub(pattern =  "HistoryOfPulmonaryDiseaseHistoryOfPulmonaryDisease"   , replacement = "History Of Pulmonary Disease" , x = tmp) 
  tmp <- gsub(pattern =  "PreviousCardiacSurgeryPreviousCardiacSurgery"   , replacement = "Previous Cardiac Surgery" , x = tmp) 
  tmp <- gsub(pattern =  "Recent.MIRecent.MI"   , replacement = "Recent MI" , x = tmp)
  tmp <- gsub(pattern =  "Albuminumin"   , replacement = "Albumin" , x = tmp)
  tmp <- gsub(pattern =  "AorticProcedure"   , replacement = "Procedure" , x = tmp)
  tmp <- gsub(pattern =  "VentilatedPre-operative erationVentilatedPre-operative eration"   , replacement = "Ventilated Pre-operative" , x = tmp)
  tmp <- gsub(pattern =  "Active.EndocarditisActive.Endocarditis"   , replacement = "Active Endocarditis" , x = tmp)
  tmp <- gsub(pattern =  "CardiogenicShock\\_pre\\_OperationCardiogenicShock\\_pre\\_Operation"   , replacement = "Cardiogenic Shock Pre-operation" , x = tmp)
  tmp <- gsub(pattern =  "Pre-operative SupportPre-operative Support"   , replacement = "Pre-operative Support" , x = tmp)
  tmp <- gsub(pattern =  "ValveProcedure Details(Valve)"   , replacement = "Procedure Details Valve(Yes)" , x = tmp)
  }
  
  rownames(modeltable) <- tmp 
  
  modeltable[,1] <- signif( as.numeric(modeltable[,1]) , 4)
  modeltable[,2] <- signif( as.numeric(modeltable[,2]) , 4)
  modeltable[1,3] <- 0
  modeltable[,3] <- signif( as.numeric(modeltable[,3])  , 4)
  modeltable[,4] <- signif( as.numeric(modeltable[,4]) , 4)
  modeltable[,4][ modeltable[,4] < 0.001] <- '<0.001'
  
  return(modeltable)
}
FB_CorrectUnivariateTableNames <- function(UnivariateTable){
  
  if(sum(UnivariateTable[,1] == " TRUE") != 0 ){
    UnivariateTable <- UnivariateTable[-which(UnivariateTable[,1] == " TRUE"),]
  if(sum(UnivariateTable[,1] == "FALSE") !=0 ){
    UnivariateTable <- UnivariateTable[-which(UnivariateTable[,1] == "FALSE"),]}
  }
  if(sum(UnivariateTable[,1] == "Not Urgent") !=0 ){
    UnivariateTable <- UnivariateTable[-which(UnivariateTable[,1] == "Not Urgent")+1,]
  UnivariateTable <- UnivariateTable[-which(UnivariateTable[,1] == "Not Urgent"),]
  }

  {tmp <- UnivariateTable[,1]
    
    
    # Categorical variables
  tmp <- gsub(pattern = 'Female' , replacement = 'Gender(Female)' , x = tmp) 
  tmp <- gsub(pattern = 'Male' , replacement = 'Gender(Male)' , x = tmp) 
#  tmp <- gsub(pattern = 'CABG and Valve' , replacement = 'Procedure Details(CABG and Valve)' , x = tmp) 
  tmp <- gsub(pattern = 'LVEF' , replacement = 'Ejection Fraction Category(LVEF)' , x = tmp) 
  tmp <- gsub(pattern = 'Mitral Valve' , replacement = 'Valve Surgery(Mitral Valve)' , x = tmp) 
  tmp <- gsub(pattern = 'Aortic Valve' , replacement = 'Valve Surgery(Aortic Valve)' , x = tmp) 
  tmp <- gsub(pattern = 'Tricuspid Valve' , replacement = 'Valve Surgery(Tricuspid Valve)' , x = tmp) 
  tmp <- gsub(pattern = 'None' , replacement = 'Valve Surgery(None)' , x = tmp) 
  tmp <- gsub(pattern = 'Multiple' , replacement = 'Valve Surgery(Multiple)' , x = tmp) 
  tmp <- gsub(pattern = 'Active.Endocarditis' , replacement = 'Active Endocarditis' , x = tmp) 
  tmp <- gsub(pattern =  "Atrial fibrillation/flutter"   , replacement = 'Pre-operative Heart Rhythm(AF)' , x = tmp) 
  tmp <- gsub(pattern =  "Complete heart block/pacing"   , replacement = 'Pre-operative Heart Rhythm(Complete heart block/pacing)' , x = tmp) 
  tmp <- gsub(pattern =  "Other abnormal rhythm"   , replacement = 'Pre-operative Heart Rhythm(Other)' , x = tmp) 
  tmp <- gsub(pattern =  "Sinus Rhythm"   , replacement = 'Pre-operative Heart Rhythm(Sinus Rhythm)' , x = tmp) 
  tmp <- gsub(pattern =  "Ventricular fibrillation or ventricular tachycardia"   , replacement = 'Pre-operative Heart Rhythm(Ventricular Fib or ventricular Tachy)' , x = tmp) 
  tmp <- gsub(pattern =  "Elective"   , replacement = 'Urgency(Elective)' , x = tmp) 
  tmp <- gsub(pattern =  "Emergency"   , replacement = 'Urgency(Emergency)' , x = tmp) 
  tmp <- gsub(pattern =  "Salvage"   , replacement = 'Urgency(Salvage)' , x = tmp) 
  tmp <- gsub(pattern =  "Urgent"   , replacement = 'Urgency(Urgent)' , x = tmp) 
  tmp <- gsub(pattern =  "CardiogenicShock_pre_Operation"   , replacement = 'Cardiogenic Shock (pre-operation)' , x = tmp) 
  tmp <- gsub(pattern =  "HistoryOfNeurologicalDysfunction"   , replacement = 'History Of Neurological Dysfunction' , x = tmp) 
  tmp <- gsub(pattern =  "HistoryOfPulmonaryDisease"   , replacement = 'History Of Pulmonary Disease' , x = tmp) 
  tmp <- gsub(pattern =  "Thoracic.Aorta"   , replacement = 'Thoracic Aorta' , x = tmp) 
  tmp <- gsub(pattern =  "PreviousCardiacSurgery"   , replacement = 'Previous Cardiac Surgery' , x = tmp) 
  tmp <- gsub(pattern =  "Recent.MI"   , replacement = 'Recent MI' , x = tmp) 
  tmp <- gsub(pattern =  "PreOpSupport"   , replacement = 'Pre-operative Support' , x = tmp) 
  tmp <- gsub(pattern =  "IntravenousInotropesPriorToAnaesthesia"   , replacement = 'Intravenous Inotropes Prior To Anaesthesia' , x = tmp) 
  tmp <- gsub(pattern =  "IntravenousNitratesOrAnyHeparin"   , replacement = 'Intravenous Nitrates Or Any Heparin' , x = tmp) 
  tmp <- gsub(pattern =  "ExtracardiacArteriopathy"   , replacement = 'Extracardiac Arteriopathy' , x = tmp) 
  tmp <- gsub(pattern =  "FALSE"   , replacement = 'AKI by UO(No)' , x = tmp) 
  tmp <- gsub(pattern =  "TRUE"   , replacement = 'AKI by UO(Yes)' , x = tmp) 
  tmp <- gsub(pattern =  "AnginaGrade"   , replacement = 'Angina Grade' , x = tmp) 
  
  
  # continuous variables
  tmp <- gsub(pattern =  "PreOp"   , replacement = 'Pre-operative ' , x = tmp) 
  tmp <- gsub(pattern =  "Preop"   , replacement = 'Pre-operative ' , x = tmp) 
  tmp <- gsub(pattern =  " Alb"   , replacement = ' Albumin' , x = tmp) 
  tmp <- gsub(pattern =  " Creat"   , replacement = ' Creatinine' , x = tmp) 
  tmp <- gsub(pattern =  " PLT"   , replacement = ' Platelets' , x = tmp) 
  tmp <- gsub(pattern =  " Bili"   , replacement = ' Bilirubin' , x = tmp) 
  tmp <- gsub(pattern =  "LogisticEUROScore"   , replacement = 'Logistic EUROScore' , x = tmp) 
  tmp <- gsub(pattern =  "ART.S"   , replacement = 'ART(Systolic)' , x = tmp) 
  tmp <- gsub(pattern =  "ART.M"   , replacement = 'ART(Mean)' , x = tmp) 
  tmp <- gsub(pattern =  "ART.D"   , replacement = 'ART(diastolic)' , x = tmp) 
  tmp <- gsub(pattern =  "ArtPO2"   , replacement = 'PaO_2' , x = tmp) 
  tmp <- gsub(pattern =  "NoradrenalineStandardised"   , replacement = 'Noradrenaline Standardised' , x = tmp) 
  tmp <- gsub(pattern =  "DopamineStandardised"   , replacement = 'Dopamine Standardised' , x = tmp) 
  tmp <- gsub(pattern =  "AdrenalineStandardised"   , replacement = 'Adrenaline Standardised' , x = tmp) 
  tmp <- gsub(pattern =  "MilrinoneStandardised"   , replacement = 'Milrinone Standardised' , x = tmp) 
  tmp <- gsub(pattern =  "DobutamineStandardised"   , replacement = 'Dobutamine Standardised' , x = tmp) 
  tmp <- gsub(pattern =  "Lac"   , replacement = "Lactic Acid" , x = tmp) 
  tmp <- gsub(pattern =  "GCSnum"   , replacement = "GCS Number" , x = tmp) 
  tmp <- gsub(pattern =  "SOFA"   , replacement = "SOFA Score" , x = tmp) 
  tmp <- gsub(pattern =  "logCASUS"   , replacement = "log CASUS Score" , x = tmp) 
  tmp <- gsub(pattern =  "FiO26num"   , replacement = "FiO_2" , x = tmp) 
  tmp <- gsub(pattern =  "PaO2OverFiO2"   , replacement = "PaO_2 over FiO_2" , x = tmp) 
  tmp <- gsub(pattern =  "SpO2"   , replacement = "SpO_2" , x = tmp) 
  }
    
  UnivariateTable[,12] <- signif(as.numeric(UnivariateTable[,12]) , 4)
  UnivariateTable[,1] <- tmp

  return(UnivariateTable)
}
FB_CorrectUnivariateTableValues <- function(UnivariateTableLatex){
  
  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,8]) > 100)  ,8] <- '>100'
  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,9]) > 100) ,9] <- '>100'
  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,10]) > 100) ,10] <- '>100'

  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,8]) <0.001)  ,8] <- '<0.001'
  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,9]) <0.001) ,9] <- '<0.001'
  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,10]) <0.001) ,10] <- '<0.001'
  
  UnivariateTableLatex[ , 12] <- signif(as.numeric(UnivariateTableLatex[,12]) , 4)
  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,11])<0.001),11]<- '<0.001'
  UnivariateTableLatex[which(as.numeric(UnivariateTableLatex[,7])<0.001),7]<- '<0.001'
  return(UnivariateTableLatex)
}
