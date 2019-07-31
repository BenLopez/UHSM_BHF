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
FB_CalulateCrossValidatedProbilities<- function( formulaformodel ,PreoperativeIndices, MasterData){
  LogisticProbility <- matrix(0,length(MasterData$AFLogical),1)
  set.seed(1)
  for(i in 1:dim(LogisticProbility)[1]){
    DataForLogistic <- POM_SampledImputation(MasterPreOpData = MasterData)
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
FC_StepWiseForwardAUC <- function(PreoperativeIndices , MasterData){
InitialAUC <- matrix(0 , length(PreoperativeIndices)[1] , length(PreoperativeIndices)[1])
for(ii in 1:dim(InitialAUC)[1]){
formulaformodel <- FB_CreateFormula('AFLogical' , PreoperativeIndices[ii] , MasterData)
InitialAUC[ii,1] <- FC_CalculateCrossValidatedROC(formulaformodel ,PreoperativeIndices, MasterData)
print(names(MasterData)[PreoperativeIndices[ii]])
print(InitialAUC[ii,1])
}

maxindex <- which.max(InitialAUC[,1])
VariabesInModel <- PreoperativeIndices[maxindex]
VariabesOutModel <- PreoperativeIndices[-maxindex]

print(names(MasterData)[VariabesInModel])
AUCTop <- max(InitialAUC[,1])
for(j in 2:length(PreoperativeIndices)){
for(ii in 1:length(VariabesOutModel)){
  formulaformodel <- FB_CreateFormula('AFLogical' , c(VariabesInModel, VariabesOutModel[ii]) , MasterData)
  InitialAUC[ii,j] <- FC_CalculateCrossValidatedROC(formulaformodel ,PreoperativeIndices, MasterData)
  print(names(MasterData)[VariabesOutModel[ii]])
  print(InitialAUC[ii,j])
}

if(  round(max(InitialAUC[,j]) , 3) >= round(AUCTop , 3) ){
  maxindex <- which.max(InitialAUC[,j])
  VariabesInModel <- c(VariabesInModel, VariabesOutModel[maxindex])
  VariabesOutModel <- VariabesOutModel[-maxindex]
  AUCTop <- max(InitialAUC[,j])
  print(names(MasterData)[VariabesInModel])
  next
  }
if(  round(max(InitialAUC[,j]) , 3) < round(AUCTop , 3) ){
  print(names(MasterData)[VariabesInModel])
  break
}
}
return(list(AUCTop , FB_CreateFormula('AFLogical' , c( VariabesInModel ) , MasterData) , InitialAUC))
}

