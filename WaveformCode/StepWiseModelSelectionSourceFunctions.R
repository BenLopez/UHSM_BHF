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
    DataForLogistic <- POM_SampledImputation(MasterPreOpData = MasterData[ , names(MasterData) %in% all.vars(formulaformodel) ])
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
  if( length(unique(PreoperativeIndices[c(which(matrixforstep[VariabesInModel,] ==1) , which(matrixforstep[VariabesOutModel[ii],] ==1) )])) == length(PreoperativeIndices[c(which(matrixforstep[VariabesInModel,] ==1) )]) ){
    next
  }
  if(j == 2){
  formulaformodel <- FB_CreateFormula('AFLogical' , unique(PreoperativeIndices[c(which(matrixforstep[VariabesInModel,] == 1) , which(matrixforstep[VariabesOutModel[ii],] ==1) )]) , MasterData)
  }else{
  formulaformodel <- FB_CreateFormula('AFLogical' , unique(PreoperativeIndices[c(which(apply(matrixforstep[VariabesInModel,] , 2 , sum) >= 1) , which(matrixforstep[VariabesOutModel[ii],] ==1) )]) , MasterData)
  }
  print(formulaformodel)
  InitialAUC[ii,j] <- FC_CalculateCrossValidatedROC(formulaformodel ,PreoperativeIndices, MasterData)
  print(InitialAUC[ii,j])
}

if(  round(max(InitialAUC[,j]) , 3) >= round(AUCTop , 3) ){
  maxindex <- which.max(InitialAUC[,j])
  VariabesInModel <- c(VariabesInModel, VariabesOutModel[maxindex])
  VariabesOutModel <- VariabesOutModel[-which(VariabesOutModel == VariabesOutModel[maxindex]) ]
  AUCTop <- max(InitialAUC[,j])
  next
  }
if(  round(max(InitialAUC[,j]) , 3) < round(AUCTop , 3) ){
  break
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
