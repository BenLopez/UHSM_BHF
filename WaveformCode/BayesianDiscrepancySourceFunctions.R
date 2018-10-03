BD_CalulateSenSpecNPVPPV<-function(ProbCalibStruct , prob_thresh){
  DecsionVector <- ProbCalibStruct[ , 1] > prob_thresh
  
  N <- sum(ProbCalibStruct[ , 2] == 0)
  P <- sum(ProbCalibStruct[ , 2] == 1)
  TP <- sum((DecsionVector == 1)*(ProbCalibStruct[ , 2] == 1))
  TN <- sum((DecsionVector == 0)*(ProbCalibStruct[ , 2] == 0))
  FP <- sum((DecsionVector == 1)*(ProbCalibStruct[ , 2] == 0))
  FN <- sum((DecsionVector == 0)*(ProbCalibStruct[ , 2] == 1))
  
  Sensitivity <- TP/P 
  Specifictity <- TN/N
  PPV <- TP /(TP + FP)
  NPV <- TN / (TN + FN)   
  
  output <- setNames(list(Sensitivity , Specifictity , PPV , NPV)  , c('Sen' , 'Spec' , 'PPV' , 'NPV'))
  return(output)
}
