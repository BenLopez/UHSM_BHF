BLBF_CalculateDensities <- function(AFModel, NAFModel , ImMatrix , AFLogical){
  ImMatrix <- as.matrix(ImMatrix)
  f_1 <- matrix(0 , dim(ImMatrix)[1] , 2)
  f_1[,1] <- predict(AFModel , x = ImMatrix)
  f_1[,2] <- predict(NAFModel , x = ImMatrix)
  return(f_1)
}
BLBF_CalculatePosteriorProbabilities <- function(AFModel, NAFModel , ImMatrix , AFLogical , alpha= c(sum(AFLogical)/length(AFLogical))){

  f_1 <- BLBF_CalculateDensities(AFModel, NAFModel , ImMatrix , AFLogical)
  PosteriorProbability <- alpha*f_1[,1]/(alpha*f_1[,1] + (1-alpha)*f_1[,2])
  return(PosteriorProbability)
} 
BLBF_ForecastingCalculateImplausabilities <- function(Data , MeanVector , CovarianceMatrix){
  return(log(apply(Data, 1 , function(X){(X - MeanVector)%*%solve(CovarianceMatrix)%*%(X - MeanVector)})))
}
BLBF_ReshapetoMatrix <- function(X){
  return(matrix(X , dim(X)[1] , dim(X)[2]*dim(X)[3]))
}
BLBF_CalculateSecondOrderStruct <- function(X){
  if(length(dim(X))>2){
  X <- BLBF_ReshapetoMatrix(X)}
  output <- BLU_CreateSOS(mu = apply(X , 2 , mean) , Sigma =  cov(X))
  return(output)
}
BLBF_CalculateForecastAdjustedBeliefs <- function(SOS , z ){
  
  lenxplusz <- length(SOS$Expectation)
  lenz <- length(z)
  lenx <- lenxplusz - length(z)

  E_X <- SOS$Expectation[(lenz + 1):lenxplusz ]
  V_X <- SOS$Covariance[(lenz + 1):lenxplusz ,  (lenz + 1):lenxplusz]
  E_z <- SOS$Expectation[1:lenz ]
  V_z <-  SOS$Covariance[1:lenz , 1:lenz]
  c_Xz <- SOS$Covariance[(lenz + 1):lenxplusz , 1:lenz]
  weights <- c_Xz%*%solve(DP_AddNugget(V_z , 1e-6*diag(diag(V_z))))     
  E_z_X <- E_X + weights%*%(z - E_z)
  V_z_X <- V_X -  weights%*%t(c_Xz)  
  
  return(setNames(list(E_z_X , V_z_X) , c('AdjustedExpectation' , 'AdjustedVariance')))
}
BLBF_CalculateAdjustedDiscrepancy<-function(AdjustedBeliefs , x){
  return(log(mahalanobis(x ,AdjustedBeliefs$AdjustedExpectation , AdjustedBeliefs$AdjustedVariance )) )
}
BLBF_CrossValidateAdjustedVersion <- function(TrainingData , TestData , LengthVector , NumberofPoints){
  
  SOS <- BLBF_CalculateSecondOrderStruct( TrainingData ) 
  AdjustedBeliefs <- BLBF_CalculateForecastAdjustedBeliefs(SOS , TestData[1:(LengthVector*(NumberofPoints))])
  AdjustedVersion <- AdjustedBeliefs$AdjustedExpectation -  TestData[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints + 1) )]
  return(AdjustedVersion)
}
BLBF_CrossValidateCalculateIm <- function(TrainingData , TestData , LengthVector , NumberofPoints){
  
  SOS <- BLBF_CalculateSecondOrderStruct( TrainingData ) 
  Im <- BLBF_CalculateAdjustedDiscrepancy(BLBF_CalculateForecastAdjustedBeliefs(SOS , TestData[1:(LengthVector*(NumberofPoints))]) , TestData[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints + 1) )])
  return(Im)
}
BLBF_CrossValidateCalculateImAll <- function(TrainingData ,  LengthVector , NumberofPoints){
  
  Implausability <- matrix(0 , dim(TrainingData)[1]  , 1)
  for(ii in 1:dim(TrainingData)[1]){
    Implausability[ii,] <- BLBF_CrossValidateCalculateIm(TrainingData[-ii , , ] , TrainingData[ii , , ] , LengthVector , NumberofPoints)
  }
  return(Implausability)  
}
BLBF_CrossValidateAdjustedVersionAll <- function(TrainingData  , LengthVector , NumberofPoints){
  
  AdjustedVersion <- matrix(0 , dim(TrainingData)[1]  , LengthVector)
  for(ii in 1:dim(TrainingData)[1]){
    AdjustedVersion[ii,] <- BLBF_CrossValidateAdjustedVersion(TrainingData[-ii , , ] , TrainingData[ii , , ] , LengthVector , NumberofPoints)
  }
  return(AdjustedVersion)  
}
BLBF_CrossValidateCalculateAdjustedVersionSOS<- function(TrainingData  , LengthVector , NumberofPoints){
  AdjustedVersions <- BLBF_CrossValidateAdjustedVersionAll(TrainingData  , LengthVector , NumberofPoints)
  SOS <- BLU_CreateSOS(mu = apply(AdjustedVersions , 2 , mean) , Sigma = cov(AdjustedVersions) )
return(SOS)
}
BLBF_CalculateForecastAdjustedVersion <- function(SOS , z , x){
  
  lenxplusz <- length(SOS$Expectation)
  lenz <- length(z)
  lenx <- lenxplusz - length(z)
  
  E_X <- SOS$Expectation[(lenz + 1):lenxplusz ]
  V_X <- SOS$Covariance[(lenz + 1):lenxplusz ,  (lenz + 1):lenxplusz]
  E_z <- SOS$Expectation[1:lenz ]
  V_z <-  SOS$Covariance[1:lenz , 1:lenz]
  c_Xz <- SOS$Covariance[(lenz + 1):lenxplusz , 1:lenz]
  weights <- c_Xz%*%solve(DP_AddNugget(V_z , 1e-6*diag(diag(V_z))))     
  E_z_X <- E_X + weights%*%(z - E_z)
  
  return(E_z_X - x)
}
BLBF_Discrepancy <- function(SOS , x){
  return(sqrt((mahalanobis(x , SOS$Expectation , SOS$Covariance))))
}
#BLBF_CalculatePosteriorProbabilityNBLA <- function( DataStructure , AFLogical , VariablesToView = -c( 13) ,  indextoview = 371  , PriorProbability =  c(sum(AFLogical)/length(AFLogical)) ){
  
 # Im <- BLBL_CrossValidateImplausabilitiesALL(TrainingData = DataStructure[, VariablesToView, indextoview] , AFLogical)
#  ImAF <- Im[,1]
#  ImNAF <- Im[,2]
  
  #AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  #NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  #PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical)
  
 # PosteriorProbability <- BLBF_CrossValidateKDEall(ImAF = ImAF , ImNAF = ImNAF , AFLogical = AFLogical  , PriorProbability = PriorProbability)

#  return(PosteriorProbability)
  
#}
BLBF_CalculatePosteriorProbabilityBLA <- function(DataStructure , AFLogical , PriorProbability , StartIndex , NumberofPoints , VariablesToView){
 
  test <- BLBF_DoubleCrossValidateInference(DataStructure ,AFLogical, StartIndex , NumberofPoints , VariablesToView  )
  test[,1]= (test[,1] - mean(test[AFLogical,1]))/sqrt(var(test[AFLogical,1]))
  test[,2]= (test[,2] - mean(test[AFLogical ==0,2]))/sqrt(var(test[AFLogical==0,2]))
  ImAF   <-  test[,1]
  ImNAF  <-  test[,2]
  
  
  ImplausibleForBoth <- ((test[,1] > 2)*(test[,2] > 2)) == 0
  PosteriorProbability <- 0*ImAF
  PosteriorProbability[ImplausibleForBoth] <- BLBF_CrossValidateKDEall(ImAF = ImAF[ImplausibleForBoth] , ImNAF = ImNAF[ImplausibleForBoth] , AFLogical = AFLogical[ImplausibleForBoth]  , PriorProbability = PriorProbability[ImplausibleForBoth])
  PosteriorProbability[ImplausibleForBoth] <- PriorProbability[ImplausibleForBoth == 0]
    
  return(PosteriorProbability)
}

BLBF_CalculatePosteriorProbabilityBLA2 <- function(DataStructure , AFLogical , PriorProbability , StartIndex , NumberofPoints , VariablesToView){
  LengthVector <- length(VariablesToView)
  tmpdata <- DataStructure[ , VariablesToView, c(StartIndex:(StartIndex+NumberofPoints) )]
  AFSOS <- BLBF_CalculateSecondOrderStruct( tmpdata[which(AFLogical == T) , , ])
  NAFSOS <- BLBF_CalculateSecondOrderStruct( tmpdata[which(AFLogical == F) , , ])
  
  ImAF <- apply(tmpdata , 1 , function(X){BLBF_CrossValidateCalculateIm(tmpdata[AFLogical == T, ,  ] , X , LengthVector  , NumberofPoints )})
  ImNAF <- apply(tmpdata , 1 , function(X){BLBF_CrossValidateCalculateIm(tmpdata[AFLogical == F, ,  ] , X , LengthVector  , NumberofPoints )})
  ImAF[AFLogical == T ] <- BLBF_CrossValidateCalculateImAll(tmpdata[AFLogical == T, ,  ] , LengthVector , NumberofPoints )
  ImNAF[AFLogical == F ] <- BLBF_CrossValidateCalculateImAll(tmpdata[AFLogical == F, ,  ] , LengthVector , NumberofPoints ) 
  
  #AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  #NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  #PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical , alpha = PriorProbability)
  
  PosteriorProbability <- BLBF_CrossValidateKDEall(ImAF = ImAF , ImNAF = ImNAF , AFLogical = AFLogical  , PriorProbability = PriorProbability)
  return(PosteriorProbability)
}
BLBF_CrossValidateKDE <- function(ImAF , ImNAF , ImAF_test , ImNAF_test , AFLogical , PriorProbability = c(sum(AFLogical)/length(AFLogical))){
  
  AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  #AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  #NAFModel <- kde(cbind(ImAF[((AFLogical ==F)*(ImAF < 3))==1] , ImNAF[((AFLogical ==F)*(AFLogical ==F)*(ImAF < 3))==1]) )
  
  PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel = AFModel, NAFModel = NAFModel , ImMatrix = cbind(ImAF_test , ImNAF_test) , AFLogical = AFLogical , alpha = PriorProbability)
  return( PosteriorProbability )
}
BLBF_CrossValidateKDEall <- function(ImAF , ImNAF , AFLogical , PriorProbability = c(sum(AFLogical)/length(AFLogical))){

  PosteriorProbability <- matrix(0 , length(AFLogical) , 1)
  NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  if(length(PriorProbability)  == 1){ 
  PosteriorProbability[AFLogical == F] <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF[AFLogical == F] , ImNAF[AFLogical == F]) , AFLogical[AFLogical == F] , alpha = PriorProbability)}else{
  PosteriorProbability[AFLogical == F] <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF[AFLogical == F] , ImNAF[AFLogical == F]) , AFLogical[AFLogical == F] , alpha = PriorProbability[AFLogical == F])}
  
  for(kk in 1:length(AFLogical)){
  if(AFLogical[kk] == 1){  
    if(length(PriorProbability)  == 1){
    PosteriorProbability[kk] <-  BLBF_CrossValidateKDE(ImAF[-kk] , ImNAF[-kk] , ImAF[kk] , ImNAF[kk]  , AFLogical[-kk] , PriorProbability)}
    if(length(PriorProbability)  > 1){
      PosteriorProbability[kk] <-  BLBF_CrossValidateKDE(ImAF = ImAF[-kk] , ImNAF = ImNAF[-kk] ,ImAF_test =  ImAF[kk] ,ImNAF_test =  ImNAF[kk]  , AFLogical = AFLogical[-kk] ,PriorProbability =  PriorProbability[kk])}
    DP_WaitBar(kk/length(AFLogical))
    }
  }
  return(PosteriorProbability)
}  

BLBF_CrossValidateCalculateAdjustedVersionSOS <- function(AdjustedVersions ){
  Discrepancy <- matrix(0 , dim(AdjustedVersions)[1] , 1)
  
  for(kk in 1:dim(Discrepancy)[1]){ 
    AdjustedVersionsSOS <- BLU_CreateSOS(mu = apply(AdjustedVersions[-kk,] , 2 , mean) , Sigma = cov(AdjustedVersions[-kk ,]) )
    Discrepancy[kk] <-  BLBF_Discrepancy(AdjustedVersionsSOS , AdjustedVersions[kk , ] )
  }  
  return(Discrepancy)
}  

BLBF_DoubleCrossValidateInference <- function(DataStructure ,AFLogical, StartIndex , NumberofPoints , VariablesToView  ){
  LengthVector <- length(VariablesToView)
  Data <- DataStructure[ , VariablesToView, c(StartIndex:(StartIndex+NumberofPoints) )]  
  ImStruct <- matrix(0 , length(AFLogical) , 2)
  
  for(kk in 1:length(AFLogical)){
    
    tmpData <- Data[-kk , , ] 
    tmpAFLogical <- AFLogical[-kk]
    
    AFSOS <-  BLBF_CalculateSecondOrderStruct( tmpData[tmpAFLogical , , ])
    NAFSOS <- BLBF_CalculateSecondOrderStruct( tmpData[tmpAFLogical ==F , , ])
    
    AdjustedVersionsAF <- t(apply( Data , 1 , function(X){BLBF_CalculateForecastAdjustedVersion(AFSOS , X[1:(LengthVector*(NumberofPoints))] , X[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints + 1) )])}))
    AdjustedVersionsNAF <- t(apply( Data , 1 , function(X){BLBF_CalculateForecastAdjustedVersion(NAFSOS , X[1:(LengthVector*(NumberofPoints))] , X[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints +1) )])}))
    AdjustedVersionsAFtest <- AdjustedVersionsAF[kk , ]
    AdjustedVersionsNAFtest <- AdjustedVersionsNAF[kk , ] 
    
    AdjustedVersionsAF <- AdjustedVersionsAF[-kk , ]
    AdjustedVersionsNAF <- AdjustedVersionsNAF[-kk , ] 
    
    AdjustedVersionsAF[tmpAFLogical == T] <- BLBF_CrossValidateAdjustedVersionAll(tmpData[tmpAFLogical == T , , ]  , LengthVector , NumberofPoints )
    AdjustedVersionsNAF[tmpAFLogical == F] <- BLBF_CrossValidateAdjustedVersionAll(tmpData[tmpAFLogical == F , , ]  , LengthVector , NumberofPoints )
    
    AdjustedVersionsSOSAF <- BLU_CreateSOS(mu = apply(AdjustedVersionsAF[tmpAFLogical,] , 2 , mean) , Sigma = cov(AdjustedVersionsNAF[tmpAFLogical == F,]) )
    AdjustedVersionsSOSNAF <- BLU_CreateSOS(mu = apply(AdjustedVersionsNAF[tmpAFLogical == F,]  , 2 , mean) , Sigma = cov(AdjustedVersionsNAF[tmpAFLogical == F,]) )
    
    ImStruct[ kk , 1] <- BLBF_Discrepancy(AdjustedVersionsSOSAF , AdjustedVersionsAFtest )
    ImStruct[ kk , 2] <- BLBF_Discrepancy(AdjustedVersionsSOSAF , AdjustedVersionsNAFtest )
    DP_WaitBar(kk/length(AFLogical))
  }
  
  return(ImStruct)
  
}


##### Bayes linear Bayes classifier ######
BLBC_FitBayesLinearBayesClassifier <- function(Data , Labels , PriorProbability =  sum(Labels)/length(Labels)){

  Im <- BLBL_CrossValidateImplausabilitiesALL(TrainingData = Data , Labels)
  ImAF <- Im[,1]
  ImNAF <- Im[,2]
  
  ImAF <- (ImAF - mean(ImAF[Labels == 1]) ) / sqrt(var(ImAF[Labels == 1]))
  ImNAF <- (ImNAF - mean(ImNAF[Labels == 0])) / sqrt(var(ImNAF[Labels == 0]))
  
  
  x11()
  par(mfrow = c(1 , 2))
  print(BC_PlotCompareSingleHists((ImAF[Labels ==F]) , (ImAF[Labels ==T]) ))
  print(BC_PlotCompareSingleHists((ImNAF[Labels ==F]) , (ImNAF[Labels ==T]) ))
  
  AFModel <- kde(cbind(ImAF[Labels == 1 ] , ImNAF[Labels == 1 ]) )
  NAFModel <- kde(cbind(ImAF[Labels == 0 ] , ImNAF[Labels == 0 ]))   
  
  x11(20,14)
  plot(NAFModel  , col ='blue' , xlab =c('Implausability AF') , ylab =c('Implausability NAF'))
  plot(AFModel, col = 'red' , add = T)
  points(ImAF[Labels == 0] , ImNAF[Labels == 0] , col = rgb(0 , 0 , 1 , alpha = 0.2) , pch = 16 )
  points(ImAF[Labels == 1] , ImNAF[Labels == 1], col = rgb(1 , 0 , 0 , alpha = 0.3), pch = 16)
  abline(v = 3)
  abline(h = 3)
  
  
  PosteriorProbability <- BLBF_CrossValidateKDEall(ImAF = ImAF ,ImNAF =  ImNAF  ,AFLogical =  (Labels == 1) , PriorProbability = PriorProbability  )
  PosteriorProbability[is.na(PosteriorProbability)] <- PriorProbability
  
  PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbability , Labels == 1))
  ROCplot2 <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1)
  NPVPPVPlot2 <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept =( 1-PriorProbability)) + geom_hline(yintercept = PriorProbability)+  ggtitle('PPV vs NPV Curves')
  
  x11(20,14)
  grid.arrange(ROCplot2 , NPVPPVPlot2 + xlim(( 1-PriorProbability) , 1) , nrow= 2)
  
  x11(20,14)
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbability , Labels == 1), BinWidth = 0.1))
  print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))
  
  return(PosteriorProbability)
  
}
BLBL_CrossValidateImplausabilitiesALL <- function(TrainingData , AFLogical){
  
  Implausabilities <- matrix(0 , dim(TrainingData)[1] , 2)
 
  for(ii in 1:dim(Implausabilities)[1]){
    Implausabilities[ii , ] <- BLBL_CrossValidateImplausabilitiesSingle(TrainingData[-ii , ] , TrainingData[ii , ],AFLogical[-ii])
  }
   return(Implausabilities)
}
BLBL_CrossValidateImplausabilitiesSingle <- function(TrainingData , TestData ,AFLogical ){
  
  mAF <- (apply(TrainingData[AFLogical == T , ]  , 2 , mean))
  CAF <- DP_AddNugget(cov(TrainingData[AFLogical == T , ]) , 1e-6*diag(diag(cov(TrainingData[AFLogical == T , ]))))
  mNAF <- (apply(TrainingData[AFLogical == F , ]  , 2 , mean))
  CNAF <- DP_AddNugget(cov(TrainingData[AFLogical == F , ]) , 1e-6*diag(diag(cov(TrainingData[AFLogical == F , ]))))
  
  return(cbind(log(mahalanobis(TestData , mAF , CAF)) , log(mahalanobis(TestData , mNAF , CNAF))) )
}
BLBF_PlotSenSpecPerformanceSweep <- function(ProbabilityUpdateStruct , AFLogical){
  
  PerformanceSweepStruct <- array(0 , c(dim(ProbabilityUpdateStruct)[2],  101 , 4) )
  
  if(exists('alpha')){
  ProbabilityUpdateStruct[is.na(ProbabilityUpdateStruct)] <- alpha[1]
  }else{
    ProbabilityUpdateStruct[is.na(ProbabilityUpdateStruct)] <- 0
  }
  for(i in 1:dim(ProbabilityUpdateStruct)[2]){
    PerformanceSweepStruct[i , , ] <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(ProbabilityUpdateStruct[, i] , AFLogical == T))
  }
  
  x11()
  plot(0*PerformanceSweepStruct[1 ,  , 1] , ylim = c(0,1) , xlim = c(0,1) , col = rgb(1,1,1), xlab ='1 - Specificity' , ylab ='Sensitivity')
  for(i in 1:dim(ProbabilityUpdateStruct)[2]){
    lines(1-PerformanceSweepStruct[i ,  , 2] ,PerformanceSweepStruct[i ,  , 1]  , col = rgb(0 , 0 , 1 , alpha = 1/i))
  }
  abline(0,1)
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  
}

BLBF_PlotProbabilityTrajectories <- function(ProbabilityUpdateStruct , AFLogical){
  
  x11()
  plot(0*ProbabilityUpdateStruct[1 , ] , ylim = c(0,1) , xlim = c(0,30) , col = rgb(1,1,1) , xlab ='Time Before AF' , ylab ='P')
  for(i in 1:dim(ProbabilityUpdateStruct)[1]){
    if(AFLogical[i] == 1){
      lines(ProbabilityUpdateStruct[i , ] , col = rgb(1,0,0 , alpha = 0.4))
    }
    if(AFLogical[i] == 0){
      lines(ProbabilityUpdateStruct[i , ] ,col = rgb(0,0,1 , alpha = 0.1))
    }
  }
  
}

BLBF_CreateNamesPWaveVariables <- function(){
  return(c('PAmp' , 
           'PAMpLocation',
           'Pstart',
           'PEnd',
           'PWidth',
           'DiffBetweenAtria',
           'PRLevel',
           'max_dp',
           'min_dp',
           'E_V',
           'max_V'))
}
BLBF_CalculatePosteriorProbabilityNBLA<- function( DataStructure , AFLogical , VariablesToView = -c( 13) ,  indextoview = 371  , PriorProbability =  c(sum(AFLogical)/length(AFLogical)) ){
  
  Im <- BLBL_CrossValidateImplausabilitiesALL(TrainingData = DataStructure[, VariablesToView, indextoview] , AFLogical)
  
  
  #AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  #NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  #PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical)
  PosteriorProbability <- matrix(0 , dim(Im)[1] , 1)
  if(length(PriorProbability) == 1){
  PosteriorProbability[ ((Im[,1] > 3)*(Im[,2] > 3)) == 1 ] <- PriorProbability
  }else{
    PosteriorProbability[ ((Im[,1] > 3)*(Im[,2] > 3)) == 1 ] <- PriorProbability[ ((Im[,1] > 3)*(Im[,2] > 3)) == 1 ]
  }
  
  if(length(PriorProbability) == 1){
  PosteriorProbability[((Im[,1] > 3)*(Im[,2] > 3)) == 0] <- BLBF_CrossValidateKDEall(ImAF = Im[((Im[,1] > 3)*(Im[,2] > 3)) == 0,1] , ImNAF = Im[((Im[,1] > 3)*(Im[,2] > 3)) == 0,2] , AFLogical = AFLogical[((Im[,1] > 3)*(Im[,2] > 3)) == 0]  , PriorProbability = PriorProbability)
  }else{
    PosteriorProbability[((Im[,1] > 3)*(Im[,2] > 3)) == 0] <- BLBF_CrossValidateKDEall(ImAF = Im[((Im[,1] > 3)*(Im[,2] > 3)) == 0,1] , ImNAF = Im[((Im[,1] > 3)*(Im[,2] > 3)) == 0,2] , AFLogical = AFLogical[((Im[,1] > 3)*(Im[,2] > 3)) == 0]  , PriorProbability = PriorProbability[((Im[,1] > 3)*(Im[,2] > 3)) == 0])  
  }
  return(PosteriorProbability)
  
}