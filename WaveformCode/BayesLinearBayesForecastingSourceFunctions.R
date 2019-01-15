BLBF_CalculateDensities <- function(AFModel, NAFModel , ImMatrix , AFLogical){
  
  f_1 <- matrix(0 , length(AFLogical) , 2)
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
  X <- BLBF_ReshapetoMatrix(X)
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
  return(log((mahalanobis(x , SOS$Expectation , SOS$Covariance))))
}
BLBF_CalculatePosteriorProbabilityNBLA <- function( DataStructure , AFLogical , VariablesToView ,  indextoview ){
  
  mAF <- (apply(DataStructure[AFLogical == T , , indextoview]  , 2 , mean))
  CAF <- cov(DataStructure[AFLogical == T , , indextoview])
  mNAF <- (apply(DataStructure[AFLogical == F , , indextoview]  , 2 , mean))
  CNAF <- cov(DataStructure[AFLogical == F , , indextoview])
  ImAF <- BLBF_ForecastingCalculateImplausabilities(DataStructure[, , indextoview] , mAF , CAF)
  ImNAF <- BLBF_ForecastingCalculateImplausabilities(DataStructure[, , indextoview] , mNAF , CNAF)
  ImAF <- (ImAF - mean(ImAF[AFLogical ==T])) / sqrt(var(ImAF[AFLogical ==T])) 
  ImNAF <- (ImNAF - mean(ImNAF[AFLogical ==F])) / sqrt(var(ImNAF[AFLogical == F])) 
  
  AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical)
  return(PosteriorProbability)
  
}
BLBF_CalculatePosteriorProbabilityBLA <- function(DataStructure , AFLogical , PriorProbability , StartIndex , NumberofPoints , VariablesToView){
  LengthVector <- length(VariablesToView)
  tmpdata <- DataStructure[ , VariablesToView, c(StartIndex:(StartIndex+NumberofPoints) )]
  AFSOS <- BLBF_CalculateSecondOrderStruct( tmpdata[which(AFLogical)[sort(sample(1:sum(AFLogical) , sum(AFLogical) -1))] , , ])
  NAFSOS <- BLBF_CalculateSecondOrderStruct( tmpdata[which(AFLogical == F)[sort(sample(1:sum(AFLogical==F) , sum(AFLogical==F) -1))] , , ])
  
  AdjustedVersionsAF <- t(apply( tmpdata , 1 , function(X){BLBF_CalculateForecastAdjustedVersion(AFSOS , X[1:(LengthVector*(NumberofPoints))] , X[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints + 1) )])}))
  AdjustedVersionsNAF <- t(apply( tmpdata , 1 , function(X){BLBF_CalculateForecastAdjustedVersion(NAFSOS , X[1:(LengthVector*(NumberofPoints))] , X[((LengthVector*(NumberofPoints)) +1):(LengthVector*(NumberofPoints +1) )])}))
  
  AdjustedVersionsAF[AFLogical == T] <- BLBF_CrossValidateAdjustedVersionAll(tmpdata[AFLogical == T , , ]  , LengthVector , NumberofPoints )
  AdjustedVersionsNAF[AFLogical == F] <- BLBF_CrossValidateAdjustedVersionAll(tmpdata[AFLogical == F , , ]  , LengthVector , NumberofPoints )
  
  # Could be over fitting - check that it is not.
  AdjustedVersionsSOSAF <- BLU_CreateSOS(mu = apply(AdjustedVersionsAF[AFLogical == T ,] , 2 , mean) , Sigma = cov(AdjustedVersionsAF[AFLogical == T ,]) )
  AdjustedVersionsSOSNAF <- BLU_CreateSOS(mu = apply(AdjustedVersionsNAF[AFLogical == F ,]  , 2 , mean) , Sigma = cov(AdjustedVersionsNAF[AFLogical == F , ]) )
  
  ImAF <- apply(AdjustedVersionsAF , 1 , function(X){BLBF_Discrepancy(AdjustedVersionsSOSAF , X )})
  ImNAF <- apply(AdjustedVersionsNAF , 1 , function(X){BLBF_Discrepancy(AdjustedVersionsSOSNAF, X)} ) 
  ImAF <- (ImAF - mean(ImAF[AFLogical ==T])) / sqrt(var(ImAF[AFLogical ==T])) 
  ImNAF <- (ImNAF - mean(ImNAF[AFLogical ==F])) / sqrt(var(ImNAF[AFLogical == F])) 
  
  AFModel <- kde(cbind(ImAF[AFLogical ==T] , ImNAF[AFLogical ==T]) )
  NAFModel <- kde(cbind(ImAF[AFLogical ==F] , ImNAF[AFLogical ==F]) )
  
  PosteriorProbability <- BLBF_CalculatePosteriorProbabilities(AFModel, NAFModel , cbind(ImAF , ImNAF) , AFLogical , alpha = PriorProbability)
  return(PosteriorProbability)
}