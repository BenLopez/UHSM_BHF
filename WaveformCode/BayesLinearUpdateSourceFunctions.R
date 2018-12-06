BLU_CreateSOS <- function(mu = 0 , Sigma = 1  , WeightMatrix  = 1){
  return(setNames(list(mu,Sigma , 1) , c('Expectation' , 'Covariance' , 'WeightMatrix')))
}
BLU_VectoriseData <- function(X){
  
  if( length(dim(X)) <= 2 ){
    return(X)
  }
  if( length(dim(X)) > 2 ){
    output <- matrix(0 , dim(X)[1] ,  dim(X)[2]*dim(X)[3]) 
    for(i in 1:dim(X)[3]){
      output[1:dim(X)[1], seq(i , (dim(X)[2]*dim(X)[3]) - (dim(X)[3] - i) , dim(X)[3] ) ] <- X[1:dim(X)[1] ,1:dim(X)[2] , i]
    }
    return(output)
  }
}
BLU_Chol <- function(X){
  return(t(chol(X)))
}
BLU_CreateWeightMatrix <- function( SOS , nugget = 0.000000001 ){
  
  dimSOS <- dim( SOS$Covariance )
  L <- BLU_Chol( DP_AddNugget(SOS$Covariance , nugget))
  
  Weightmatrix = array(0 , c(dimSOS) )
  for(ii in 2:dimSOS[1]){
    covXD = as.matrix(SOS$Covariance[ ii , 1:(ii-1) ])
    Weightmatrix[1:(ii-1) , ii] <- t(covXD)%*%backsolve(t(L[1:(ii-1) , 1:(ii-1)]) , backsolve((L[1:(ii-1) , 1:(ii-1)]) , diag(ii-1)  , upper.tri = FALSE), upper.tri = TRUE)
  #DP_WaitBar(ii/dimSOS[1])
  }
  
  return(Weightmatrix)
}
BLU_CalculateAdjustedExpectation <- function( SOS , D ,  variableindex){
  if(!is.matrix(SOS$WeightMatrix) == 1){
    SOS$WeightMatrix <- BLU_CreateWeightMatrix(SOS)
  }
  
  diff <- D - SOS$Expectation
  AdjustedExpectation <- SOS$Expectation +  t(SOS$WeightMatrix)%*%diff
  return(AdjustedExpectation)         
}
BLU_CalulateAdjustedCovariance <- function(AdjustedVersions , variableindex , rangeoftimes){
  
  if( length(variableindex) > 1){
  AdjustedCovariance <- array(0 , c(length(rangeoftimes)  , length(variableindex) , length(variableindex) )  )
  
  for(ii in 1:length(rangeoftimes)){
    AdjustedCovariance[ii , , ] <- cov(t(AdjustedVersions[ ((length(variableindex)*(ii - 1)) + 1) :(length(variableindex)*ii), ]) )  
  }  }
  if( length(variableindex) == 1){
    AdjustedCovariance <- array(0 , c(length(rangeoftimes)  , length(variableindex) )  )
    for(ii in 1:length(rangeoftimes)){
    AdjustedCovariance[ii , ] <- var(AdjustedVersions[ ((length(variableindex)*(ii - 1)) + 1) :(length(variableindex)*ii), ] )  
    }                                 
  }
  return(AdjustedCovariance)
}
BLU_ReshapeVectorToMatrix <- function(D, rangeoftimes , variableindex){
  output <- matrix( 0, length(rangeoftimes) , length(variableindex) )
  for(i in 1:length(variableindex)){
    output[  , i ] <- D[seq(i  , length(rangeoftimes)*length(variableindex) - (length(variableindex) - i) , length(variableindex) )  ]
  }
  return(output)
}
BLU_CalulateMultivariateDiscrepancyTimeSeries <- function(D , rangeoftimes ,variableindex , AdjustedCovariance ){
  if(length(variableindex)>1 ){
  D2 <- BLU_ReshapeVectorToMatrix( D , rangeoftimes , variableindex)
  Discrepancies <- matrix(0 , dim(D2)[1] , 1)
  for(ii in 1:length(rangeoftimes)){
    Discrepancies[ ii , 1 ] <- t(D2[ii , ])%*%solve(AdjustedCovariance[ii,,])%*%D2[ii , ]
  }
  }
  if(length(variableindex) == 1 ){
    D2 <- BLU_ReshapeVectorToMatrix( D , rangeoftimes , variableindex)
    Discrepancies <- matrix(0 , dim(D2)[1] , 1)
    for(ii in 1:length(rangeoftimes)){
      Discrepancies[ ii , 1 ] <- t(D2[ii , ])%*%solve(AdjustedCovariance[ii,])%*%D2[ii , ]
    }
  }
  return(Discrepancies)
}
