KDE_GaussianKernel <- function(X , Xstar , H){
  
  
  # Expenential family correlation length family
  if(dim(as.matrix(Xstar))[1] != 1){
  if(dim(as.matrix(Xstar))[1] != 1 & dim(as.matrix(Xstar))[1] == 1){
    Xstar <- t(as.matrix(Xstar))
  }else{
    stop('Xstar must be 1xd')
  }  
  }
  
  # Turn into matrices
  d <- dim(X)[2]
  det_H = det(H)
  inv_H = solve(H)
  
  distmatrix = as.matrix(apply(X , 1 , function(Y){(Y - Xstar)%*%inv_H%*%t(Y - Xstar)}))
  distmatrix = log((2*pi)^(-d/2)) + log(det_H) + (-0.5*distmatrix)
  return(exp(distmatrix))
} 
KDE_DensityEstimate <- function(X , Xstar , H , KernelFUN = KDE_GaussianKernel){
  
  return( mean(KernelFUN(X , Xstar , H)) )
  
}
KDE_RuleofThumbH <- function(X){
  d <- size(X)[2]
  n <- size(X)[1]
  
  V <- apply(X , 2 , var)
  H <-  diag(((4/(d+2))^(1/(d+4)))*(n^(-1/(d+4)))*sqrt(V))
  return(H)
}
KDE_CalulateHistImplausability <- function(X , X2){
  
  output <- matrix(0 , size(X)[2] , 1)
  
  for(variabletoview in c(1:size(X)[2])){
    
    tmp  <-  hist(rbind(X[ , variabletoview] , X2[ , variabletoview]) , breaks = 20  , plot = FALSE)
    tmp2 <-  hist(X[ , variabletoview]  , breaks = tmp$breaks, plot = FALSE)
    tmp3 <-  hist(X2[ , variabletoview]  , breaks = tmp$breaks, plot = FALSE)
    
    output[variabletoview] <- sum(abs(tmp2$density - tmp3$density))
    
  }  
  return(output)
}
KDE_CalulateMeanMinDistances <- function( X , X2 ){

  if((size(X)[1]*size(X)[2]) >  (10000^2)){
    stop('A Matrix is too large, sample before appluing function.')
  } 
  
  V_sample <- sqrt(apply(X , 2 , var))
  tmpX <- t(apply(X , 1 , function(Y){Y/V_sample} ))
  tmpX2 <- t(apply(X2 , 1 , function(Y){Y/V_sample} ))
  
 
  return(mean(apply(abs(as.matrix(pdist(tmpX , tmpX2))) , 2 , min)))
}
KDE_CalulateJointImplausability <- function(Trainingset , validationsample){
# Wrapper function
  return(rbind(KDE_CalulateHistImplausability(Trainingset , validationsample) , KDE_CalulateMeanMinDistances(Trainingset , validationsample)) )
}
KDE_CreateSecondOrderSpecificationImp <- function(Trainingset , Validationset ,  numberofsamples = 1000 , numberofrepitions = 100 ){
  
  Immatrix <- matrix(0 , numberofrepitions ,  size(Trainingset)[2] +1  )
  
  for( ii in 1:numberofrepitions ){
    
    sampleindexes <- sample( 1:dim(Validationset)[1] , numberofsamples )
    validationsample <- Validationset[sampleindexes , ]
    
    Immatrix[ii , 1:( size(Trainingset)[2] + 1 )] <- KDE_CalulateJointImplausability( Trainingset , validationsample) 
    
  }
  output <- list( mu = apply(Immatrix , 2 , mean) , Sigma = cov(Immatrix) , Inv_Sigma = solve(cov(Immatrix)) )
  return(output)
}
KDE_CalulateImplausabilityVectorforKDE <- function( Trainingset , MclustDistributionStruct , numberofsamples = 1000 , numberofrepitions = 1 ){
  
  Immatrix <- matrix(0 , numberofrepitions ,  size(Trainingset)[2] +1  )
  
  for( ii in 1:numberofrepitions ){
    ValidationSample <- BC_SampleGMM(MclustDistributionStruct ,  numberofsamples )
    Immatrix[ii , 1:( size(Trainingset)[2] + 1 )] <- KDE_CalulateJointImplausability(Trainingset , ValidationSample) 
  }
  return(apply(Immatrix , 2 , mean ))  
}
KDE_CreateMclustClassFromSample <- function(X , H ){
  MclustDistributionStruct <- BC_CreateDefaultmclustStruct() 
  MclustDistributionStruct$parameters$pro <- rep(1 , size(X)[1])/size(X)[1]
  MclustDistributionStruct$parameters$mean <- t(X)
  MclustDistributionStruct$parameters$variance$sigma <- array(H , c( size(X)[2] ,  size(X)[2] , size(X)[1]) )
  return( MclustDistributionStruct )  
}
KDE_MaxIm <- function(Im , SecondOrderImSpecifictaion){
  return(max(abs((Im - SecondOrderImSpecifictaion$mu) /  sqrt(diag(SecondOrderImSpecifictaion$Sigma)))))
}
KDE_JointIm <- function(Im , SecondOrderImSpecifictaion){
  diff <- as.matrix(Im - SecondOrderImSpecifictaion$mu)
  return( t(diff)%*%SecondOrderImSpecifictaion$Inv_Sigma%*%diff )
}
