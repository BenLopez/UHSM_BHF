
##### Kernel density estimation #####
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
  #inv_H = solve(H)
  
  distmatrix = mahalanobis(x = X , center = Xstar , cov = H)
  distmatrix = -0.5*( d*log(2*pi) + log(det_H) + distmatrix  )
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
KDE_HistoryMatchBandWidth <- function( Trainingset , Validationset , numberofsamples = 1000 , EmulatorSettings = BE_CreateDefaultEmulationClass() , HistoryMatchSettings = BE_CreateDefaultHistoryMatchClass() , PriorRange = c(0,0.1)){
  print('Calculating Second Order Specification for Diagnostics.')  
  SecondOrderImSpecifictaion <- KDE_CreateSecondOrderSpecificationImp(Trainingset = Trainingset 
                                                                      , Validationset = Validationset  
                                                                      , numberofsamples = numberofsamples
                                                                      , numberofrepitions = 1000)  
  print('Second Order Specification for Diagnostics Calculated.')  
  
  print('History Matching')
  {
    SizechiStariMinus1 <- 1
    WaveCounter <- 1
    chi_star <<- matrix(0 , 1 , 1)
    while( length(chi_star) >=  SizechiStariMinus1 ){
      
      SizechiStariMinus1 <- length(chi_star)
      
      tmp <- BE_CalulateNonImplausibleSets(EmulatorSettings = EmulatorSettings , 
                                           HistoryMatchSettings = HistoryMatchSettings , 
                                           PriorRange = PriorRange,
                                           chi_star = chi_star , WaveCounter)
      chi_star <<- tmp$chi_star
      EmulatorSettings <- tmp$EmulatorSettings
      
      
      PriorRange <- DP_CalculateLimits(chi_star)
      print(paste0('Number of Non-Implausible Points = ', length(chi_star)))
      print(paste0('Wave ' , WaveCounter , 'completed.'))
      abline(v = PriorRange[1])
      abline(v = PriorRange[2])
      
      if(is.na(chi_star[1])){
        break
      }
      
      WaveCounter <- WaveCounter + 1
      
    }
    HistoryMatchSettings$AddPointsToEmulator = 1
    chi_star <<- sample( chi_star , SizechiStariMinus1 , replace = TRUE)
    
    SizechiStariMinus1 <- 10
    while(var(chi_star) <=  SizechiStariMinus1){
      
      SizechiStariMinus1 <- var(chi_star)
      
      tmp <- BE_CalulateNonImplausibleSets(EmulatorSettings = EmulatorSettings , 
                                           HistoryMatchSettings = HistoryMatchSettings , 
                                           PriorRange = PriorRange,
                                           chi_star = chi_star , WaveCounter )
      chi_star <<- tmp$chi_star
      EmulatorSettings <- tmp$EmulatorSettings
      
      PriorRange <- DP_CalculateLimits(chi_star)
      print(paste0('Number of Non-Implausible Points = ', length(chi_star)))
      print(paste0('Wave ' , WaveCounter , 'completed.'))
      WaveCounter <- WaveCounter + 1
      abline(v = PriorRange[1])
      abline(v = PriorRange[2])
      
      
      if(is.na(chi_star[1])){
        break
      }
      
    }
  }
  
  print('History Match Complete.')
  return(chi_star)
}
##### pusedo Kernel density estimation #####

##### One Dimensional Kernel Density Estimation #####
KDE_OneDDenistyEstimateBatch <- function(X , Xstar , H){
  X <- as.matrix(X)
  Xstar <- as.matrix(Xstar)
  
  return(apply(Xstar , 1 , function(Y){ KDE_OneDDenistyEstimate(X , Y , H ) } ))
}
KDE_OneDDenistyEstimate <- function(X , Xstar , H){
  X <- as.matrix(X)
  Xstar <- as.matrix(Xstar)
  
  f_i <- apply(X , 1 , function(Y){dnorm(x = Xstar , mean = Y , sd = sqrt(H) ) })
  
  return(mean(f_i))
}