CD_CalculateActualUpdatedProbabilities <- function(SampleGP , alpha , KXX , mu1 , mu2 , v1 , v2 ){
  ActualProbabilities <- matrix(0, dim(SampleGP)[1] + 1 , 1)
  ActualProbabilities[1 , ] <- alpha
  
  for(jj in 2:(dim(ActualProbabilities)[1])){
    ii = jj - 1
    if(ii == 1){
      f_1 = dnorm(SampleGP[ii] , mean = mu1 , sd = sqrt( v1 ))
      f_2 = dnorm(SampleGP[ii] , mean = mu2 , sd = sqrt( v2 ))
    }
    if(ii > 1){
      # Recursive conditioning
      invD <- solve(KXX[1:(ii-1) , 1:(ii-1) ] + 0.000000001*diag(dim(as.matrix(KXX[1:(ii-1) ,1:(ii-1) ]))[1]) )
      E_D1 <- mu1 + (rev(KXX[2:(ii) , 1]) )%*%invD%*%(SampleGP[1:(ii-1)] - mu1 )
      V_D1 <- v1*(1  - (rev(KXX[2:(ii) , 1]) )%*%invD%*%(rev(KXX[2:(ii) , 1]) ))  
      
      E_D2 <- mu2 + (rev(KXX[2:(ii) , 1]) )%*%invD%*%(SampleGP[1:(ii-1)] - mu2 )
      V_D2 <- v2*(1  - (rev(KXX[2:(ii) , 1]) )%*%invD%*%(rev(KXX[2:(ii) , 1]) ))  
      
      f_1 = dnorm( SampleGP[ii] , mean = E_D1 , sd = sqrt(V_D1) )
      f_2 = dnorm( SampleGP[ii] , mean = E_D2 , sd = sqrt(V_D2) )
    }
    ActualProbabilities[jj , 1] <-  (ActualProbabilities[jj-1 , 1]*f_1) / ((1 - ActualProbabilities[jj-1 , 1])*f_2 + ActualProbabilities[jj-1 , 1]*f_1)
    #DP_WaitBar(ii/(dim(ActualProbabilities)[1]))
  }
  return(ActualProbabilities)
}
CD_CalculateIndividualDensities <- function(SampleGP ,mu1 , mu2 , v1,v2 ){
  SampleGP <- as.matrix(SampleGP)
  f_i = matrix( 0, dim(SampleGP)[1] , 2 )
  f_i[ , 1] <- dnorm( SampleGP , mean = mu1  , sd = sqrt(v1) )
  f_i[ , 2] <- dnorm( SampleGP , mean = mu2  , sd = sqrt(v2) )
  return(f_i)
}
CD_CalculateIndividualDensitiesMOPG <- function( SampleGP , mu1 , mu2 , Sigma1 , Sigma2 ){
  SampleGP <- as.matrix(SampleGP)
  f_i = matrix( 0, dim(SampleGP)[1] , 2 )
  f_i[ , 1] <- dmvnorm( SampleGP , mean = mu1  ,  sigma =  Sigma1 )
  f_i[ , 2] <- dmvnorm( SampleGP , mean = mu2  ,  sigma =  Sigma2 )
  return( f_i )
}
CD_CalulateIndenpendentEstimatedProbabilities <- function(SampleGP , alpha ,mu1 , mu2 , v1,v2 ){
  
  f_i = CD_CalculateIndividualDensities(SampleGP , mu1 , mu2 , v1 , v2 )
  
  JointProbabilties <- matrix(0, dim(SampleGP)[1] + 1 , 1)
  JointProbabilties[1, 1] = alpha 
  
  for(ii in 2:dim(JointProbabilties)[1]){
    JointProbabilties[ii , 1] <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])
  }
  return(JointProbabilties)
  
}
CD_CalulateIndenpendentEstimatedProbabilitiesMO <- function(SampleGP , alpha , mu1 , mu2 , v1 , v2 , weight = 0 ){
  
  f_i = CD_CalculateIndividualDensitiesMOPG(SampleGP , mu1 , mu2 , v1 , v2 )
  
  JointProbabilties <- matrix(0, dim(SampleGP)[1] + 1 , 1)
  JointProbabilties[1, 1] = alpha 
  
  for(ii in 2:dim(JointProbabilties)[1]){
    JointProbabilties[ii , 1] <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])
    if(weight != 0){
    JointProbabilties[ii , 1] <-  weight*JointProbabilties[ii-1 , 1] + (1-weight)*JointProbabilties[ii , 1] 
    }
    }
  return(JointProbabilties)
  
}
CD_CreateSecondorderSpecifiction <- function(X , Y){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  samplecovariance = cov(DP_RemoveNaRows(cbind(X , Y)))
  samplemean = colMeans(DP_RemoveNaRows(cbind(X , Y)) )
  
  
  output <- setNames(list( as.matrix(samplemean[1:dim(X)[2]]) ,
                           as.matrix(samplemean[(dim(X)[2] + 1):length(samplemean)]) ,
                           as.matrix(samplecovariance[1:dim(X)[2] , 1:dim(X)[2]]) ,
                           as.matrix(samplecovariance[(dim(X)[2]+1):dim(samplecovariance)[2] ,(dim(X)[2]+1):dim(samplecovariance)[2] ]) ,
                           as.matrix(samplecovariance[1:dim(X)[2] , (dim(X)[2]+1):dim(samplecovariance)[2]]  )) ,
                     c('E_X' , 'E_Y' , 'V_X' , 'V_Y' , 'Cov_XY') )
  return(output)
}
CD_BayesLinearAdjustmentSingleObs <- function(SOS , Y ){
  
  invV_D <- solve(SOS$V_Y)
  E_y_X = SOS$E_X + t(SOS$Cov_XY)%*%invV_D%*%( Y - SOS$E_Y )
  V_y_X = SOS$V_X - t(SOS$Cov_XY)%*%invV_D%*%(SOS$Cov_XY)
  
  output <- setNames( list( E_y_X , V_y_X ) , c('E_y_X' , 'V_y_X') )
  return( output )
  
}
CD_CalculateAdjustedUpdate <- function(AdjustedBeliefsStruct , SampleGP , mu1 , mu2 , v1 , v2 ){
  SampleProbabilities <- BE_SampleLHSinab( a = (AdjustedBeliefsStruct$E_y_X - 2*sqrt(AdjustedBeliefsStruct$V_y_X)) , b = (AdjustedBeliefsStruct$E_y_X + 2*sqrt(AdjustedBeliefsStruct$V_y_X)) , numbersamples = 100)
  f_1 = dnorm( SampleGP , mean = mu1 , sd = sqrt( v1 ) )
  f_2 = dnorm( SampleGP , mean = mu2 , sd = sqrt( v2 ) )
  
  SampleProbabilities <- apply(SampleProbabilities , 1 , function(X){ (f_1*X)/ (f_1*X   +   f_2*(1-X)) } )
  
  output <- setNames( list( mean(SampleProbabilities) , var(SampleProbabilities))  , c('E_B' , 'V_B') )
  return(output)  
}
CD_CalculateActualUpdatedProbabilitiesMOGP <- function(SampleGP , alpha = 0.1 , KXX , invD = 0 , mu1 , mu2 , Sigma1 , Sigma2 ){
  Probabilities <- matrix(0, dim(SampleGP)[1] + 1 , 1)
  Probabilities[1 , ] <- alpha
  
  E_D1 <- matrix(0, dim(SampleGP)[1]  , length(mu1))
  V_D1 <- array(0, c(dim(SampleGP)[1]  , length(mu1),length(mu1) ))
  E_D2 <- matrix(0, dim(SampleGP)[1]  , length(mu1))
  V_D2 <- array(0, c(dim(SampleGP)[1]  , length(mu1),length(mu1) ) )
  mahal  <- matrix(0, dim(SampleGP)[1]  , 2)
    f_1  <- matrix(0, dim(SampleGP)[1]  , 1)
  f_2  <- matrix(0, dim(SampleGP)[1]  , 1)
  if(is.list(invD) == FALSE){
    invD = CD_CalculateInverseVarDStack(KXX)
  }
  
  for(jj in 2:(dim(ActualProbabilities)[1])){
    ii = jj - 1
    if(ii == 1){
      f_1[ii,] = dmvnorm(x = SampleGP[ii,] , mean = mu1 , sigma = Sigma1)
      f_2[ii,] = dmvnorm(x = SampleGP[ii,] , mean = mu2 , sigma = Sigma2)
    }
    if(ii == 2){
      # Maniplulation to deal with R's horror handeling of matrices.
      tmp <- invD[[ii]]
      E_D1[ii,] <- t(mu1) + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1),] - t(mu1) )
      V_D1[ii,,] <- kronecker(Sigma1 , (1  - ( rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) )))  
      
      E_D2[ii,] <- t(mu2) + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1),] - t(mu2) )
      V_D2[ii,,] <- kronecker(Sigma2,(1  - (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) )))
      
      f_1[ii,] = dmvnorm( SampleGP[ii,] , mean = t(E_D1[ii,]) , sigma  = V_D1[ii,,] )
      f_2[ii,] = dmvnorm( SampleGP[ii,] , mean = t(E_D2[ii,]) , sigma =  V_D2[ii,,] ) 
      mahal[ii,1] = t(E_D1[ii,] - SampleGP[ii,])%*%solve(V_D1[ii,,])%*%(E_D1[ii,] - SampleGP[ii,])
      mahal[ii,2] = t(E_D2[ii,] - SampleGP[ii,])%*%solve(V_D2[ii,,])%*%(E_D2[ii,] - SampleGP[ii,])
      
            }
    
    if(ii > 2){
      # Recursive conditioning
      # Gaussian process update.
      tmp <- invD[[ii]]
      E_D1[ii,] <- t(mu1) + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1),] - t(matrix(mu1 , dim(SampleGP[1:(ii-1),])[2] ,dim(SampleGP[1:(ii-1),])[1] )) )
      V_D1[ii,,] <- kronecker(Sigma1 , (1  - ( rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) )))  
      
      E_D2[ii,] <- t(mu2) + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1),] - t(matrix(mu2 , dim(SampleGP[1:(ii-1),])[2] ,dim(SampleGP[1:(ii-1),])[1] )) )
      V_D2[ii,,] <- kronecker(Sigma2 , (1  - (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) )))
      
      f_1[ii,] = dmvnorm( SampleGP[ii,] , mean = t(E_D1[ii,]) , sigma  = V_D1[ii,,] )
      f_2[ii,] = dmvnorm( SampleGP[ii,] , mean = t(E_D2[ii,]) , sigma =  V_D2[ii,,] )
      
      mahal[ii,1] = t(E_D1[ii,] - SampleGP[ii,])%*%solve(V_D1[ii,,])%*%(E_D1[ii,] - SampleGP[ii,])
      mahal[ii,2] = t(E_D2[ii,] - SampleGP[ii,])%*%solve(V_D2[ii,,])%*%(E_D2[ii,] - SampleGP[ii,])
      
      }
    Probabilities[jj , 1] <-  (Probabilities[jj-1 , 1]*f_1[ii,]) / ((1 - Probabilities[jj-1 , 1])*f_2[ii,] + Probabilities[jj-1 , 1]*f_1[ii,])
    #DP_WaitBar(ii/(dim(ActualProbabilities)[1]))
  }
  return(Probabilities)
}
CD_CreateDefaultSpecification1D <- function(l =10 , p=2 , mu=0.02 , v=0.4 , n=500){
  
  Specification <- setNames( list(l , p , mu , v , n) , c('l' , 'p' , 'mu' , 'v' , 'n') )
  return(Specification)
  
}
CD_SampleDataFromSpecification1D <- function( Specification = CD_CreateDefaultSpecification1D()  ,  numberofsamples = 100 ){
  
  X <- matrix( c( 1:Specification$n ) ,Specification$n , 1)
  KXX <- CF_ExponentialFamily(X , X , Specification$l , Specification$p )
  SetofSamples <- matrix(0 , Specification$n , numberofsamples )
  
  for( kk in 1:numberofsamples ){
    SetofSamples[ , kk] <-  Specification$mu + BE_SampleGP( Specification$v*KXX[1:Specification$n , 1:Specification$n] )
  }
  return(SetofSamples)
}
CD_SampleDataFromSpecification1DMO <- function( Specification = CD_CreateDefaultSpecification1D()  ,  numberofsamples = 100 ){
  
  X <- matrix( c( 1:Specification$n ) ,Specification$n , 1)
  KXX <- CF_ExponentialFamily(X , X , Specification$l , Specification$p )
  SetofSamples <- array(0 , c(Specification$n , length(Specification$mu) ,  numberofsamples) )
  
  for( kk in 1:numberofsamples ){
    SetofSamples[ , , kk] <-  t(apply(BE_SampleSeparableMVGP(KXX , Specification$v)   , 1 , function(X){X + Specification$mu}))
  }
  return(SetofSamples)
}
CD_SampleDataFromSpecificationTP <- function( Specification = CD_CreateDefaultSpecification1D()  ,  numberofsamples = 100 , df = 5 ){
  
  X <- matrix( c( 1:Specification$n ) ,Specification$n , 1)
  KXX <- CF_ExponentialFamily(X , X , Specification$l , Specification$p )
  SetofSamples <-  Specification$mu + t(rmvt(n = numberofsamples , sigma = ((df - 2)/(df))*(Specification$v*KXX + 0.0000000000001 * diag(dim(KXX)[1])) , df = df))
  
  return(SetofSamples)
}
CD_CalculateUpdatedProbabilitiesFromSpecifications <- function(SetofSamples , Specification1 , Specification2 , alpha = 0.1){
  X <- matrix( c( 1:Specification1$n ) ,Specification1$n , 1)
  KXX <- CF_ExponentialFamily(X , X , Specification1$l , Specification1$p )
  invD <- CD_CalculateInverseVarDStack(KXX)
  
  Probabilities <- matrix(0 , dim(SetofSamples)[1]  + 1 ,  dim(SetofSamples)[2])
  for(kk in 1:dim(SetofSamples)[2] ){
    Probabilities[ , kk] <- CD_CalculateActualUpdatedProbabilitiesFromStack( SampleGP = as.matrix(SetofSamples[,kk]) , alpha = alpha ,KXX =  KXX , invD = invD , mu1 =  Specification1$mu , mu2 = Specification2$mu  , v1 = Specification1$v, v2 = Specification2$v )
  }
  
  return(Probabilities)
  
}
CD_CalculateSecondOrderSpecification <- function( Specification1 , Specification2 , alpha = 0.1  , numberofsamples = 100, numberinupdate = 10){
  
  SetofSamples <- cbind(CD_SampleDataFromSpecification1D(Specification = Specification1 , numberofsamples = round(alpha*numberofsamples)) , CD_SampleDataFromSpecification1D(Specification = Specification2 , numberofsamples = round((1-alpha)*numberofsamples)))
  ActualProbabilities <- CD_CalculateUpdatedProbabilitiesFromSpecifications(SetofSamples , Specification1 , Specification2 , alpha)
  
  numberofsamples <- round(alpha*numberofsamples) + round((1-alpha)*numberofsamples)
  AdjustedProbabilities <- array(0 , c(Specification1$n + 1 , 2 , numberofsamples) )
  AdjustedProbabilities[1 , 1 , ] <- alpha
  AdjustedProbabilities[1 , 2 , ] <- 0
  IndependentApproximation <- array(0 , c(Specification1$n + 1 , 2 , numberofsamples) )
  
  SOS <- list()
  
  for( kk in 1:Specification1$n ){
    
    # Calculate indepedent approximation
    for( jj in 1:numberofsamples){
      
      AdjustedBeliefsStruct <-  setNames( list( AdjustedProbabilities[kk ,1 ,  jj] , 0 ) , c('E_y_X' , 'V_y_X') )
      IndendepntApproxStruct <- CD_CalculateAdjustedUpdate(AdjustedBeliefsStruct , SetofSamples[ kk , jj] , Specification1$mu , Specification2$mu , Specification1$v , Specification2$v )
      
      if(jj ==1){
        IndependentApproximation[ kk , 1 , jj ]  <-  IndendepntApproxStruct$E_B
        IndependentApproximation[ kk , 2 , jj ]  <-  IndendepntApproxStruct$V_B
      }else{
        IndependentApproximation[ kk , 1 , jj ]  <-  IndendepntApproxStruct$E_B
        IndependentApproximation[ kk , 2 , jj ]  <-  IndendepntApproxStruct$V_B 
      }
    }
    
    # Calculate SOS   
    if(kk == 2){
      SOS[[kk]] <- CD_CreateSecondorderSpecifiction( as.matrix(ActualProbabilities[kk + 1 , ])  , cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  ]) , as.matrix(IndependentApproximation[kk , 1 ,  ] )) )
    }
    if(kk > 2){
      SOS[[kk]] <- CD_CreateSecondorderSpecifiction( as.matrix(ActualProbabilities[kk + 1 , ])  , cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  ])) , as.matrix(IndependentApproximation[kk , 1 ,  ] )) )
    }
    
    # Calculate Adjusted Probabilities
    
    for( jj in 1:numberofsamples){
      if(kk == 1){ 
        AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  IndependentApproximation[kk,1 , jj] 
        AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  IndependentApproximation[kk,2 , jj] 
      }
      if(kk == 2){ 
        
        tmpSOS <- SOS[[kk]]
        
        AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS , t(cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  jj]) , as.matrix(IndependentApproximation[kk , 1 ,  jj] ))) )
        if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
        if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
        
        AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  AdjustedProbabilityStruct[[1]]
        AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  AdjustedProbabilityStruct[[2]]
        
      }
      
      if(kk>2){
        
        tmpSOS <- SOS[[kk]]
        if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
        if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
        
        AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS ,  t(cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  jj])) , as.matrix(IndependentApproximation[kk , 1 ,  jj] ))) )
        AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  AdjustedProbabilityStruct[[1]]
        AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  AdjustedProbabilityStruct[[2]]
        
      }  
    }  
    DP_WaitBar( kk/Specification1$n )
  }
  
  return( SOS )
}
CD_CalulateAdjustedBeliefsForPrevision1D <- function( SOS , Sample , Specification1 , Specification2,  alpha = 0.1 ){
  
  n <- length(SOS)
  numberinupdate <- length(SOS[[n]]$E_Y) - 2
  
  AdjustedProbabilities <- array(0 , c(n + 1 , 2 ) )
  AdjustedProbabilities[1 , 1  ] <- alpha
  AdjustedProbabilities[1 , 2  ] <- 0
  IndependentApproximation <- array(0 , c(n + 1 , 2 ) )
  
  for( kk in 1:length(SOS) ){
    
    # Calculate independent approximation
    {AdjustedBeliefsStruct <-  setNames( list( AdjustedProbabilities[kk ,1] , 0 ) , c('E_y_X' , 'V_y_X') )
    IndendepntApproxStruct <-  CD_CalculateAdjustedUpdate(AdjustedBeliefsStruct , Sample[ kk] , Specification1$mu , Specification2$mu , Specification1$v , Specification2$v )
    IndependentApproximation[ kk , 1 ]  <-  IndendepntApproxStruct$E_B
    IndependentApproximation[ kk , 2 ]  <-  IndendepntApproxStruct$V_B}
    
    
    if( kk == 1 ){ 
      AdjustedProbabilities[kk + 1 , 1 ]  <-  IndependentApproximation[kk,1] 
      AdjustedProbabilities[kk + 1 , 2 ]  <-  IndependentApproximation[kk,2] 
    }
    if( kk == 2 ){ 
      
      tmpSOS <- SOS[[kk]]
      
      AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS , t(cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ]) , as.matrix(IndependentApproximation[kk , 1 ] ))) )
      if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
      if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
      
      AdjustedProbabilities[kk + 1 , 1 ]  <-  AdjustedProbabilityStruct[[1]]
      AdjustedProbabilities[kk + 1 , 2 ]  <-  AdjustedProbabilityStruct[[2]]
      
    }
    if( kk > 2 ){
      
      tmpSOS <- SOS[[kk]]
      
      AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS ,  t(cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1])) , as.matrix(IndependentApproximation[kk , 1 ] ))) )
      if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
      if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
      
      AdjustedProbabilities[kk + 1 , 1 ]  <-  AdjustedProbabilityStruct[[1]]
      AdjustedProbabilities[kk + 1 , 2 ]  <-  AdjustedProbabilityStruct[[2]]
      
    } 
    
  }
  
  return(AdjustedProbabilities)
  
}
CD_CalculateInverseVarDStack <- function(KXX){
  
  invD <- list()
  for(jj in 2:(dim(KXX)[1] + 1)){
    ii = jj - 1
    # Recursive conditioning
    invD[[ii]] <- solve(KXX[1:(ii-1) , 1:(ii-1) ] + 0.000000001*diag(dim(as.matrix(KXX[1:(ii-1) ,1:(ii-1) ]))[1]) )
  }
  return(invD)
}
CD_CalculateActualUpdatedProbabilitiesFromStack <- function(SampleGP , alpha , KXX , invD = 0 , mu1 , mu2 , v1 , v2 ){
  
  if(is.list(invD) == FALSE){
    invD = CD_CalculateInverseVarDStack(KXX)
  }
  
  ActualProbabilities <- matrix(0, dim(SampleGP)[1] + 1 , 1)
  ActualProbabilities[1 , ] <- alpha
  

  for(jj in 2:(dim(ActualProbabilities)[1])){
    ii = jj - 1
    if(ii == 1){
      f_1 = dnorm(SampleGP[ii] , mean = mu1 , sd = sqrt( v1 ))
      f_2 = dnorm(SampleGP[ii] , mean = mu2 , sd = sqrt( v2 ))
    }
    if(ii > 1){
      # Recursive conditioning
      tmp <- invD[[ii]]
      E_D1 <- mu1 + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1)] - mu1 )
      V_D1 <- v1*(1  - (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) ))  
      
      E_D2 <- mu2 + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1)] - mu2 )
      V_D2 <- v2*(1  - (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) ))  
      
      f_1 = dnorm( SampleGP[ii] , mean = E_D1 , sd = sqrt(V_D1) )
      f_2 = dnorm( SampleGP[ii] , mean = E_D2 , sd = sqrt(V_D2) )
      rm(tmp)
    }
    ActualProbabilities[jj , 1] <-  (ActualProbabilities[jj-1 , 1]*f_1) / ((1 - ActualProbabilities[jj-1 , 1])*f_2 + ActualProbabilities[jj-1 , 1]*f_1)
    #DP_WaitBar(ii/(dim(ActualProbabilities)[1]))
  }
  return(ActualProbabilities)
}
CD_CalculateActualUpdatedProbabilitiesFromStackTP <- function(SampleGP , alpha , KXX , invD = 0 , mu1 , mu2 , v1 , v2 , df = 5 ){
  
  if(is.list(invD) == FALSE){
    invD = CD_CalculateInverseVarDStack(KXX)
  }
  
  ActualProbabilities <- matrix(0, dim(SampleGP)[1] + 1 , 1)
  ActualProbabilities[1 , ] <- alpha
  E_D1 <- matrix(0, dim(SampleGP)[1]  , 1)
  V_D1 <- matrix(0, dim(SampleGP)[1]  , 1)
  E_D2 <- matrix(0, dim(SampleGP)[1]  , 1)
  V_D2 <- matrix(0, dim(SampleGP)[1]  , 1)
  f_1  <- matrix(0, dim(SampleGP)[1]  , 1)
  f_2  <- matrix(0, dim(SampleGP)[1]  , 1)
  aa  <- matrix(0, dim(SampleGP)[1]  , 1)
  bb  <- matrix(0, dim(SampleGP)[1]  , 1)
  
  for(jj in 2:(dim(ActualProbabilities)[1])){
    {  
    ii = jj - 1
      
    if(ii == 1){
      f_1[ii] = CD_Tdensity( x=SampleGP[ii] , mean = mu1, sd = sqrt(v1) , df =df )
      f_2[ii] = CD_Tdensity( x=SampleGP[ii] , mean = mu2, sd = sqrt(v2) , df =df )
    }
    if(ii > 1){
      # Recursive conditioning
      tmp <- invD[[ii]]
      # Varaince adjustemnt for t process
      aa[ii] <- (df  + t(SampleGP[1:(ii-1)] - mu1 )%*%((1/v1)*tmp)%*%(SampleGP[1:(ii-1)] - mu1 ) -2) /(df + dim(invD[[ii]])[1] - 2)
      bb[ii] <- (df  + t(SampleGP[1:(ii-1)] - mu2 )%*%((1/v2)*tmp)%*%(SampleGP[1:(ii-1)] - mu2 ) -2) /(df + dim(invD[[ii]])[1] - 2)
      #aa[ii] <- 1
      #bb[ii] <- 1
      
      E_D1[ii] <- mu1 + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1)] - mu1 )
      V_D1[ii] <- aa[ii]*(v1*(1  - (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) )))  
      
      E_D2[ii] <- mu2 + (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(SampleGP[1:(ii-1)] - mu2 )
      V_D2[ii] <- bb[ii]*(v2*(1  - (rev(KXX[2:(ii) , 1]) )%*%tmp%*%(rev(KXX[2:(ii) , 1]) )))  
       
      f_1[ii] = CD_Tdensity( x = SampleGP[ii] , mean = E_D1[ii], sd = sqrt( V_D1[ii] ) , df = df + dim( invD[[ii]] )[1] )
      f_2[ii] = CD_Tdensity( x = SampleGP[ii] , mean = E_D2[ii], sd = sqrt( V_D2[ii] ) , df = df + dim( invD[[ii]] )[1] )
      rm(tmp)
    }
    ActualProbabilities[jj , 1] <-  (ActualProbabilities[jj-1 , 1]*f_1[ii]) / ((1 - ActualProbabilities[jj-1 , 1])*f_2[ii] + ActualProbabilities[jj-1 , 1]*f_1[ii])
  }
     #DP_WaitBar(ii/(dim(ActualProbabilities)[1]))
  }
  return(ActualProbabilities)
}
CD_Tdensity <- function(x , mean = 0 , sd = 1 , df = 5){
  sd <- sqrt(((df )/(df)))*sd
  return( (1/sd)*(dt((x - mean)/sd , df)) )
}
CD_CalculateUpdatedProbabilitiesFromSpecificationsTP <- function(SetofSamples , Specification1 , Specification2 , alpha = 0.1 , df = 5){
  X <- matrix( c( 1:Specification1$n ) ,Specification1$n , 1)
  KXX <- CF_ExponentialFamily(X , X , Specification1$l , Specification1$p )
  invD <- CD_CalculateInverseVarDStack(KXX)
  
  Probabilities <- matrix(0 , dim(SetofSamples)[1]  + 1 ,  dim(SetofSamples)[2])
  for(kk in 1:dim(SetofSamples)[2] ){
    Probabilities[ , kk] <- CD_CalculateActualUpdatedProbabilitiesFromStackTP( SampleGP = as.matrix(SetofSamples[,kk]) , alpha = alpha ,KXX =  KXX , invD = invD , mu1 =  Specification1$mu , mu2 = Specification2$mu  , v1 = Specification1$v, v2 = Specification2$v , df = df )
  }
  
  return(Probabilities)
  
}
CD_CalulateAdjustedBeliefsForPrevision1DTP <- function( SOS , Sample , Specification1 , Specification2,  alpha = 0.1 , df = 5){
  
  n <- length(SOS)
  numberinupdate <- length(SOS[[n]]$E_Y) - 2
  
  AdjustedProbabilities <- array(0 , c(n + 1 , 2 ) )
  AdjustedProbabilities[1 , 1  ] <- alpha
  AdjustedProbabilities[1 , 2  ] <- 0
  IndependentApproximation <- array(0 , c(n + 1 , 2 ) )
  
  for( kk in 1:length(SOS) ){
    
    # Calculate independent approximation
    {AdjustedBeliefsStruct <-  setNames( list( AdjustedProbabilities[kk ,1] , 0 ) , c('E_y_X' , 'V_y_X') )
    IndendepntApproxStruct <-  CD_CalculateAdjustedUpdateTP(AdjustedBeliefsStruct , Sample[ kk] , Specification1$mu , Specification2$mu , Specification1$v , Specification2$v  , df = df)
    IndependentApproximation[ kk , 1 ]  <-  IndendepntApproxStruct$E_B
    IndependentApproximation[ kk , 2 ]  <-  IndendepntApproxStruct$V_B}
    
    
    if( kk == 1 ){ 
      AdjustedProbabilities[kk + 1 , 1 ]  <-  IndependentApproximation[kk,1] 
      AdjustedProbabilities[kk + 1 , 2 ]  <-  IndependentApproximation[kk,2] 
    }
    if( kk == 2 ){ 
      
      tmpSOS <- SOS[[kk]]
      
      AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS , t(cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ]) , as.matrix(IndependentApproximation[kk , 1 ] ))) )
      if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
      if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
      
      AdjustedProbabilities[kk + 1 , 1 ]  <-  AdjustedProbabilityStruct[[1]]
      AdjustedProbabilities[kk + 1 , 2 ]  <-  AdjustedProbabilityStruct[[2]]
      
    }
    if( kk > 2 ){
      
      tmpSOS <- SOS[[kk]]
      
      AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS ,  t(cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1])) , as.matrix(IndependentApproximation[kk , 1 ] ))) )
      if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
      if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
      
      AdjustedProbabilities[kk + 1 , 1 ]  <-  AdjustedProbabilityStruct[[1]]
      AdjustedProbabilities[kk + 1 , 2 ]  <-  AdjustedProbabilityStruct[[2]]
      
    } 
    
  }
  
  return(AdjustedProbabilities)
  
}
CD_CalculateAdjustedUpdateTP <- function(AdjustedBeliefsStruct , SampleGP , mu1 , mu2 , v1 , v2 , df = 5 ){
  SampleProbabilities <- BE_SampleLHSinab( a = (AdjustedBeliefsStruct$E_y_X - 2*sqrt(AdjustedBeliefsStruct$V_y_X)) , b = (AdjustedBeliefsStruct$E_y_X + 2*sqrt(AdjustedBeliefsStruct$V_y_X)) , numbersamples = 100)
  f_1 = CD_Tdensity(x= SampleGP , mean = mu1 , sd = sqrt( v1 ) , df )
  f_2 = CD_Tdensity(x= SampleGP , mean = mu2 , sd = sqrt( v2 )  , df)
  
  SampleProbabilities <- apply(SampleProbabilities , 1 , function(X){ (f_1*X)/ (f_1*X   +   f_2*(1-X)) } )
  
  output <- setNames( list( mean(SampleProbabilities) , var(SampleProbabilities))  , c('E_B' , 'V_B') )
  return(output)  
}
CD_CalulateAdjustedBeliefsforPrevisonFromSpecifications <- function(SOS , SetofSamples , Specification1 , Specification2,  alpha = alpha){
  
  AdjustedProbabilties <- matrix(0 , dim(SetofSamples)[1] + 1, dim(SetofSamples)[2])
  for(i in 1:dim( AdjustedProbabilties )[2]){
    AdjustedProbabilties[,i] <- CD_CalulateAdjustedBeliefsForPrevision1D( SOS , as.matrix(SetofSamples[ , i]) , Specification1 , Specification2,  alpha = alpha )[,1]
  }
  return(AdjustedProbabilties)
  
}
CD_CalulateAdjustedBeliefsforPrevisonFromSpecificationsTP <- function(SOS , SetofSamples , Specification1 , Specification2,  alpha = 0.1 , df = 5){
  
  AdjustedProbabilties <- matrix(0 , dim(SetofSamples)[1] + 1, dim(SetofSamples)[2])
  for(i in 1:dim( AdjustedProbabilties )[2]){
    AdjustedProbabilties[,i] <- CD_CalulateAdjustedBeliefsForPrevision1DTP( SOS , as.matrix(SetofSamples[ , i]) , Specification1 , Specification2,  alpha = alpha , df = df)[,1]
  }
  return(AdjustedProbabilties)
  
}
CD_CalculateUpdatedProbabilitiesFromSpecificationsMO <- function(SetofSamples , Specification1 , Specification2 , alpha = 0.1){
  X <- matrix( c( 1:Specification1$n ) ,Specification1$n , 1)
  KXX <- CF_ExponentialFamily(X , X , Specification1$l , Specification1$p )
  invD <- CD_CalculateInverseVarDStack(KXX)
  
  Probabilities <- matrix(0 , dim(SetofSamples)[1]  + 1 ,  dim(SetofSamples)[3])
  for(kk in 1:dim(SetofSamples)[3] ){
    Probabilities[ , kk] <- CD_CalculateActualUpdatedProbabilitiesMOGP( SampleGP = as.matrix(SetofSamples[,,kk]) , alpha = alpha ,KXX =  KXX , invD = invD , mu1 =  Specification1$mu , mu2 = Specification2$mu  , Sigma1 = Specification1$v, Sigma2 = Specification2$v )
  DP_WaitBar(kk/dim(SetofSamples)[3])
    }
  
  return(Probabilities)
  
}
CD_CalculateAdjustedUpdateMO <- function(AdjustedBeliefsStruct , SampleGP , Specification1 , Specification2 ){
  SampleProbabilities <- BE_SampleLHSinab( a = (AdjustedBeliefsStruct$E_y_X - 2*sqrt(AdjustedBeliefsStruct$V_y_X)) , b = (AdjustedBeliefsStruct$E_y_X + 2*sqrt(AdjustedBeliefsStruct$V_y_X)) , numbersamples = 100)
  f_1 = dmvnorm( SampleGP , mean = Specification1$mu , sigma = Specification1$v )
  f_2 = dmvnorm( SampleGP , mean = Specification2$mu , sigma = Specification2$v  )
  
  SampleProbabilities <- apply(as.matrix(SampleProbabilities) , 1 , function(X){ (f_1*X)/ (f_1*X   +   f_2*(1-X)) } )
  
  output <- setNames( list( mean(SampleProbabilities) , var(SampleProbabilities))  , c('E_B' , 'V_B') )
  return(output)  
}
CD_CalculateSecondOrderSpecificationMO <- function( Specification1 , Specification2 , alpha = 0.1  , numberofsamples = 100, numberinupdate = 10){
  
  # Sample dataset from Gaussian Process
  SetofSamples <- array(0 , c( Specification1$n , length(Specification1$mu) , round(alpha*numberofsamples) + round((1-alpha)*numberofsamples)) )
  SetofSamples[,,1:round(alpha*numberofsamples)] <- CD_SampleDataFromSpecification1DMO( Specification = Specification1 ,  numberofsamples = round(alpha*numberofsamples) )
  SetofSamples[,,(round(alpha*numberofsamples) + 1):size(SetofSamples)[3]] <- CD_SampleDataFromSpecification1DMO( Specification = Specification2 ,  numberofsamples = round((1-alpha)*numberofsamples) ) 
  
  # Calculate Actual Probabilities
  ActualProbabilities <- CD_CalculateUpdatedProbabilitiesFromSpecificationsMO(SetofSamples , Specification1 , Specification2 , alpha)
  
  numberofsamples <- round(alpha*numberofsamples) + round((1-alpha)*numberofsamples)
  AdjustedProbabilities <- array(0 , c(Specification1$n + 1 , 2 , numberofsamples) )
  AdjustedProbabilities[1 , 1 , ] <- alpha
  AdjustedProbabilities[1 , 2 , ] <- 0
  IndependentApproximation <- array(0 , c(Specification1$n + 1 , 2 , numberofsamples) )
  
  SOS <- list()
  
  for( kk in 1:Specification1$n ){
    
    # Calculate indepedent approximation
    for( jj in 1:numberofsamples){
      
      AdjustedBeliefsStruct <-  setNames( list( AdjustedProbabilities[kk ,1 ,  jj] , AdjustedProbabilities[kk ,2 ,  jj] ) , c('E_y_X' , 'V_y_X') )
      IndendepntApproxStruct <- CD_CalculateAdjustedUpdateMO(AdjustedBeliefsStruct , SampleGP = SetofSamples[ kk , , jj] , Specification1 , Specification2 )
      
      if(jj ==1){
        IndependentApproximation[ kk , 1 , jj ]  <-  IndendepntApproxStruct$E_B
        IndependentApproximation[ kk , 2 , jj ]  <-  IndendepntApproxStruct$V_B
      }else{
        IndependentApproximation[ kk , 1 , jj ]  <-  IndendepntApproxStruct$E_B
        IndependentApproximation[ kk , 2 , jj ]  <-  IndendepntApproxStruct$V_B 
      }
    }
    
    # Calculate SOS   
    if(kk == 2){
      SOS[[kk]] <- CD_CreateSecondorderSpecifiction( as.matrix(ActualProbabilities[kk + 1 , ])  , cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  ]) , as.matrix(IndependentApproximation[kk , 1 ,  ] )) )
    }
    if(kk > 2){
      SOS[[kk]] <- CD_CreateSecondorderSpecifiction( as.matrix(ActualProbabilities[kk + 1 , ])  , cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  ])) , as.matrix(IndependentApproximation[kk , 1 ,  ] )) )
    }
    
    # Calculate Adjusted Probabilities
    
    for( jj in 1:numberofsamples){
      if(kk == 1){ 
        AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  IndependentApproximation[kk,1 , jj] 
        AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  IndependentApproximation[kk,2 , jj] 
      }
      if(kk == 2){ 
        
        tmpSOS <- SOS[[kk]]
        
        AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS , t(cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  jj]) , as.matrix(IndependentApproximation[kk , 1 ,  jj] ))) )
        if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
        if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
        
        AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  AdjustedProbabilityStruct[[1]]
        AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  AdjustedProbabilityStruct[[2]]
        
      }
      
      if(kk>2){
        
        tmpSOS <- SOS[[kk]]
        if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
        if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
        
        AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS ,  t(cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  jj])) , as.matrix(IndependentApproximation[kk , 1 ,  jj] ))) )
        AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  AdjustedProbabilityStruct[[1]]
        AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  AdjustedProbabilityStruct[[2]]
        
      }  
    }  
    DP_WaitBar( kk/Specification1$n )
  }
  
  return( SOS )
}
CD_CalulateAdjustedBeliefsForPrevision1DMO <- function( SOS , Sample , Specification1 , Specification2,  alpha = 0.1 ){
  
  n <- length(SOS)
  numberinupdate <- length(SOS[[n]]$E_Y) - 2
  
  AdjustedProbabilities <- array(0 , c(n + 1 , 2 ) )
  AdjustedProbabilities[1 , 1  ] <- alpha
  AdjustedProbabilities[1 , 2  ] <- 0
  IndependentApproximation <- array(0 , c(n + 1 , 2 ) )
  
  for( kk in 1:length(SOS) ){
    
    # Calculate independent approximation
    {AdjustedBeliefsStruct <-  setNames( list( AdjustedProbabilities[kk ,1] , 0 ) , c('E_y_X' , 'V_y_X') )
    IndendepntApproxStruct <-  CD_CalculateAdjustedUpdateMO(AdjustedBeliefsStruct , (Sample[ kk , ]) , Specification1, Specification2 )
    IndependentApproximation[ kk , 1 ]  <-  IndendepntApproxStruct$E_B
    IndependentApproximation[ kk , 2 ]  <-  IndendepntApproxStruct$V_B}
    
    
    if( kk == 1 ){ 
      AdjustedProbabilities[kk + 1 , 1 ]  <-  IndependentApproximation[kk,1] 
      AdjustedProbabilities[kk + 1 , 2 ]  <-  IndependentApproximation[kk,2] 
    }
    if( kk == 2 ){ 
      
      tmpSOS <- SOS[[kk]]
      
      AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS , t(cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ]) , as.matrix(IndependentApproximation[kk , 1 ] ))) )
      if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
      if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
      
      AdjustedProbabilities[kk + 1 , 1 ]  <-  AdjustedProbabilityStruct[[1]]
      AdjustedProbabilities[kk + 1 , 2 ]  <-  AdjustedProbabilityStruct[[2]]
      
    }
    if( kk > 2 ){
      
      tmpSOS <- SOS[[kk]]
      
      AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS ,  t(cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1])) , as.matrix(IndependentApproximation[kk , 1 ] ))) )
      if(AdjustedProbabilityStruct[[1]] < 0){AdjustedProbabilityStruct[[1]] = 0}
      if(AdjustedProbabilityStruct[[1]] > 1 ){AdjustedProbabilityStruct[[1]] = 1}
      
      AdjustedProbabilities[kk + 1 , 1 ]  <-  AdjustedProbabilityStruct[[1]]
      AdjustedProbabilities[kk + 1 , 2 ]  <-  AdjustedProbabilityStruct[[2]]
      
    } 
    
  }
  
  return(AdjustedProbabilities)
  
}