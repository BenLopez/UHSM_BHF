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

CD_CalulateIndenpendentEstimatedProbabilities <- function(SampleGP , alpha ,mu1 , mu2 , v1,v2 ){
  
  f_i = CD_CalculateIndividualDensities(SampleGP , mu1 , mu2 , v1 , v2 )
  
  JointProbabilties <- matrix(0, dim(SampleGP)[1] + 1 , 1)
  JointProbabilties[1, 1] = alpha 
  
  for(ii in 2:dim(JointProbabilties)[1]){
    JointProbabilties[ii , 1] <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])
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