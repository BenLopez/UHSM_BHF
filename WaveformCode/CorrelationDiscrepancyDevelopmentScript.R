# Script to analyse the affect of discrepancy on the calulations

{mu1 <- sum(LocalDistributionStruct[[2]]$parameters$mean[2 , ]*LocalDistributionStruct[[2]]$parameters$pro)
v1 <-  sum((LocalDistributionStruct[[2]]$parameters$mean[2 , ]^2 + LocalDistributionStruct[[2]]$parameters$variance$sigma[2 , 2 , ])*LocalDistributionStruct[[2]]$parameters$pro) - mu1^2

mu2 <- mu1 + 0.8
v2 <-  1.1*v1
#mu2 <- sum(LocalDistributionStruct[[3]]$parameters$mean[2 , ]*LocalDistributionStruct[[3]]$parameters$pro)  
#v2 <-  sum((LocalDistributionStruct[[3]]$parameters$mean[2 , ]^2 + LocalDistributionStruct[[3]]$parameters$variance$sigma[2 , 2 , ])*LocalDistributionStruct[[3]]$parameters$pro) - mu2^2

X <- matrix(c(1:1000) ,1000 , 1)
l <- 100
p <- 1
KXX <- CF_ExponentialFamily(X , X , l , p)

SampleGP1 <- mu1 + BE_SampleGP( v1*KXX )
SampleGP2 <- mu2 + BE_SampleGP( v2*KXX )
  

plot(SampleGP1 , type='l' , col = rgb(1,0,0,alpha = 0.5))
lines(SampleGP2 , type='l' , col = rgb(0,0,1,alpha = 0.5))

alpha <- 0.1

# Individual probabilities

f_i = matrix(0, dim(SampleGP1)[1] , 2)
f_i[ , 1] <- dnorm( SampleGP1 , mean = mu1  , sd = sqrt(v1))
f_i[ , 2] <- dnorm( SampleGP1 , mean = mu2  , sd = sqrt(v2))

IndividualProbabilities = (alpha*f_i[,1]) / ((1 - alpha)*f_i[,2] + alpha*f_i[,1])
plot(IndividualProbabilities , type = 'l' ,  col = rgb(1,0,0 , alpha = 0.5) , ylim = c(0,1))
}
# Seqentially Updated
JointProbabilties <- matrix(0, dim(SampleGP1)[1] + 1 , 1)
JointProbabilties[1, 1] = alpha 

for(ii in 2:(dim(JointProbabilties)[1] - 1)){
  JointProbabilties[ii , 1] <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])
}

lines(JointProbabilties , type = 'l' , col = rgb(0,1,0 , alpha = 0.5))


# Actual probabilities

ActualProbabilities <- matrix(0, dim(SampleGP1)[1] + 1 , 2)
ActualProbabilities[1 , ] <- alpha

for(jj in 2:(dim(ActualProbabilities)[1])){
  ii = jj - 1
  if(ii > 1){
  invD <- solve(KXX[1:(ii-1) , 1:(ii-1) ] + 0.000000001*diag(dim(as.matrix(KXX[1:(ii-1) ,1:(ii-1) ]))[1]) )
  E_D1 <- mu1 + (rev(KXX[2:(ii) , 1]) )%*%invD%*%(SampleGP1[1:(ii-1)] - mu1 )
  V_D1 <- v1*(1  - (rev(KXX[2:(ii) , 1]) )%*%invD%*%(rev(KXX[2:(ii) , 1]) ))  
  
  E_D2 <- mu2 + (rev(KXX[2:(ii) , 1]) )%*%invD%*%(SampleGP1[1:(ii-1)] - mu2 )
  V_D2 <- v2*(1  - (rev(KXX[2:(ii) , 1]) )%*%invD%*%(rev(KXX[2:(ii) , 1]) ))  

  f_1 = dnorm( SampleGP1[ii] , mean = E_D1 , sd = sqrt(V_D1) )
  f_2 = dnorm( SampleGP1[ii] , mean = E_D2 , sd = sqrt(V_D2) )
  f_3 = dnorm( SampleGP2[ii] , mean = E_D1 , sd = sqrt(V_D1) )
  f_4 = dnorm( SampleGP2[ii] , mean = E_D2 , sd = sqrt(V_D2) )
  }
  if(ii == 1){
  f_1 = dnorm(SampleGP1[ii] , mean = mu1 , sd = sqrt(v1))
  f_2 = dnorm(SampleGP1[ii] , mean = mu2 , sd = sqrt(v1))
  f_3 = dnorm(SampleGP2[ii] , mean = mu1 , sd = sqrt(v2))
  f_4 = dnorm(SampleGP2[ii] , mean = mu2 , sd = sqrt(v2))
  }
  ActualProbabilities[jj , 1] <-  (ActualProbabilities[jj-1 , 1]*f_1) / ((1 - ActualProbabilities[jj-1 , 1])*f_2 + ActualProbabilities[jj-1 , 1]*f_1)
  ActualProbabilities[jj , 2] <-  (ActualProbabilities[jj-1 , 2]*f_3) / ((1 - ActualProbabilities[jj-1 , 2])*f_4 + ActualProbabilities[jj-1 , 2]*f_3)
  
  DP_WaitBar(ii/(dim(ActualProbabilities)[1]))
}

ActualProbabilities[dim(ActualProbabilities)[1] , ] <-  0

plot( ActualProbabilities[ , 1] , type = 'l' , col = rgb(0,0,1 , alpha = 0.5) , ylim = c(0,1))
lines( JointProbabilties[ , 1 ], col = rgb(0,1,0 , alpha = 0.5))
lines( IndividualProbabilities, col = rgb(1,0,0 , alpha = 0.5))


beta = c(seq(0,1,0.001))
meansquarederror <- matrix(0 , length(beta ) , 1)
for(jj in 1:length(beta)){
for(ii in 2:(dim(JointProbabilties)[1] - 1)){
  if(ii == 2){ JointProbabilties[ii , 1] <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])}
  if(ii > 2){
  tmp <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])
  JointProbabilties[ii , 1] <-  beta[jj]*JointProbabilties[ii-1 , 1] + (1-beta[jj])*tmp
}
}
  meansquarederror[jj,1] <- mean((ActualProbabilities[ , 1] - JointProbabilties[ , 1 ])^2)
}

beta = beta[ which.min(meansquarederror) ]
for(ii in 2:(dim(JointProbabilties)[1] - 1)){
  if(ii == 2){ JointProbabilties[ii , 1] <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])}
  if(ii > 2){
    tmp <-  (JointProbabilties[ii-1 , 1]*f_i[ii-1,1]) / ((1 - JointProbabilties[ii-1 , 1])*f_i[ii-1,2] + JointProbabilties[ii-1 , 1]*f_i[ii-1,1])
    JointProbabilties[ii , 1] <-  beta*JointProbabilties[ii-1 , 1] + (1-beta)*tmp
  }
}

plot( ActualProbabilities[ , 1] - JointProbabilties[ , 1 ] , type = 'l' , col = rgb(0,0,1 , alpha = 0.5) , ylim = c(-1,1))
lines( ActualProbabilities[ , 1] , type = 'l' , col = rgb(0,1,0 , alpha = 0.5) )
lines( JointProbabilties[ , 1 ], col = rgb(1,0,0 , alpha = 0.5)) 
title(paste0('Mean = ' , mean(ActualProbabilities[ , 1] - JointProbabilties[ , 1 ]) , ' Var = '  , var(ActualProbabilities[ , 1] - JointProbabilties[ , 1 ])))
lines( JointProbabilties[ , 1 ] + 2*sqrt(var(ActualProbabilities[ , 1] - JointProbabilties[ , 1 ])) , col = rgb(0,1,1 , alpha = 0.5) ) 
lines( JointProbabilties[ , 1 ] - 2*sqrt(var(ActualProbabilities[ , 1] - JointProbabilties[ , 1 ])) , col = rgb(0,1,1 , alpha = 0.5) ) 



numberofsamples <- 1000
ActualProbabilities <- array( 0 , c(numberofsamples , 2) )
IndEstProbabilities <- array( 0 , c(numberofsamples , 2) )
n = 2


for(kk in 1:numberofsamples){
  
  sampleunif <- runif(1)
  if(sampleunif < alpha){SampleGP <- mu1 + BE_SampleGP( v1*KXX[1:2 , 1:2] )}
  if(sampleunif > alpha){SampleGP <- mu2 + BE_SampleGP( v2*KXX[1:2 , 1:2] )}
  
  ActualProbabilities[kk , ] <- CD_CalculateActualUpdatedProbabilities(as.matrix(SampleGP[1:2]) , alpha , KXX[1:2 , 1:2] , mu1 , mu2 , v1 , v2 )[2:3]
  IndEstProbabilities[kk , ] <- CD_CalulateIndenpendentEstimatedProbabilities(as.matrix(SampleGP[1:2]) , alpha , mu1 , mu2 , v1 , v2 )[2:3]
 
}

SOS <- CD_CreateSecondorderSpecifiction(X = ActualProbabilities[ , 2] , Y=IndEstProbabilities)


Adjustedbeliefs <- array( 0 , c(numberofsamples , 2)  )

for(kk in 1:numberofsamples){
  
Adjustedbeliefs[kk,1] <- CD_BayesLinearAdjustmentSingleObs(SOS , IndEstProbabilities[kk,])[[1]]
Adjustedbeliefs[kk,2] <- CD_BayesLinearAdjustmentSingleObs(SOS , IndEstProbabilities[kk,])[[2]]

}


n <- 1000
numberofsamples <- 1000
SetofSamples <- matrix(0 , n , numberofsamples)
ActualProbabilities <- matrix(0 , n + 1 , numberofsamples)
SampleClass <-  matrix(0 ,  numberofsamples , 1)

for( kk in 1:numberofsamples ){
  
  sampleunif <- runif(1)
  if( sampleunif < alpha ){
    SampleGP <- mu1 + BE_SampleGP( v1*KXX[1:n , 1:n] )
    SampleClass[kk , 1] <- 1 }
  if(sampleunif > alpha){SampleGP <- mu2 + BE_SampleGP( v2*KXX[1:n , 1:n] )}
  SetofSamples[ , kk] <- SampleGP

  ActualProbabilities[ , kk] <- CD_CalculateActualUpdatedProbabilities( SampleGP , alpha , KXX , mu1 , mu2 , v1 , v2 )
  DP_WaitBar(kk/numberofsamples)
}

AdjustedProbabilities <- array(0 , c(n + 1 , 2 , numberofsamples) )
AdjustedProbabilities[1 , 1 , ] <- alpha
AdjustedProbabilities[1 , 2 , ] <- 0

IndependentApproximation <- array(0 , c(n + 1 , 2 , numberofsamples) )

SOS <- list()
numberinupdate = 100

for( kk in 1:n ){

# Calculate indepedent approximation

for( jj in 1:numberofsamples){
    
    AdjustedBeliefsStruct <-  setNames( list( AdjustedProbabilities[kk ,1 ,  jj] , 0 ) , c('E_y_X' , 'V_y_X') )
    IndendepntApproxStruct <- CD_CalculateAdjustedUpdate(AdjustedBeliefsStruct , SetofSamples[ kk , jj] , mu1 , mu2 , v1 , v2 )
  
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
    SOS[[kk]] <-CD_CreateSecondorderSpecifiction(as.matrix(ActualProbabilities[kk + 1 , ])  , cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  ]) , as.matrix(IndependentApproximation[kk , 1 ,  ] )) )
  }
  if(kk > 2){
  SOS[[kk]] <- CD_CreateSecondorderSpecifiction(as.matrix(ActualProbabilities[kk + 1 , ])  , cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  ])) , as.matrix(IndependentApproximation[kk , 1 ,  ] )) )
  }
  
  #SOS[[kk]]$E_X <- as.matrix(0.1)
  #SOS[[kk]]$E_Y <- as.matrix(c(0.1,0.1))
  
# Calculate Adjusted Probabilities
  
for( jj in 1:numberofsamples){
  if(kk == 1){ 
    AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  IndependentApproximation[kk,1 , jj] 
    AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  IndependentApproximation[kk,2 , jj] 
  }
  if(kk == 2){ 
    tmpSOS <- SOS[[kk]]
    #tmpSOS$E_X <- AdjustedProbabilities[ kk , 1 ,  jj ]
    #tmpSOS$V_X <- AdjustedProbabilities[ kk , 2 ,  jj ]
    
    AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS , t(cbind(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  jj]) , as.matrix(IndependentApproximation[kk , 1 ,  jj] ))) )
    AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  AdjustedProbabilityStruct[[1]]
    AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  AdjustedProbabilityStruct[[2]]
  }
  if(kk>2){
    tmpSOS <- SOS[[kk]]
    #tmpSOS$E_X <- AdjustedProbabilities[ kk , 1 ,  jj ]
    #tmpSOS$V_X <- AdjustedProbabilities[ kk , 2 ,  jj ]
      
    AdjustedProbabilityStruct <-   CD_BayesLinearAdjustmentSingleObs( tmpSOS ,  t(cbind(t(as.matrix(AdjustedProbabilities[max(2 , (kk- numberinupdate)):kk , 1 ,  jj])) , as.matrix(IndependentApproximation[kk , 1 ,  jj] ))) )
    AdjustedProbabilities[kk + 1 , 1 ,  jj]  <-  AdjustedProbabilityStruct[[1]]
    AdjustedProbabilities[kk + 1 , 2 ,  jj]  <-  AdjustedProbabilityStruct[[2]]
    }  
  }  
  DP_WaitBar(kk/n)
}


{par(mfrow = c(1,1))
indiex <- 49
plot((AdjustedProbabilities[  , 1,indiex] ) , type ='l' , col = rgb(0,0,1 , alpha = 0.5) , ylim = c(0,1))
lines((AdjustedProbabilities[  , 1,indiex] ) - 3*sqrt(AdjustedProbabilities[  , 2 , indiex] ) , type ='l' , col = rgb(0,1,1 , alpha = 0.5))
lines((AdjustedProbabilities[  , 1,indiex] ) + 3*sqrt(AdjustedProbabilities[  , 2 , indiex] ) , type ='l' , col = rgb(0,1,1 , alpha = 0.5))

lines((ActualProbabilities[  ,indiex] ) , type ='l' , col = rgb(1,0,0, alpha = 0.5))}


ExpectationX <- matrix(0 , n - (numberinupdate+2)  , 1)
ExpectationY <- matrix(0 , n - (numberinupdate+2)  , numberinupdate+2)
Weights <- matrix(0 , n - (numberinupdate+2)  , numberinupdate+2)
CorrelationM <- matrix(1 , n - (numberinupdate+2)  , numberinupdate+2)
VarianceX <- matrix(0 , n - (numberinupdate+2)  , 1)
VarianceY <- array(0 , c( n - (numberinupdate+2)  , numberinupdate+2 , numberinupdate+2))
CovarianceXY <- matrix(0 , n - (numberinupdate+2) , numberinupdate+2)

for(i in (numberinupdate+2):n){
  ExpectationX[i- (numberinupdate+2)] <- SOS[[i]]$E_X
  ExpectationY[i- (numberinupdate+2) , ] <- SOS[[i]]$E_Y
  if(i > numberinupdate + 2){
  Weights[i- (numberinupdate+2),] <- t(SOS[[i]]$Cov_XY)%*%solve(SOS[[i]]$V_Y)
  CorrelationM[ i- (numberinupdate+2) ] <- SOS[[i]]$V_Y[1 , 2] / sqrt(prod(diag(SOS[[i]]$V_Y)))
  VarianceX[ i- (numberinupdate+2) ]  <-  SOS[[i]]$V_X
  VarianceY[i- (numberinupdate+2),,]  <-  SOS[[i]]$V_Y
  CovarianceXY[i- (numberinupdate+2),] <- SOS[[i]]$Cov_XY
  }
}


{
x11(20,14)
par(mfrow = c(4,1))
plot(ExpectationX , type ='l', col = rgb(0,0,1 , alpha = 0.5) , ylim = c(0.08,0.11)  , xlab = 'time' )
lines(ExpectationY[,1], col = rgb(1,0,0 , alpha = 0.5))
lines(ExpectationY[,2], col = rgb(0,1,0 , alpha = 0.5))
title('SOS Expectations')
plot( Weights[ , 1] ,  type ='l', col = rgb(0,0,1 , alpha = 0.5)   , xlab = 'time' )
lines( Weights[,2]  ,  col = rgb(1,0,0 , alpha = 0.5))
lines( Weights[,3]  ,  col = rgb(0,1,0 , alpha = 0.5))
lines( Weights[,25]  ,  col = rgb(1,1,0 , alpha = 0.5))


title('SOS Weights')
plot(CorrelationM  , col = rgb(0,0,1 , alpha = 0.5) , ylim = c(0,1)  , xlab = 'time' )
title('SOS Correlation')
plot(AdjustedProbabilities[  , 2,1] , col = rgb(0,0,1 , alpha = 0.5)  , xlab = 'time' )
title('Adjusted Variance')
}


