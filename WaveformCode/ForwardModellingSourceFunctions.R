FM_PlotPWaveAnalysisSingle <- function( ECG  , Beats  , QSwidth  = 0){
  
  regiontocheck <- max(1,min(dim(Beats)[1] - 500 , 1000)):min(1500 , dim(Beats)[1])
  lengthcheck <- 0
  while( lengthcheck < 250 ){
    if(regiontocheck[length(regiontocheck)] > dim(Beats)[1]){break}
    ECGBeats <- AFD_ExtractAllSQ(ECG = ECG , RPeaks = Beats[regiontocheck,] , QSwidth  = QSwidth)
    lengthcheck <- length(ECGBeats$numvalues)
    if(lengthcheck < 250 ){regiontocheck <- regiontocheck + 500}
  }
  
  
  options(expressions = 10000)
  output = ggplot( )
  if(length(ECGBeats$numvalues) > 100  ){
    for(i in 1:(dim(ECGBeats$Date)[1] -50)){
      if(ECGBeats$numvalues[i ] > as.numeric(quantile( ECGBeats$numvalues  , 0.95)) ){next}
      if(ECGBeats$numvalues[i ] < as.numeric(quantile( ECGBeats$numvalues  , 0.5)) ){next}
      tmp = 1:(min(dim(ECGBeats$Date)[2] , ECGBeats$numvalues[i ]) -1)
      output <-   output + geom_line(data = data.frame(time = (1:length(ECGBeats$Date[i,tmp]))/length(ECGBeats$Date[i,tmp]),
                                                       V =  ECGBeats$Value[i,tmp] - mean(ECGBeats$Value[i,tmp])) ,
                                     aes(time , V) , color =  rgb(0 , 0 , 1 , alpha = 0.05) ) 
    }  
  }
  return(output+ ylim(-30,30) + xlim(0.5, 1))
}
FM_GMMDensity <- function(x , Pi =1 , mu =0 , sigma =1){
  f <- matrix(0 , length(x) , 1)
  for(i in 1:length(Pi)){
  f <- f + Pi[i]*dnorm(x = x , mean = mu[i] , sd = sigma[i])  
  }
  return(f)
}
FM_dGMMDensity <- function(x , Pi =1 , mu =0 , sigma =1){
  f <- matrix(0 , length(x) , 1)
  for(i in 1:length(Pi)){
    f <- f + Pi[i]*ddnorm(x = x , mean = mu[i] , sd = sigma[i])  
  }
  return(f)
}
FM_HistoryMatchRRDensity <- function( PriorNonImplausibleSet , x , f_x , MD , RRtimes  , imthreshold = 0.6){
  
  # Build kdeestimate for history matching  
  m <- rollmedian( RRtimes , k = 21 , na.pad = T )
  mm <- median( RRtimes )
  
  RRtimes <- RRtimes - m + mm + rnorm(length(RRtimes) , 0 , 0.0025 )
  RRtimes <- RRtimes[!is.na(RRtimes)]
  
  RPeakKDEEstmate <- kde( RRtimes )
  y <- predict(RPeakKDEEstmate , x = x)
  
  NonZeroLogical <- (y/max(y)) > 0.001
  NonZeroLogical[1:which(NonZeroLogical)[length(which(NonZeroLogical))]] <- T
  
  
  # Im <- colMeans( apply( f_x[ , NonZeroLogical] , 1 , function(X){abs(X - y[ NonZeroLogical])} )/ t(MD[ , NonZeroLogical]) )
    Im <- colMeans( apply( f_x[ , NonZeroLogical] , 1 , function(X){abs(X - y[ NonZeroLogical])} )/ t(MD[ , NonZeroLogical]) )
    maxIm <- apply( apply( f_x[ , NonZeroLogical] , 1 , function(X){abs(X - y[ NonZeroLogical])} )/ t(MD[ , NonZeroLogical]) , 2, max )

  if(sum(Im <= imthreshold) > 1){
    return( setNames(list( PriorNonImplausibleSet[Im < imthreshold,] , Im[Im < imthreshold] , f_x[Im < imthreshold ,  ] , y ) , c('NonImplausibleSets' , 'Implausability' , 'f_x' , 'y' )  ) )
  }else{
    return( setNames(list( PriorNonImplausibleSet[which.min(Im),] , min(Im) , f_x[which.min(Im) ,  ] , y ) , c('MinImplausiblePoint' , 'Implausability' , 'f_x' , 'y' ) ) )
  }
}
FM_HistoryMatchRRCulmativeDensity <- function( PriorNonImplausibleSet , x , F_x ,f_x ,  MD , RRtimes  , Corr_sdhat , imthreshold = 3  ){
  set.seed(1)
  # Build kdeestimate for history matching  
  m <- rollmedian( RRtimes , k = 21 , na.pad = T )
  mm <- median( RRtimes )
  
  RRtimes <- RRtimes - m + mm + rnorm(length(RRtimes) , 0 , 0.0025 )
  RRtimes <- RRtimes[!is.na(RRtimes)]
  
  y <- FM_CalculateCDFS( RRtimes = RRtimes , xx = x )

  NonZeroLogical <- c( ((y > 0.03)*(y < 0.97))==1 ) ==1
  
  if(sum(NonZeroLogical) <= 1){
    NonZeroLogical <- c( ((y > 0.01)*(y < 0.99))==1 ) ==1
  }
  
  
  Im2 <- colMeans( apply( as.matrix(F_x[ , NonZeroLogical]) , 1 , function(X){abs(X - y[ NonZeroLogical])} )/ t(MD[ , NonZeroLogical] + 0.0045)  , na.rm = T)
  Im <- apply( apply( as.matrix(F_x[ , NonZeroLogical]) , 1 , function(X){abs(X - y[ NonZeroLogical])} ) / t(MD[ , NonZeroLogical] + 0.0045) , 2 , function(X){max(X , na.rm = T)} )
  RPeakKDEEstmate <- kde( RRtimes )
  z <- predict(RPeakKDEEstmate , x = x - 0.005)
  
  if(sum(Im <= imthreshold) > 1){
  
  LogicalVector <- (Im < imthreshold)
  Im3 <- matrix(0 , sum(LogicalVector) , 1)
  
  for(i in 1:dim(Im3)[1]){
    #covMattmp <- DP_AddNugget( ((MD[which(LogicalVector)[i] , NonZeroLogical] + 0.0045) %*% t(MD[which(LogicalVector)[i] , NonZeroLogical]+ 0.0045)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[which(LogicalVector)[i] , NonZeroLogical]+ 0.0045) )
    Im3[i, ] <- (abs((z[NonZeroLogical] - f_x[which(LogicalVector)[i] , NonZeroLogical])/f_x[which(LogicalVector)[i] , NonZeroLogical]  ))[which.max(z[NonZeroLogical]) ]
    #Im3[i, ] <- mahalanobis(x = y[NonZeroLogical] , center = F_x[which(LogicalVector)[i] , NonZeroLogical] , cov = covMattmp  ) 
  }
  }
  if(sum(Im <= imthreshold) <= 1){
    Im3<- (abs( (z[NonZeroLogical] - f_x[which.min(Im) , NonZeroLogical])/f_x[which.min(Im) , NonZeroLogical] ) )[which.max(z[NonZeroLogical])]
    #covMattmp <- DP_AddNugget( ((MD[which.min(Im) , NonZeroLogical] + 0.0045) %*% t(MD[which.min(Im) , NonZeroLogical]+ 0.0045)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[which.min(Im)  , NonZeroLogical]+ 0.0045) )
    #Im3 <- mahalanobis(x = y[NonZeroLogical] , center = F_x[which.min(Im) , NonZeroLogical] , cov = covMattmp  ) 
  }
  
  if(sum(Im <= imthreshold) > 1){
    return( setNames(list( PriorNonImplausibleSet[Im < imthreshold,] , cbind(Im[Im < imthreshold] , Im2[Im < imthreshold] , Im3) , F_x[Im < imthreshold ,  ] , y  ) , c('NonImplausibleSets' , 'Implausability' , 'f_x' , 'y'  )  ) )
  }else{
    return( setNames(list( PriorNonImplausibleSet[which.min(Im),] , cbind(min(Im ) , Im2[which.min(Im)] , Im3) , F_x[which.min(Im) ,  ] , y  ) , c('MinImplausiblePoint' , 'Implausability' , 'f_x' , 'y'  ) ) )
  }
}
FM_SampleRealisationsSet <- function( PriorNonImplausibleSet , N , x ){
  output <- matrix( 0 , dim(PriorNonImplausibleSet)[1] , length(x) )
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    output[i , ] <-  FM_CalculateCDFS( RRtimes = FM_SampleGMM( X = PriorNonImplausibleSet[i,] , N = N )  , xx = x) 
    #output[i , ] <-  FM_CalculateCDFS( RRtimes = FM_EvaluateDenistyEstimate( x = x , X = PriorNonImplausibleSet[i,] )  , xx = x) 
  }
  return(output)
}
FM_CalulateImForGroundTruth <- function(x , F_x , f_x , PriorNonImplausibleSet , MD , Corr_sdhat , N = 231  ){
  
  ImplausabilityMatrix <- matrix(0 , dim(PriorNonImplausibleSet)[1] , 3)
  
  for( ii in 1:dim(PriorNonImplausibleSet)[1] ){
    
    RRtimes <- FM_SampleGMM(X = PriorNonImplausibleSet[ii , ] ,  N)
    y <- FM_CalculateCDFS( RRtimes = RRtimes , xx = x )
    NonZeroLogical <- c( ((y > 0.03)*(y < 0.97))==1 ) ==1
    
    if(sum(NonZeroLogical) <= 2){
      NonZeroLogical <- c( ((y > 0.01)*(y < 0.99))==1 ) ==1
    }
    if(sum(NonZeroLogical) <= 2){
      NonZeroLogical <- c( ((y > 0.001)*(y < 0.999))==1 ) ==1
    }
    
    
    #mean( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical] + 0.003) , na.rm = T)
    #max( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical]+ 0.003) , na.rm = T)
    #covMattmp <- DP_AddNugget( ((MD[ii , NonZeroLogical]+ 0.003) %*% t(MD[ii , NonZeroLogical]+ 0.003)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[ii , NonZeroLogical]+ 0.00000001) )
    #mahalanobis(x = y[NonZeroLogical] , center = F_x[ii , NonZeroLogical] , cov = covMattmp  ) 
    
    ImplausabilityMatrix[ii , 1] <- mean( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical] + 0.0045) , na.rm = T)
    ImplausabilityMatrix[ii , 2] <- max( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical]+ 0.0045) , na.rm = T)
    #covMattmp <- DP_AddNugget( ((MD[ii , NonZeroLogical] + 0.0045) %*% t(MD[ii , NonZeroLogical]+ 0.0045)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[ii , NonZeroLogical]+ 0.0045) )
    #ImplausabilityMatrix[ii , 3] <- mahalanobis(x = y[NonZeroLogical] , center = F_x[ii , NonZeroLogical] , cov = covMattmp  ) 
    
    RPeakKDEEstmate <- kde( RRtimes )
    z <- predict(RPeakKDEEstmate , x = x - 0.005)
    
    
    ImplausabilityMatrix[ii , 3] <- (abs( (z[NonZeroLogical] - f_x[ii , NonZeroLogical])/f_x[ii , NonZeroLogical] ) )[which.max(z[NonZeroLogical])]
    
  }
  return(ImplausabilityMatrix)
  
}
FM_dlaplace <- function( x , mu = 0 , sigma = 1 , alpha = 0 ){
  return(alpha*dnorm(x , mu , sigma) + (1-alpha)*dlaplace(x , mu , sigma/sqrt(2)) )
}
FM_rlaplace <- function( n , mu = 0 , sigma = 1  , alpha ){
  if( alpha == 0 ){
  return(rlaplace( n , mu , sigma/sqrt(2) ))
  }
  if(alpha == 1){
  return(rnorm(n , mu , sigma))
  }
  if(alpha != 1 || alpha != 0){
    SampleMultinomial <- t(rmultinom(n  , size = 1 , c( alpha , 1-alpha ) ))
    output <- matrix(0 , dim(SampleMultinomial)[1] , 1)
    output[SampleMultinomial[,1] == 1, ] <- rnorm(n = sum(SampleMultinomial[,1])  , mean  = mu , sd = sigma )
    output[SampleMultinomial[,2] == 1, ] <- rlaplace( n = sum(SampleMultinomial[,2]) , m =  mu , s  = sigma/sqrt(2))
  }
  return(output)
}
FM_EvaluateDenistyEstimate <- function(x , X){
 a <-  FM_GMMDensity(x = x , Pi = X[c(1,3) ] , mu = X[c(4,6) ] , sigma = X[c(7,9) ] ) 
 b <- X[2]*FM_dlaplace(x =x , mu = X[5] , sigma = X[8] , alpha = X[10] ) 
 return(a + b)
}
FM_CalculateCDFS  <- function(RRtimes , xx = seq(0.25 , 2 , 0.01)){
  output <- matrix(0 , length(xx) , 1)
  for(i in 1:length(xx)){
    output[i] <- sum(RRtimes <= xx[i])/length(RRtimes)   
  }
  return(output)
}
FM_SampleGMM <- function( X , N = 250 ){

  # sample multinomial distribution for mixtures model
  SampleMultinomial <- t(rmultinom(N  , size = 1 , X[1:3]))
  
  
  output <- matrix(0 , N , 1)
  for(i in 1:dim(SampleMultinomial)[2]){
    if(sum(SampleMultinomial[,i] == 1) ==0){next} # If not sample for componet i: skip
    if( X[10] == 1 ){
    # If GMM 
    output[SampleMultinomial[,i] == 1 , ] <- rnorm(n = sum(SampleMultinomial[,i]) , mean  = X[3+i] , sd = X[6 + i])
    }else{
      # If mixture of GMM and laplace 
      output[SampleMultinomial[,i] == 1 , ] <- FM_rlaplace(n = sum(SampleMultinomial[,i] ) , mu  = X[3+i] , sigma = X[6 + i] , alpha = X[10])
    }
    
  }
  return(output)
}
FM_SampleRealisations <- function(X = c(1 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1) , x ,  N = 250 , NN = 1000){
  output <- matrix(0 , NN , length(x))
  for(i in 1:NN){
    output[i , ] <-  FM_CalculateCDFS( RRtimes = FM_SampleGMM(X , N = N)  , xx = x) 
  }
  return(output)
}
FM_CreateDefaultFeatureTable <- function( ReHmOutput , ReIreHmOutput ){
  
  d <- matrix(' ' , 6 , 5)
  
  colnames(d) <- c('Feature' , 'Note' , 'Non-implausible' , 'Expectation' , 'Variance' )
  
  d[1 , 1] <- 'Heart Rhythum (bpm)'
  d[2 , 1] <- 'P-wave Amplitude (mV)'
  d[3 , 1] <- 'P-wave Durations (ms)'
  d[4 , 1] <- 'P-wave Interval (ms)'

  d[1 , 2] <- 'Regular'
  
  return(d)
}
FMPWaveHM_EmulateTQSegment <- function( QS_Struct , EmulatorParameters = PWaveHM_CreateDefaultEmulationclass() , Xstar  , upper = 0.9 , lower = 0.1){
  
  tmp = 1:(min(dim(QS_Struct$Date)[2] , QS_Struct$numvalues[1 ]) -1)
  EmulatorParameters$X <- (1:length(QS_Struct$Date[1,tmp]))/length(QS_Struct$Date[1,tmp])
  EmulatorParameters$Y <- QS_Struct$Value[1,tmp] - mean(QS_Struct$Value[1,tmp])
  
  EmulatedQS <- matrix( 0 , dim(QS_Struct$Date)[1] , length(Xstar) )
  #EmulatedQS <- matrix( 0 , 100 , length(Xstar) )
  counter <- 1
  print( 'Emulating TQ Segments.' )
  for(i in 1:dim(QS_Struct$Date)[1]){
    #if(QS_Struct$numvalues[i ] > as.numeric( quantile( QS_Struct$numvalues  , upper)) ){next}
    #if(QS_Struct$numvalues[i ] < as.numeric( quantile( QS_Struct$numvalues  , lower)) ){next}
    tmp = 1:(min(dim(QS_Struct$Date)[2] , QS_Struct$numvalues[i ]) -1)
    EmulatorParameters$X <- (1:length(QS_Struct$Date[i,tmp]))/length(QS_Struct$Date[i,tmp])
    EmulatorParameters$Y <- QS_Struct$Value[i,tmp] - mean(QS_Struct$Value[i,tmp])
    if(length(  EmulatorParameters$X ) < 10){next}
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(Xstar , EmulatorParameters  , meanonly = 1)
    EmulatedQS[counter,] <- EmulatorOutput$E_D_fX
    counter <- counter +1
    #if(counter > 100){break}
  }  
  print('TQ Segments Emulated.' )
  #EmulatedQS <- EmulatedQS[-(counter:dim(EmulatedQS)[1] ), ]
  return(EmulatedQS)
}
FMPWaveHM_BLUBetas <- function( H , z  , E_Beta , V_Beta  , V_me , alpha = 0 ){
  
  z <- as.matrix(z)
  
  if(alpha == 0){
    if(length(V_me) == 1){
      alpha <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + V_me*diag(dim(H)[1])  )}else{
        alpha <- t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + diag(V_me)  )  
      }
  }
  E_z_Beta <- E_Beta + alpha%*%(z - H%*%E_Beta) 
  V_z_Beta <- V_Beta - alpha%*%H*V_Beta
  
  return( setNames( list( E_z_Beta ,  V_z_Beta ) , c('E_z_Beta' , 'V_z_Beta') )  )
}
