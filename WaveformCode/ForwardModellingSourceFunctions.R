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

FM_QuadraticEstimate <- function(Y){
  H = cbind(matrix(1 , length(Y) , 1) , matrix(1:length(Y) , length(Y) , 1) , matrix(c(1:length(Y))^2 , length(Y) , 1)  )
  Beta = solve(t(H)%*%H)%*%t(H)%*%Y
  return(H%*%Beta)
}

FM_EmulatorEstimate <- function(Y){
 X <- 1:length(Y)
 EmulatorParameters <- BE_CreateDefaultEmulationClass()
 EmulatorParameters$X <- X
 EmulatorParameters$Y <- Y
 EmulatorParameters$MeanFunction <- function(X){
   X <- as.matrix(X)
   H = cbind(as.matrix(1 + 0*X) )
   return(H)
 }
 EmulatorParameters$CorrelationLength <- function(X , n){
   return(100)
 }
 EmulatorParameters$w <- function(X){
   return(0.75*var(Y)*diag(length(Y)))
 }
 EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = X , EmulatorSettings = EmulatorParameters , meanonly = 1 )
return(EmulatorOutput$E_D_fX)
 }

FM_HistoryMatchRRCulmativeDensity <- function( PriorNonImplausibleSet , x , F_x ,f_x ,  MD , RRtimes  , Corr_sdhat , imthreshold = 3 , imthreshold2 = 1.6  ){
  # Build kdeestimate for history matching  
  #m <- rollmedian( RRtimes , k = 21 , na.pad = T )
  mm <- FM_EmulatorEstimate( RRtimes )
  m <- median(RRtimes)
  RRtimes <- RRtimes + m - mm + rnorm(length(RRtimes) , 0 , 0.01 )
  RRtimes <- RRtimes[!is.na(RRtimes)]
  
  y <- matrix(FM_CalculateCDFS( RRtimes = RRtimes , xx = x ) , dim(x)[1] , dim(x)[2])

  NonZeroLogical <- matrix(T , length(RRtimes) , 1)
  #NonZeroLogical <- c( ((y > 0.03)*(y < 0.97))==1 ) ==1
  
  #if(sum(NonZeroLogical) <= 2){
  #  NonZeroLogical <- c( ((y > 0.01)*(y < 0.99))==1 ) ==1
  #}
  
  #if(sum(NonZeroLogical) <= 2){
  #  NonZeroLogical <- c( ((y > 0.001)*(y < 0.999))==1 ) ==1
  #}
  
  StdResidVector <- abs((F_x - y) / MD) 
  Im2 <- rowMeans(StdResidVector)
  Im <- apply(StdResidVector , 1 , max)
  
  #Im2 <- colMeans( apply( as.matrix(F_x[ , NonZeroLogical]) , 1 , function(X){abs(X - y[ NonZeroLogical])} )/ t(MD[ , NonZeroLogical] + 0.0045)  , na.rm = T)
  #Im <- apply( apply( as.matrix(F_x[ , NonZeroLogical]) , 1 , function(X){abs(X - y[ NonZeroLogical])} ) / t(MD[ , NonZeroLogical] + 0.0045) , 2 , function(X){max(X , na.rm = T)} )
  
  RPeakKDEEstmate <- kde( RRtimes )
  z <- predict(RPeakKDEEstmate , x = x[which.min(Im2) , ])
  
  LogicalVector <- (Im <= imthreshold)&(Im2 <  imthreshold2)

  #if(sum(LogicalVector) > 1){
  
  #Im3 <- matrix(0 , sum(LogicalVector) , 1)
  #Im4 <- matrix(0 , sum(LogicalVector) , 1)
  
  #for(i in 1:dim(Im3)[1]){
    #covMattmp <- DP_AddNugget( ((MD[which(LogicalVector)[i] , NonZeroLogical] + 0.0045) %*% t(MD[which(LogicalVector)[i] , NonZeroLogical]+ 0.0045)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[which(LogicalVector)[i] , NonZeroLogical]+ 0.0045) )
  #  Im3[i , ] <- (abs((z[NonZeroLogical] - f_x[which(LogicalVector)[i] , NonZeroLogical])/f_x[which(LogicalVector)[i] , NonZeroLogical]  ))[which.max(z[NonZeroLogical]) ]
  #  Im4[i , ] <- sum(log(FM_EvaluateDenistyEstimate(RRtimes , PriorNonImplausibleSet[which(LogicalVector)[i] , ])   ))
    #Im3[i, ] <- mahalanobis(x = y[NonZeroLogical] , center = F_x[which(LogicalVector)[i] , NonZeroLogical] , cov = covMattmp  ) 
  #}
  #}
  #if(sum(LogicalVector) <= 1){
  #  Im3<- (abs( (z[NonZeroLogical] - f_x[which.min(Im) , NonZeroLogical])/f_x[which.min(Im) , NonZeroLogical] ) )[which.max(z[NonZeroLogical])]
   # Im4<- sum(log(FM_EvaluateDenistyEstimate(RRtimes , PriorNonImplausibleSet[which.min(Im) , ])   ))
    
    #covMattmp <- DP_AddNugget( ((MD[which.min(Im) , NonZeroLogical] + 0.0045) %*% t(MD[which.min(Im) , NonZeroLogical]+ 0.0045)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[which.min(Im)  , NonZeroLogical]+ 0.0045) )
    #Im3 <- mahalanobis(x = y[NonZeroLogical] , center = F_x[which.min(Im) , NonZeroLogical] , cov = covMattmp  ) 
  #}
  
  if(sum(LogicalVector) > 1){
    return( setNames(list( PriorNonImplausibleSet[LogicalVector,] , cbind(Im[LogicalVector] , Im2[LogicalVector] ) , F_x[LogicalVector ,  ] , y[which.min(Im2) , ]  ) , c('NonImplausibleSets' , 'Implausability' , 'f_x' , 'y'  )  ) )
  }else{
    return( setNames(list( PriorNonImplausibleSet[which.min(Im),] , cbind(min(Im ) , Im2[which.min(Im)] ) , F_x[which.min(Im) ,  ] , y[which.min(Im2),]  ) , c('MinImplausiblePoint' , 'Implausability' , 'f_x' , 'y'  ) ) )
  }
}
FM_SampleRealisationsSet <- function( PriorNonImplausibleSet , N , x ){
  output <- matrix( 0 , dim(PriorNonImplausibleSet)[1] , dim(x)[2] )
  for(i in 1:dim(PriorNonImplausibleSet)[1]){
    output[i , ] <-  FM_CalculateCDFS( RRtimes = FM_SampleGMM( X = PriorNonImplausibleSet[i,] , N = N )  , xx = x[i,]) 
    #output[i , ] <-  FM_CalculateCDFS( RRtimes = FM_EvaluateDenistyEstimate( x = x , X = PriorNonImplausibleSet[i,] )  , xx = x) 
  }
  return(output)
}
FM_CalulateImForGroundTruth <- function(x , F_x , f_x , PriorNonImplausibleSet , MD , Corr_sdhat , N = 231  ){
  
  ImplausabilityMatrix <- matrix(0 , dim(PriorNonImplausibleSet)[1] , 4)
  
  for( ii in 1:dim(PriorNonImplausibleSet)[1] ){
    {  
    #RRtimes <- FM_SampleGMM(X = PriorNonImplausibleSet[ii , ] ,  N)
    G_0 <- function(N){ FM_SampleGMM( X = PriorNonImplausibleSet[ii,] , N) }
   
    RRtimes <- FM_SampleDP(c , l , N , G_0)
    y <- FM_CalculateCDFS( RRtimes = RRtimes , xx = x[ii,] )
    #lines(x , y)
    
    NonZeroLogical = matrix(T , length(y) , 1) 
   # NonZeroLogical <- c( ((y > 0.03)*(y < 0.97))==1 ) ==1
  #  if(sum(NonZeroLogical) <= 2){
   #   NonZeroLogical <- c( ((y > 0.01)*(y < 0.99))==1 ) ==1
  #  }
   # if(sum(NonZeroLogical) <= 2){
  #    NonZeroLogical <- c( ((y > 0.001)*(y < 0.999))==1 ) ==1
  #  }
    
    
    #mean( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical] + 0.003) , na.rm = T)
    #max( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical]+ 0.003) , na.rm = T)
    #covMattmp <- DP_AddNugget( ((MD[ii , NonZeroLogical]+ 0.003) %*% t(MD[ii , NonZeroLogical]+ 0.003)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[ii , NonZeroLogical]+ 0.00000001) )
    #mahalanobis(x = y[NonZeroLogical] , center = F_x[ii , NonZeroLogical] , cov = covMattmp  ) 
    
    ImplausabilityMatrix[ii , 1] <- mean( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical] + 0.0045) , na.rm = T)
    ImplausabilityMatrix[ii , 2] <- max( abs(y[NonZeroLogical] - F_x[ii , NonZeroLogical]) / (MD[ii , NonZeroLogical]+ 0.0045) , na.rm = T)
    #ImplausabilityMatrix[ii , 4] <- sum(log(FM_EvaluateDenistyEstimate(RRtimes ,PriorNonImplausibleSet[ii,] )))
    #covMattmp <- DP_AddNugget( ((MD[ii , NonZeroLogical] + 0.0045) %*% t(MD[ii , NonZeroLogical]+ 0.0045)) * Corr_sdhat[NonZeroLogical , NonZeroLogical]  , 0.0001*diag(MD[ii , NonZeroLogical]+ 0.0045) )
    #ImplausabilityMatrix[ii , 3] <- mahalanobis(x = y[NonZeroLogical] , center = F_x[ii , NonZeroLogical] , cov = covMattmp  ) 
    
    #RPeakKDEEstmate <- kde( RRtimes )
    #z <- predict(RPeakKDEEstmate , x = x )
    
    #ImplausabilityMatrix[ii , 3] <- (abs( (z[NonZeroLogical] - f_x[ii , NonZeroLogical])/f_x[ii , NonZeroLogical] ) )[which.max(z[NonZeroLogical])]
    }
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
  # Look up table for large input length(xx)
  if(length(xx) > 10000){
    minxx = min(xx)
    maxxx = max(xx)
    rangexx = maxxx - minxx
    Lookupinputs <- seq(minxx , maxxx , (rangexx)/9999 )
    LookupValues <- FM_CalculateCDFS( RRtimes = RRtimes , xx = Lookupinputs  )
    tmp <-  ( (rangexx)/9999 )
    LookupPoints <- round( (xx - minxx) / tmp) + 1
    output <- LookupValues[LookupPoints]
    return(output)
  }
  # Loop for small input length(xx)
  if(length(xx) <= 10000){
    output <- matrix(0 , length(xx) , 1)
    for(i in 1:length(xx)){
      output[i] <- sum(RRtimes <= xx[i])/length(RRtimes)   
    }
    return(output)
  }
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
  
  colnames(d) <- c('Feature' , 'Note' , 'Non-implausible' , 'Median' , 'IQR' )
  
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
  
  Start_Xstar <- Xstar[1]
  #EmulatedQS <- matrix( 0 , 100 , length(Xstar) )
  counter <- 1
  #print( 'Emulating TQ Segments.' )
  for(i in 1:dim(QS_Struct$Date)[1]){
    #if(QS_Struct$numvalues[i ] > as.numeric( quantile( QS_Struct$numvalues  , upper)) ){next}
    #if(QS_Struct$numvalues[i ] < as.numeric( quantile( QS_Struct$numvalues  , lower)) ){next}
    tmp = 1:(min(dim(QS_Struct$Date)[2] , QS_Struct$numvalues[i ]) -1)
    EmulatorParameters$X <- (1:length(QS_Struct$Date[i,tmp]))/length(QS_Struct$Date[i,tmp])
    
    EmulatorParameters$Y <- QS_Struct$Value[i,tmp]
    EmulatorParameters$Y <- EmulatorParameters$Y[EmulatorParameters$X >= max(Start_Xstar[1] , EmulatorParameters$X[which.min(abs(QS_Struct$Date[i,tmp] - 0.275) ) ] ) ]
    EmulatorParameters$X <- EmulatorParameters$X[EmulatorParameters$X >= max(Start_Xstar[1] , EmulatorParameters$X[which.min(abs(QS_Struct$Date[i,tmp] - 0.275) ) ] ) ]
    
    if(length(  EmulatorParameters$X ) < 5){next}
    Xstar = seq(from = min(EmulatorParameters$X) , to = max(EmulatorParameters$X)   , length.out = length(Xstar) )
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar , EmulatorSettings = EmulatorParameters  , meanonly = 1)
    EmulatedQS[counter,] <- EmulatorOutput$E_D_fX
    counter <- counter +1
    #if(counter > 100){break}
  }  
  #print('TQ Segments Emulated.' )
  EmulatedQS <- EmulatedQS[-(counter:dim(EmulatedQS)[1] ), ]
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
FM_SampleDP <- function(c , l , n , G_0){
  
  b <- rbeta(n, 1, c)
  p <- numeric(n)
  p[1] <- b[1]
  p[2:n] <- sapply(2:n, function(i) b[i] * prod(1 - b[1:(i-1)]))
  y <- G_0(n)
  theta <- sample(y, prob = p, replace = TRUE) 
  A <- sqrt(var(theta)/(var(theta) + l^2))
  theta <- A*(theta + rnorm(n = n  , mean = 0 , l )) + (1 - A)*mean(theta)
  return(theta)
}
PWaveHM_CreateDefaultEmulationclass<- function(){ 
  EmulatorParameters <- BE_CreateDefaultEmulationClass()
  EmulatorParameters$w <- function(X){
    return( 0.5*diag(length(X) ) )
  }
  EmulatorParameters$CorrelationLength<- function(X , n){
    return(0.01)
  }
  return(EmulatorParameters)
}
FM_GMMCulDenisty <- function(x , Pi =1 , mu =0 , sigma =1){
  f <- matrix(0 , length(x) , 1)
  for(i in 1:length(Pi)){
    f <- f + Pi[i]*pnorm(x , mean = mu[i] , sd = sigma[i])  
  }
  return(f)
}
FM_cdflaplace <- function(x , mu = 0 , sigma = 1 , alpha = 0){
  return(alpha*pnorm(x , mu , sigma) + (1-alpha)*plaplace(x , mu , sigma/sqrt(2)) )
}
FM_EvalulateCDFEstimate <-function(x , X){
  a <-  FM_GMMCulDenisty(x = x , Pi = X[c(1,3) ] , mu = X[c(4,6) ] , sigma = X[c(7,9) ] ) 
  b <- X[2]*FM_cdflaplace(x =x , mu = X[5] , sigma = X[8] , alpha = X[10] ) 
  return(a + b)
}
FM_ExtractActiveOutputs <- function(y , x){
  
  NonZeroLogical <- c( ((y > 0.01)*(y < 0.99))==1 ) ==1
  if(sum(NonZeroLogical) <= 2){
    NonZeroLogical <- c( ((y > 0.001)*(y < 0.999))==1 ) ==1
  }
  
  seq(min(x[NonZeroLogical]) , max(x[NonZeroLogical]) , (max(x[NonZeroLogical]) - min(x[NonZeroLogical]))/length(x) )
  return(seq(min(x[NonZeroLogical]) , max(x[NonZeroLogical]) , (max(x[NonZeroLogical]) - min(x[NonZeroLogical]))/length(x) )[1:201])
  
}
FM_SampleGMMBigeminy <- function(X , N = 250 ){
  output <- matrix(0 , N , 1)
  for(i in 1:N){
    if(mod(i , 2) == 1){
      output[i , ] <- rnorm(1 , mean = X[4] , sd = (X[7]) )
    }  
    if(mod(i , 2) == 0){
      output[i , ] <- rnorm(1 , mean = X[6] , sd = (X[9]) )
    }  
  }
  return(output)
}
FM_GetNonImplausibleSetsFromlabelledDataset <- function(){
  
  path <- list()
  path <- choose.dir( caption  = paste0( "Select folder containing NSR data repository" ))
  listAllPatients <- as.matrix(list.dirs(path = path[[1]], full.names = FALSE, recursive = FALSE))
  
  
  SetOfNonImplausibleSets <- matrix(0 , 0 , 10)
  
  for(j in 1:length(listAllPatients)){
    load(paste0("D:\\nsrdb_outputs\\HMOutput" , listAllPatients[[j]] , '.RData'))
    for( i in 1:length(outputstruct)){
      SetOfNonImplausibleSets <-  rbind(SetOfNonImplausibleSets , outputstruct[[i]][[2]]$NonImplausibleSets )
    }
  }
  return(unique( SetOfNonImplausibleSets ))   
}
FM_ltafdbExtractStartandEndAF <- function(RPeakData , AFlocations , minutethreshold =5){
  
  AFAnotation <- BC_CreateAnnotationFromInference(RPeakData$RRCombined$t , data.frame(Start = AFlocations[,1] , End = AFlocations[,2]))
  StartandEndAF <- AFD_GetStartEndAF(t = RPeakData$RRCombined$t , logicaltimeseries = (AFAnotation == 1) , minutethreshold = minutethreshold)
  return(StartandEndAF)
}
FM_ltafCheckSegmentisAF <- function(timestamp ,  AFLogical , RPeakData , n = 500){
  
  return(AFLogical[ min(which.min(abs(timestamp - RPeakData$RRCombined$t) ) + n , length(AFLogical)) ] ==1)
  
}