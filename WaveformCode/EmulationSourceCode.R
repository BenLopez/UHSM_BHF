##### Covariance functions #####

CF_ExponentialFamily <- function(x , xstar , l , p){
  # Expenential family correlation length family
  options(warn=-1)  
  # Turn into matrices
x <- as.matrix(x)
xstar <- as.matrix(xstar)
l <- as.matrix(l^p[1])
p <- as.matrix(p)

sizex <- dim(x)
sizexstar <- dim(xstar)


distmatrix = array(0, c(sizex[1] , sizexstar[1] , sizex[2]))

if(sizex[2] > 1){

  for (i in 1:sizex[2])
{
  distmatrix[  ,  , i] <- (( as.matrix(pdist( as.matrix((x[,i])) , as.matrix( (xstar[,i]) )  )) )^p[1])/l[i]
}

} else
{
  distmatrix <- (( as.matrix(pdist( as.matrix((x[,1])) , as.matrix( (xstar[,1]) )  )) )^p[1])/l[1]
}
if(sizex[2] > 1)
{
distmatrix <- rowSums(distmatrix , dims = 2)
}

KXXstar <- exp( -0.5 * distmatrix )
return(KXXstar)
options(warn=0)
}
CF_NeuralNetworkSingle <- function(x , xstar , Sigma){
  tildex <- cbind( matrix( 1 , dim(x)[1] , 1) , x )
  tildexstar <- cbind( matrix( 1 , dim(xstar)[1] , 1), xstar )
  return((2/pi)*asin(  ((2*(tildex%*%Sigma%*%t(tildexstar)) ) / sqrt(( 1 + 2*tildex%*%Sigma%*%t(tildex) )*( 1 + 2*tildexstar%*%Sigma%*%t(tildexstar)))) ))
  
}
CF_NeuralNetwork<- function(x , xstar , Sigma = diag(c(10,10))){
  x <- as.matrix(x)
  xstar <- as.matrix(xstar)
  KXX <- matrix(0 , dim(x)[1] ,dim(xstar)[1] )
  # Vectorise this horror loop!
  for(i in 1:dim(x)[1]){
    for(j in 1:dim(xstar)[1]){
      KXX[i,j] <- CF_NeuralNetworkSingle(x = as.matrix(x[i,]) ,
                                         xstar = as.matrix(xstar[j,]) ,
                                         Sigma = Sigma )
      
    }
  }
  return(KXX)
}

##### Emulation functions #####

BE_CreateDefaultEmulationClass <- function(){
  EmulatorParameters <- setNames(list(1 , 1 , 1, 1 , 1 , 1 ) , c('MeanFunction' , 'CorrFunction' , 'CorrelationLength' , 'w' , 'X' , 'Y'))
  EmulatorParameters$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X) , X )
    return(H)
  }
  EmulatorParameters$CorrFunction <- function(x , xstar , l){
    return(CF_ExponentialFamily(x , xstar , l , 2))
  }
  EmulatorParameters$CorrelationLength <- function(X , n){
    return((DP_FindMeanMindistances(X)))
  }
  EmulatorParameters$w <- function(X){
    return(1)
  }
  return(EmulatorParameters)
} 
BE_EmulationClassAddData <- function(y , x , EmulatorSettings = BE_CreateDefaultEmulationClass() , HistoryMatchSettings = BE_CreateDefaultHistoryMatchClass()){
  if(HistoryMatchSettings$AddPointsToEmulator == 1 & length( EmulatorSettings$X ) > 1){
    EmulatorSettings$X <- rbind(as.matrix(EmulatorSettings$X) , as.matrix(x))
    EmulatorSettings$Y <- rbind(as.matrix(EmulatorSettings$Y) , as.matrix(y))
  }
  if(HistoryMatchSettings$AddPointsToEmulator == 0 || length( EmulatorSettings$X ) == 1 ){
  EmulatorSettings$X <- x
  EmulatorSettings$Y <- y}
  return(EmulatorSettings)
}

BE_BayesLinearEmulatorGLSEstimates <- function( xstar , EmulatorSettings = BE_CreateDefaultEmulationClass() ){
  
  x <- as.matrix(EmulatorSettings$X)
  y <- as.matrix(EmulatorSettings$Y)
  n <- dim(x)[1]
  l <- EmulatorSettings$CorrelationLength( x , 1 )
  w <- EmulatorSettings$w(x)
  
  KXX <- EmulatorSettings$CorrFunction(x , x , l )
  KXstarX <-  EmulatorSettings$CorrFunction( xstar , x , l )
  KXstarXstar <-  EmulatorSettings$CorrFunction(xstar , xstar , l )
  H <- EmulatorSettings$MeanFunction(x)
  sizeh = dim(H)
  Hstar <- EmulatorSettings$MeanFunction(xstar)
  sizex <- dim(x)[1]
  
  
  if(length(w) > 1){
  D <- KXX + w%*%diag(sizex[1])
  }else{
  D <- KXX + w*diag(sizex[1])
  }
  
  L = t(chol(D))
  
  alpha <- solve(L , H)
  Q <- t(chol(t(alpha)%*%alpha))
  # Calulate approximated adjusted expectation of Beta
  Betastar <- solve(t(Q) , solve(Q))%*%t(alpha)%*%solve( L , y)
  
  zeta <- solve(L , y - H%*%Betastar)
  
  # Calulate approximated adjusted expectation of sigma
  Sigmastar <- 1/(sizex[1] - sizeh[2]) * ( t(zeta)%*%zeta )
  
  omega <-solve(L, t(KXstarX))
  kappa <-  solve( Q , t( Hstar - t(omega)%*%alpha ))
  
  E_z_y <- Hstar%*%Betastar + KXstarX %*% solve(t(L) , zeta) 
  V_z_y <- Sigmastar[1]*(  KXstarXstar  - t(omega)%*%omega + t(kappa)%*%kappa )
  
  return(list(E_D_fX = E_z_y , V_D_fX = V_z_y))

}
BE_BayesLinearEmulatorLSEstimates<- function(xstar , EmulatorSettings = BE_CreateDefaultEmulationClass() , meanonly = 0 ){
  
  x <- as.matrix(EmulatorSettings$X)
  y <- as.matrix(EmulatorSettings$Y)
  n <- dim(x)[1]
  l <- EmulatorSettings$CorrelationLength( x , 1 )
  w <- EmulatorSettings$w(x)
  
  KXX <- EmulatorSettings$CorrFunction(x , x , l )
  KXstarX <-  EmulatorSettings$CorrFunction( xstar , x , l )
  if(meanonly == 0){
  KXstarXstar <-  EmulatorSettings$CorrFunction(xstar , xstar , l )}
  H <- EmulatorSettings$MeanFunction(x)
  sizeh = dim(H)
  Hstar <- EmulatorSettings$MeanFunction(xstar)
  sizex <- dim(x)[1]
  
  BetaHat <- solve(t(H)%*%H)%*%t(H)%*%y
  epsilon <- y - H%*%BetaHat
  
  SigmaHat <- 1/( sizeh[1] - sizeh[2] -1 ) * t(epsilon)%*%epsilon
  SigmaHat <- as.numeric(SigmaHat - median(EmulatorSettings$w(x)[EmulatorSettings$w(x) !=0] ))

  Var_D <- solve(as.numeric(SigmaHat)*KXX + EmulatorSettings$w(x)%*%diag(sizeh[1]))
  Cov_XstarD <- as.numeric(SigmaHat)*KXstarX

  E_z_y <- Hstar%*%BetaHat +  Cov_XstarD%*%Var_D%*%epsilon
  
  if(meanonly == 0){
  V_z_y <- SigmaHat*KXstarXstar - (Cov_XstarD%*%Var_D%*%t(Cov_XstarD))
  }
  if(meanonly == 1){
  V_z_y <- 0  
  }
  
  return(list(E_D_fX = E_z_y , V_D_fX = V_z_y , E_D_MX = Hstar%*%BetaHat))  
}

BE_BayesLinearEmulatorLSEstimatesBatchMode<- function(xstar , EmulatorSettings = BE_CreateDefaultEmulationClass()){
  
  
  if(BE_CheckifPreCalulationhasBeenPerformed(EmulatorSettings) == FALSE){
    EmulatorSettings <-  BE_PerformPreCalulationForLSEEmulator(EmulatorSettings)
  }
  
  x <- as.matrix(EmulatorSettings$X)
  y <- as.matrix(EmulatorSettings$Y)
  n <- dim(x)[1]
  l <- EmulatorSettings$CorrelationLength( x , 1 )
  w <- EmulatorSettings$w(x)
  KXX <- EmulatorSettings$KXX
  H <- EmulatorSettings$H
  sizeh <- EmulatorSettings$sizeh
  sizex <-  EmulatorSettings$sizex
  
  BetaHat <- EmulatorSettings$BetaHat
  epsilon <- EmulatorSettings$epsilon
  SigmaHat <-EmulatorSettings$SigmaHat
  Var_D <- EmulatorSettings$Var_D
  E_z_y <- matrix(0 , dim(xstar)[1] , 1)
  V_z_y <- matrix(0 , dim(xstar)[1] , 1)
  batches <- seq(1 , dim(xstar)[1]  , 1000)
  if(batches[length(batches)] != dim(xstar)[1]){ batches = c( batches , dim(xstar)[1])}
  
  for(i in 1:(length(batches) - 1)){
  
  Hstar <- EmulatorSettings$MeanFunction(xstar[batches[i]:batches[i+1] , ] )
  KXstarX <-  EmulatorSettings$CorrFunction( xstar[batches[i]:batches[i+1] , ] , x , l )
  Cov_XstarD <- SigmaHat*KXstarX
  
  E_z_y[batches[i]:batches[i+1] , ] <- Hstar%*%BetaHat +  Cov_XstarD%*%Var_D%*%epsilon
  V_z_y[batches[i]:batches[i+1] , ] <- (SigmaHat - diag(Cov_XstarD%*%Var_D%*%t(Cov_XstarD)) )
  
  if(length(batches) > 1000){
  DP_WaitBar(i/(length(batches) - 1))
  }
  }
  
  return(list(E_D_fX = E_z_y , V_D_fX = V_z_y))  
}
BE_BayesLinearEmulatorWithInputUncertainty <- function(xstar , v_xstar , EmulatorSettings , numbersamples = 1000 ){
  SampleofPoints <- BE_SampleLHSinab(numbersamples , a = xstar - 2*sqrt(v_xstar) , b = xstar + 2*sqrt(v_xstar)   )
  emulatoroutput <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar = SampleofPoints , EmulatorSettings = LocalEmulationClass)
  output <- setNames(list( mean(emulatoroutput$E_D_fX) , var(emulatoroutput$E_D_fX) + mean(emulatoroutput$V_D_fX) ) ,  c('E_D_fX' , 'V_D_fX' ))
  return(output)
}

BE_PerformPreCalulationForLSEEmulator <- function(EmulatorSettings){
  x <- as.matrix(EmulatorSettings$X)
  y <- as.matrix(EmulatorSettings$Y)
  n <- dim(x)[1]
  l <- EmulatorSettings$CorrelationLength( x , 1 )
  w <- EmulatorSettings$w(x)
  
  KXX <- EmulatorSettings$CorrFunction(x , x , l )
  EmulatorSettings$KXX <- KXX
  H <- EmulatorSettings$MeanFunction(x)
  EmulatorSettings$H <- H
  sizeh = dim(H)
  EmulatorSettings$sizeh <- sizeh
  sizex <- dim(x)[1]
  EmulatorSettings$sizex <- sizex
  
  BetaHat <- solve(t(H)%*%H)%*%t(H)%*%y
  EmulatorSettings$BetaHat <- BetaHat
  epsilon <- y - H%*%BetaHat
  EmulatorSettings$epsilon <- epsilon
  SigmaHat <- 1/( sizeh[1] - sizeh[2] -1 ) * t(epsilon)%*%epsilon
  SigmaHat <- as.numeric(SigmaHat - median(EmulatorSettings$w(x)[EmulatorSettings$w(x) !=0]))
  EmulatorSettings$SigmaHat <- SigmaHat
  
  Var_D <- solve(SigmaHat*KXX + EmulatorSettings$w(x)%*%diag(sizeh[1]))
  EmulatorSettings$Var_D <- Var_D
  
  
  if(SigmaHat < 0 ){warning('Negative variance! Please check sigmaHat calulation')}
  
  return(EmulatorSettings)
  
}
BE_CheckifPreCalulationhasBeenPerformed <- function(EmulatorSettings){
  setoffields <- c("MeanFunction",
                   "CorrFunction",
                   "CorrelationLength",
                   "w",
                   "X",
                   "Y",
                   "KXX",
                   "H",
                   "sizeh",
                   "sizex",
                   "BetaHat",
                   "epsilon",
                   "Var_D",
                   "SigmaHat") 
  
  if(sum(names(EmulatorSettings) %in% setoffields) == length(setoffields) ){output <- TRUE}
  if(sum(names(EmulatorSettings) %in% setoffields) != length(setoffields) ){output <- FALSE}
  return(output)
}

#### History Matching Functions #####
BE_CalculateImEmulatorTrainingSet <- function(X , Simulator , ImMeasure){
  X <- as.matrix(X)
  fX <- matrix(0 , dim(X)[1] , 1)
  
  for(ii in 1:dim(X)[1]){
    fX[ii,] <- ImMeasure( Simulator(X[ii,])  )   
    #DP_WaitBar(ii/dim(X)[1])  
  }
  return(fX)
}
BE_CreateDefaultHistoryMatchClass <- function( ){
  HistoryMatchSettings <- setNames(list( 1 , 1 , 1 , 1 , 1 ,1) , c('n_emulatortrainingpoints' , 'n_emulatorsamplepoints' , 'p' , 'AddPointsToEmulator' , 'Im_Thresh' , 'KDE'))
  HistoryMatchSettings$n_emulatortrainingpoints <- 10
  HistoryMatchSettings$n_emulatorsamplepoints <- 10000
  HistoryMatchSettings$p <- 1
  HistoryMatchSettings$AddPointsToEmulator <- 0
  HistoryMatchSettings$Im_Thresh <- 3
  HistoryMatchSettings$KDE <- 0
  return( HistoryMatchSettings )
}
BE_CalulateNonImplausibleSets <- function(EmulatorSettings = BE_CreateDefaultEmulationClass() , HistoryMatchSettings , PriorRange , chi_star , WaveCounter){
  
  n_emulatortrainingpoints = HistoryMatchSettings$n_emulatortrainingpoints
  n_emulatorsamplepoints = HistoryMatchSettings$n_emulatorsamplepoints
  p = HistoryMatchSettings$p
  
  if(HistoryMatchSettings$KDE == 0 || length(chi_star) ==1 ){
  X <- DP_RescaleZeroOneToab(randomLHS(n_emulatortrainingpoints , p) , PriorRange[1] , PriorRange[2])
  }
  
  #if(HistoryMatchSettings$KDE == 1 & length(chi_star) > 1){
  #X <- BC_SampleGMM( MclustDistributionStruct =  KDE_CreateMclustClassFromSample(X = as.matrix(chi_star) , H = 0.00001*sqrt( var(chi_star)) ) ,numberofsamples =  n_emulatortrainingpoints)
  #X <- X[ X >0 ]
  #}
  
  print('Evaluating simulator.')
  Im <- BE_CalculateImEmulatorTrainingSet(X , Simulator , ImMeasure)
  #if(WaveCounter  == 1){
  #  X = rbind(0,X)
  #  Im = rbind(0,Im)}
  print('Simulator evaluated.')
  
  EmulatorSettings <- BE_EmulationClassAddData( y = Im , x = X , EmulatorSettings = EmulatorSettings , HistoryMatchSettings = HistoryMatchSettings)
 
  if(HistoryMatchSettings$KDE == 0|| length(chi_star) ==1){
  Xstar <- sort(rbind(DP_RescaleZeroOneToab(randomLHS(n_emulatorsamplepoints , p) , PriorRange[1] , PriorRange[2]) , PriorRange[1] ))
  }
  #if(HistoryMatchSettings$KDE == 1 & length(chi_star) > 1){
  #Xstar <- BC_SampleGMM( MclustDistributionStruct =  KDE_CreateMclustClassFromSample(X = as.matrix(chi_star) , H = 0.00001*sqrt( var(chi_star)) ) ,numberofsamples =  n_emulatorsamplepoints)
  #Xstar <- sort(Xstar[Xstar > 0])
  #}
  print( 'Evaluating emulator.' )
  Em_output <- BE_BayesLinearEmulatorLSEstimates(xstar  = Xstar , EmulatorSettings = EmulatorSettings)
  print( 'Emulator evaluated.' )
  
  BE_PlotOneDOutput(Em_output = Em_output ,  EmulatorSettings = EmulatorSettings , Xstar = Xstar )
  abline(3 , 0)
  if( sum((Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh )  > 0  || sum( (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh )  > 0 ){
    if( sum((Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh )  > 0){
      chi_star <- Xstar[ (Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)))  < HistoryMatchSettings$Im_Thresh  ]    }
    if(sum( (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX)))<HistoryMatchSettings$Im_Thresh )  == 0 & sum(( Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)))<HistoryMatchSettings$Im_Thresh )  == 0){  
      chi_star <- NA}
  
    }
  return(list('chi_star' = chi_star , 'EmulatorSettings' = EmulatorSettings))  
}
BE_HistoryMatch <- function(TrainingSet , TrainingSet2, EmulatorSettings = BE_CreateDefaultEmulationClass() , HistoryMatchSettings = BE_CreateDefaultHistoryMatchClass() , PriorRange = c(0,1) ){
  
  print('History Matching')
  SizechiStariMinus1 <- 1
  WaveCounter <- 1
  chi_star <<- matrix(0 , 1 , 1)
  while( length(chi_star) >=  SizechiStariMinus1 ){
    
    SizechiStariMinus1 <- length(chi_star)
    
    tmp <- BE_CalulateNonImplausibleSets(EmulatorSettings = EmulatorSettings , 
                                         HistoryMatchSettings = HistoryMatchSettings , 
                                         PriorRange = PriorRange,
                                         chi_star = chi_star , 
                                         WaveCounter = WaveCounter)
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
  EmulatorSettings$CorrelationLength <- function(X , n){
    return(DP_FindMeanMindistances(EmulatorSettings$X))
  }
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
  
  return(list('chi_star' = chi_star , 'EmulatorSettings' = EmulatorSettings))  
}

##### Plotting function #####

BE_PlotOneDOutput <- function(Em_output ,  EmulatorSettings , Xstar ){
  x11(20, 14)
  plot(EmulatorSettings$X , EmulatorSettings$Y  ,
       xlim = c(min(Xstar) , max(Xstar)) ,
       ylim = c(min(Em_output$E_D_fX) , max(Em_output$E_D_fX)),
       xlab = 'x' ,
       ylab = 'f(x)', 
       main = 'Emulator Performance')
  points(EmulatorSettings$X , EmulatorSettings$Y + 2*sqrt(diag(EmulatorSettings$w(X))) , col = 'blue')
  points(EmulatorSettings$X , EmulatorSettings$Y - 2*sqrt(diag(EmulatorSettings$w(X))) , col = 'blue')
  lines(Xstar ,  Em_output$E_D_fX  , col = 'red')
  lines(Xstar ,  Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  lines(Xstar ,  Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  if(exists('HistoryMatchSettings')){
  abline(HistoryMatchSettings$Im_Thresh,0)
  }
}
BE_PlotStdResiduals <- function(zstar , emulatoroutput , e=0){
  
  StdResids = (zstar - emulatoroutput$E_D_fX)/sqrt(diag(emulatoroutput$V_D_fX) + e)
  
  output <- ggplot( data.frame(Im_Crit = emulatoroutput$E_D_fX , StdResids = StdResids) , aes(Im_Crit , StdResids) ) +
    geom_point( color = 'blue') + ggtitle('Standardised Residuals') + geom_hline( yintercept  = 3)+ geom_hline( yintercept  = -3)
  
  x11()
  print(output)
  return(output)
}

#BE_PlotOneDOutput <- function(Em_output ,  EmulatorSettings , Xstar ){
 # x11(20, 14)
  #par(mfrow = c(1 , 2))
  #plot(EmulatorSettings$X , EmulatorSettings$Y  ,
  #     xlab = 'x' ,
  #     ylab = 'f(x)', 
  #     main = 'Emulator Performance')
  #lines(Xstar ,  Em_output$E_D_fX  , col = 'red')
  #lines(Xstar ,  Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  #lines(Xstar ,  Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  #plot(EmulatorSettings$X , EmulatorSettings$Y  ,
   #    xlim = c(min(Xstar) , max(Xstar)) ,
  #     ylim = c(min(Em_output$E_D_fX) , max(Em_output$E_D_fX)),
   #    xlab = 'x' ,
  #     ylab = 'f(x)', 
    #   main = 'Emulator Performance')
  #lines(Xstar ,  Em_output$E_D_fX  , col = 'red')
  #lines(Xstar ,  Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  #lines(Xstar ,  Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  
#}
 

BE_CreateOpticalDensityplot2D <- function(PriorSample , ImplausabilityLogical ){
  
  output <- ggplot(  ) +
    geom_point(data = data.frame(x = PriorSample[ImplausabilityLogical == 0 ,1],
                                 y = PriorSample[ImplausabilityLogical == 0 ,2]),
               aes(x , y) , colour = 'red' , alpha = 0.05) +
    geom_point(data = data.frame(x = PriorSample[ImplausabilityLogical == 1  ,1],
                                 y = PriorSample[ImplausabilityLogical == 1  ,2]),
               aes(x , y) , colour = 'blue' , alpha = 0.05) 
  return(output)            
}

##### Sampling functions #####
BE_SampleGP <- function(KXX){
  L = chol(KXX + 0.0000000000001 * diag(dim(KXX)[1]) )
  return(t(L) %*% rnorm(dim(KXX)[1]))
}
BE_SampleSeparableMVGP <- function( KXX , Sigma , L = NA ){
  if( is.na(L) ){
    L1 <-  chol( KXX + 0.0000000001*diag(dim(KXX)[1]) ) 
    L2 <-  chol( Sigma  ) 
    L <-   kronecker(L1 , L2)
  }
  output <- t(matrix(t(L) %*% rnorm( dim(L)[1] ) , dim(Sigma)[1] , dim(KXX)[1] ))
  return( output )
}
BE_SampleLHSinab <- function(a , b, numbersamples = 1000 ){
  # Function to sample a latin hypercube on ranges definned by (a) and (b)
  output <- randomLHS(numbersamples , length(a))
  for(i in 1:length(a)){
    output[ , i] <- DP_RescaleZeroOneToab(X = output[ , i] , a = a[i] , b = b[i])  
  }
  return(output)
}
BE_SampleTP <- function(KXX , df = 5,mean = 0 ){
output <-t(rmvt(n = 1 , sigma = ((df - 2)/(df))*(KXX + 0.0000000000001 * diag(dim(KXX)[1])) , df = df , delta = as.matrix(rep(mean , dim(KXX)[1]))))
return(output)
}
BE_BayesLinearEmulatorLSEstimatesMO <- function(xstar , EmulatorSettings = BE_CreateDefaultEmulationClass() , meanonly = 0 ){
  
  x <- as.matrix(EmulatorSettings$X)
  y <- as.matrix(EmulatorSettings$Y)
  n <- dim(x)[1]
  l <- EmulatorSettings$CorrelationLength( x , 1 )
  w <- EmulatorSettings$w(x)
  
  KXX <- EmulatorSettings$CorrFunction(x , x , l )
  KXstarX <-  EmulatorSettings$CorrFunction( xstar , x , l )
  if(meanonly == 0){
  KXstarXstar <-  EmulatorSettings$CorrFunction(xstar , xstar , l )}
  
  H <- EmulatorSettings$MeanFunction(x)
  sizeh = dim(H)
  Hstar <- EmulatorSettings$MeanFunction(xstar)
  sizex <- dim(x)[1]
  
  BetaHat <- solve(t(H)%*%H)%*%t(H)%*%y
  epsilon <- y - H%*%BetaHat
  
  SigmaHat <- 1/( sizeh[1] - sizeh[2] -1 ) * t(epsilon)%*%epsilon
  
  Var_D <- solve(KXX + EmulatorSettings$w(x)%*%diag(sizeh[1]))
  Cov_XstarD <- KXstarX
  
  E_z_y <- Hstar%*%BetaHat +  Cov_XstarD%*%Var_D%*%epsilon
  
  if(meanonly == 0){
    V_z_y <- KXstarXstar - (Cov_XstarD%*%Var_D%*%t(Cov_XstarD))
  }
  if(meanonly == 1){
    V_z_y <- 0  
  }
  
  return(list(E_D_fX = E_z_y , C_D_fX = V_z_y , SigmaHat = SigmaHat))  
}
