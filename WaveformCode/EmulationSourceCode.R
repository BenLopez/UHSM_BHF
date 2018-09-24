##### Covariance functions #####

CF_ExponentialFamily <- function(x , xstar , l , p){
  # Expenential family correlation length family
  
  # Turn into matrices
x <- as.matrix(x)
xstar <- as.matrix(xstar)
l <- as.matrix(l^p[1])
p <- as.matrix(p)

sizex <- dim(x)
sizexstar <- dim(xstar)


distmatrix = matrix(0, sizex[1] , sizexstar[1] , sizex[2])

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
distmatrix < - apply(distmatrix , 3 , sum)
}

KXXstar <- exp( -0.5 * distmatrix )
return(KXXstar)
}

##### Inference functions #####

BE_CreateDefaultEmulationClass <- function(){
  EmulatorParameters <- setNames(list(1 , 1 , 1, 1 , 1 , 1 ) , c('MeanFunction' , 'CorrFunction' , 'CorrelationLength' , 'w' , 'X' , 'Y'))
  EmulatorParameters$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X) , X , X^2)
    return(X)
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
  n <- dim(X)[1]
  l <- EmulatorSettings$CorrelationLength( x , 1 )
  w <- EmulatorSettings$w(X)
  
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
BE_CalulateNonImplausibleSets <- function(EmulatorSettings = BE_CreateDefaultEmulationClass() , HistoryMatchSettings , PriorRange , chi_star){
  
  n_emulatortrainingpoints = HistoryMatchSettings$n_emulatortrainingpoints
  n_emulatorsamplepoints = HistoryMatchSettings$n_emulatorsamplepoints
  p = HistoryMatchSettings$p
  
  if(HistoryMatchSettings$KDE == 0 || length(chi_star) ==1 ){
  X <- DP_RescaleZeroOneToab(randomLHS(n_emulatortrainingpoints , p) , PriorRange[1] , PriorRange[2])
  }
  if(HistoryMatchSettings$KDE == 1 & length(chi_star) > 1){
  X <- BC_SampleGMM( MclustDistributionStruct =  KDE_CreateMclustClassFromSample(X = as.matrix(chi_star) , H = 0.00001*sqrt( var(chi_star)) ) ,numberofsamples =  n_emulatortrainingpoints)
  X <- X[ X >0 ]
  }
  
  print('Evaluating simulator.')
  Im <- BE_CalculateImEmulatorTrainingSet(X , Simulator , ImMeasure)
  print('Simulator evaluated.')
  
  EmulatorSettings <- BE_EmulationClassAddData( y = Im , x = X , EmulatorSettings = EmulatorSettings , HistoryMatchSettings = HistoryMatchSettings)
 
  if(HistoryMatchSettings$KDE == 0|| length(chi_star) ==1){
  Xstar <- sort(rbind(DP_RescaleZeroOneToab(randomLHS(n_emulatorsamplepoints , p) , PriorRange[1] , PriorRange[2]) , PriorRange[1] ))
  }
  if(HistoryMatchSettings$KDE == 1 & length(chi_star) > 1){
  Xstar <- BC_SampleGMM( MclustDistributionStruct =  KDE_CreateMclustClassFromSample(X = as.matrix(chi_star) , H = 0.00001*sqrt( var(chi_star)) ) ,numberofsamples =  n_emulatorsamplepoints)
  Xstar <- sort(Xstar[Xstar > 0])
  }
  print( 'Evaluating emulator.' )
  Em_output <- BE_BayesLinearEmulatorGLSEstimates(xstar  = Xstar , EmulatorSettings = EmulatorSettings)
  print( 'Emulator evaluated.' )
  
  BE_PlotOneDOutput(Em_output = Em_output ,  EmulatorSettings = EmulatorSettings , Xstar = Xstar )
  abline(3 , 0)
  if( sum((Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh )  > 0  || sum( (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh )  > 0 ){
    if( sum((Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh )  > 0){
      chi_star <- Xstar[ (Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)))  < HistoryMatchSettings$Im_Thresh  ]    }
    if(sum( (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh )  > 0){
      chi_star <- Xstar[ (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh  ]   }  
    if(sum( (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX)))<HistoryMatchSettings$Im_Thresh )  > 0 & sum(( Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)))<HistoryMatchSettings$Im_Thresh )  > 0){
      chi_star <- rbind(as.matrix(Xstar[ (Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)))  < HistoryMatchSettings$Im_Thresh  ])  ,as.matrix( Xstar[ (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX))) < HistoryMatchSettings$Im_Thresh  ] ))
    }
    if(sum( (Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX)))<HistoryMatchSettings$Im_Thresh )  == 0 & sum(( Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)))<HistoryMatchSettings$Im_Thresh )  == 0){  
      chi_star <- NA}
  
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
  points(EmulatorSettings$X , EmulatorSettings$Y + 2*sqrt(EmulatorSettings$w(X)) , col = 'blue')
  points(EmulatorSettings$X , EmulatorSettings$Y - 2*sqrt(EmulatorSettings$w(X)) , col = 'blue')
  lines(Xstar ,  Em_output$E_D_fX  , col = 'red')
  lines(Xstar ,  Em_output$E_D_fX + 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  lines(Xstar ,  Em_output$E_D_fX - 2*sqrt(diag(Em_output$V_D_fX)) , col = 'blue')
  abline(HistoryMatchSettings$Im_Thresh,0)
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
