{
  {
    if(file.exists('CheckforDefaultsScript.R')){
      source('CheckforDefaultsScript.R')
    }else{
      pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
      source("LibrariesAndSettings.R" , print.eval  = TRUE )
      DP_LoadPatientIndex()
      DP_ChooseDataReps()
      FilestoProcess <- DP_ChooseECGstoProcess() 
      HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
    }
    listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
    set.seed( 1 )
  }

}


source('FM_CreateRhythumPriors.R')

GPDen_Logistic <- function(X){
  return(1/(1+exp(-X)))
}
GPDS_ConditionalUpdate <- function(X , Xstar, Y , CorrFun , sigma){
  
  KXX <- DP_AddNugget(CorrFun(X , X) , 1e-5)
  KXstarX <- CorrFun(Xstar , X)
  
  L = t(chol(KXX))
  zeta <- solve(L , Y)
  omega <- solve(L, t(KXstarX))
  
  E_Y <- KXstarX%*%solve(t(L) , zeta)
  V_Y <- sigma - t(omega)%*%omega
  
  return(setNames(list(E_Y , V_Y) , c('E_Y' , 'V_Y')) )
}
GPDS_SampleGPDS <- function(G_0 , CorrFun , N = 500){
 
  X <- matrix(0,10*N,1)
  G <- matrix(0,10*N,1)
  D <- matrix(0,N,1)
  r <- 0
  counter <- 0
  
  while(counter < N){
    
    x_r = G_0(1)
    if(r == 0){
      gx_r <- rnorm(1 , 0 , sqrt(CorrFun(1,1)))
    }else{
      updateedGp <- GPDS_ConditionalUpdate(X = X[1:r,] , Xstar = x_r, Y = G[1:r,] , CorrFun = CorrFun , sigma = sigma)
      gx_r <- rnorm(1 , updateedGp$E_Y , updateedGp$V_Y)
    }
    
    u_r <- runif(1)
    if(u_r < GPDen_Logistic(gx_r)){
      counter <- counter + 1
      D[counter , ] <- x_r
    }
    r <- r+1
    X[r,] <-  x_r
    G[r,] <-  gx_r

  }
   return(D) 
}




# Sample from base denisty.
sigma <- 1
el = 0.2
index <- 2
numberofsamples = 200
CorrFun <- function(X , Xstar){sigma*CF_ExponentialFamily(X , Xstar, el, 2)}
G_0 <- function(N){ FM_SampleGMM( X = PriorNonImplausibleSetTotal[index,] , N) }

x <- seq(0.6,2,0.01)
plot(x , FM_EvalulateCDFEstimate(x  , PriorNonImplausibleSetTotal[index,] ) , type ='l' , xlab = 'r' , ylab = 'Culmulative Probability')
for(i in 1:50){
lines(x , FM_CalculateCDFS(FM_SampleGMM(PriorNonImplausibleSetTotal[index,] , N = numberofsamples) , x) , col = rgb(0,0,1,alpha = 0.1))
}
SampleMatrix <- matrix(0 , 50 , numberofsamples)
for(i in 1:50){
SampleMatrix[i,] <- GPDS_SampleGPDS(G_0 , CorrFun , N = numberofsamples)
lines(x , FM_CalculateCDFS(SampleMatrix[i,]  , x) , col = rgb(1,0,0,alpha = 0.1))
}


tmp <- apply(SampleMatrix , 2 , function(X){FM_CalculateCDFS(X  , x)})
tmp <- tmp - matrix(rep(FM_EvalulateCDFEstimate(x  , PriorNonImplausibleSetTotal[index,] ) , numberofsamples) , length(x) , numberofsamples )
plot(apply(tmp , 1 , mean) )

plot(apply(tmp , 1 , var) )
blah <- (FM_EvalulateCDFEstimate(x  , PriorNonImplausibleSetTotal[index,] )*(1 - FM_EvalulateCDFEstimate(x  , PriorNonImplausibleSetTotal[index,] )))
lines(blah/(numberofsamples ) + blah/( numberofsamples ) )
