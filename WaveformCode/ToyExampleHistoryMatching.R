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
  set.seed(1)
}


#Functions
{
CalculateConjugatePosterior <- function(E_mu , V_mu , X , n , sigma){
  
  x_bar <- mean(X)
  mu_X <- ((n/sigma)*x_bar + (1/V_mu)*E_mu)/((n/sigma) +(1/V_mu)) 
  V_X <-((n/sigma) +(1/V_mu))^(-1)
  return( data.frame( mu_X = mu_X , V_X = V_X ) )
}
CalculateAdjustedExpectations <- function(E_mu , V_mu , X , n , sigma){
  E_D <- E_mu
  V_D <- V_mu + 1/n*sigma
  E_D_mu <- E_mu + (V_mu/V_D)*(mean(X) - E_D)
  V_D_mu <- V_mu - (V_mu/V_D)*(V_mu)
  return(data.frame(E_D_mu = E_D_mu , V_D_mu = V_D_mu))
}
SampleBayesLinear <- function(PriorDist = rnorm , ObsDist =rnorm , numbersample = 10000 , doplot = 1){
  
  BayesLinearAdjustedVersions<- matrix(0,numbersample , 1)
  
  for(i in 1:numbersample ){
    
    mu <- PriorDist( n = 1 , mean = E_mu , sd = sqrt(V_mu) )
    X <- ObsDist( n = n , mean = mu , sd = sqrt(sigma) )
    BayesLinearAdjustedVersions[i]<-CalculateAdjustedExpectations(E_mu , V_mu , X , n , sigma)$E_D_mu - mu
    
  }  
  
  if(doplot == 1){
  x11()
  print(hist(BayesLinearAdjustedVersions, breaks = 100, freq = F , col = rgb(0,0,1,alpha = 0.5) , ylim = c(0,1)))
  print(abline( v = 0 ))
  print(abline( v = 3*sqrt(CalculateAdjustedExpectations(E_mu , V_mu , X , n , sigma)$V_D_mu)  , col = 'red' ))
  print(abline( v = -3*sqrt(CalculateAdjustedExpectations(E_mu , V_mu , X , n , sigma)$V_D_mu) , col = 'red' ))
  }
  return(BayesLinearAdjustedVersions)
  
}
SampleBayes <- function(PriorDist = rnorm , ObsDist =rnorm , numbersimulations = 1000000, doplot = 1){
  mu <- PriorDist( n = 1 , mean = E_mu , sd = sqrt(V_mu) )
  musaved <- mu
  X <- ObsDist( n = n , mean = mu , sd = sqrt(sigma) )
  
  #X <- SampleTdis( n = n , mu = mu , V = sigma )
  X_bar <- mean(X)
  
  PosteriorMoments <- CalculateConjugatePosterior(E_mu , V_mu , X , n , sigma)
  
  if(doplot == 1){
    x11()
    x <- seq(-5 , 5 , 0.01)
    plot(x ,  dnorm(x = x , mean = PosteriorMoments$mu_X , sd = sqrt(PosteriorMoments$V_X))  , type='l', col = 'red', xlab = 'x' , ylab = 'Density')
    lines(x ,   dnorm(x = x , mean = E_mu , sd = sqrt(V_mu) ) , type = 'l' )
    title('Prior and Posterior')
    abline(v = mu)
  }
  
  
  #MulitpleMumuSamples <- as.matrix(SampleTdis( n = numbersimulations , mu = E_mu , V = V_mu ))
  MulitpleMumuSamples <- as.matrix(PriorDist( n = numbersimulations , mean = E_mu , sd = sqrt(V_mu) ))
  
  
  MulitpleDataSamples <- rowMeans(t(apply(MulitpleMumuSamples , 1 , function(X){ObsDist(n , X , sqrt(sigma) ) } )))
  #MulitpleDataSamples <- rowMeans(t(apply(MulitpleMumuSamples , 1 , function(X){SampleTdis(n , X , sigma ) } )))}
  
  if(doplot == 1){
    hist(MulitpleMumuSamples[abs(MulitpleDataSamples - X_bar) < 0.01 ], breaks = 100, freq = F , add = T )
  }
  return( data.frame(edf=FM_CalculateCDFS(MulitpleMumuSamples[abs(MulitpleDataSamples - X_bar) < 0.01 ] , seq(-4,4,0.1)) , cdf = pnorm(seq(-4,4,0.1) , mean = PosteriorMoments$mu_X , sd = sqrt(PosteriorMoments$V_X))  ))
}
CalculateKolmorgorovtextstatistic <- function(A,B){
  return(max(abs(A - B)))
}

}

E_mu <- 0
V_mu <- 1
sigma <- 5
n <- 20

KewnessKurtosisPoints <- BE_SampleLHSinab(a = c(1.8,-10) , b = c(40,10) , 1000)
KewnessKurtosisPoints <- KewnessKurtosisPoints[KewnessKurtosisPoints[,2] > -sqrt((KewnessKurtosisPoints[,1] -0.9) - 0.01) , ]
KewnessKurtosisPoints <- KewnessKurtosisPoints[KewnessKurtosisPoints[,2] < +sqrt((KewnessKurtosisPoints[,1] -0.9) - 0.01) , ]


# Full Bayes Correct
test <- SampleBayes(PriorDist = rnorm , ObsDist =rnorm , numbersimulations = 1000000, doplot = 1)
CalculateKolmorgorovtextstatistic(test$edf , test$cdf)

BayesPerformance <- matrix(0 , 100,1 )
for(i in 1:100){
PearsonsSampler <- function(n , mean , sd){rpearson(n , moments = c(mean,sd^2,KewnessKurtosisPoints[i,2],KewnessKurtosisPoints[i,1]))}
PearsonsSampler2 <- function(n , mean , sd){rpearson(n , moments = c(mean,sd^2,KewnessKurtosisPoints[100+i,2],KewnessKurtosisPoints[100+i,1]))}
test <- SampleBayes(PriorDist = PearsonsSampler , ObsDist = PearsonsSampler2  , doplot = 0)
BayesPerformance[i] <- CalculateKolmorgorovtextstatistic(test$edf , test$cdf)
DP_WaitBar(i/100)
}
BayesPerformance[is.na(BayesPerformance)] <- mean(BayesPerformance, na.rm = T)
# Bayes Linear

DataforRegression <- data.frame(X1 = KewnessKurtosisPoints[1:100,1], X2 =KewnessKurtosisPoints[1:100,2] , X3 =KewnessKurtosisPoints[101:200,1] , X4 =KewnessKurtosisPoints[101:200,2] , Y =BayesPerformance )
model = lm( Y~ X1+I(X1^2) + X2 + I(X2^2) +X3+I(X3^2)+X4+I(X4^2),data = DataforRegression)
summary(model)

x11()
par(mfrow = c(2,4))

plot(KewnessKurtosisPoints[1:100,1] , BayesPerformance , xlab = 'Prior Kurtosis', ylab = 'Max EDF Distance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))
plot(KewnessKurtosisPoints[1:100,2] , BayesPerformance , xlab = 'Prior Skewness', ylab = 'Max EDF Distance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))
plot(KewnessKurtosisPoints[101:200,1] , BayesPerformance, xlab = 'Observation Kurtosis', ylab = 'Max EDF Distance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))
plot(KewnessKurtosisPoints[101:200,2] , BayesPerformance , xlab = 'Observation Skewness', ylab = 'Max EDF Distance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))


BayesLinearPerformance <- matrix(0 , 100,1 )
for(i in 1:100){
PearsonsSampler <- function(n , mean , sd){rpearson(n , moments = c(mean,sd^2,KewnessKurtosisPoints[i,2],KewnessKurtosisPoints[i,1]))}
PearsonsSampler2 <- function(n , mean , sd){rpearson(n , moments = c(mean,sd^2,KewnessKurtosisPoints[100+i,2],KewnessKurtosisPoints[100+i,1]))}
test <- SampleBayesLinear(PriorDist = PearsonsSampler, ObsDist = PearsonsSampler2 , doplot = 0)
BayesLinearPerformance[i] <- var(test , na.rm = T) 
DP_WaitBar(i/100)
}

plot(KewnessKurtosisPoints[1:100,1] , BayesLinearPerformance , xlab = 'Prior Kurtosis', ylab = 'Sample Adjusted Variance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))
plot(KewnessKurtosisPoints[1:100,2] , BayesLinearPerformance , xlab = 'Prior Skewness', ylab = 'Sample Adjusted Variance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))
plot(KewnessKurtosisPoints[101:200,1] , BayesLinearPerformance, xlab = 'Observation Kurtosis', ylab = 'Sample Adjusted Variance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))
plot(KewnessKurtosisPoints[101:200,2] , BayesLinearPerformance , xlab = 'Observation Skewness', ylab = 'Sample Adjusted Variance',pch = 16 , col = rgb(0,0,1,alpha = 0.5))


DataforRegression <- data.frame(X1 = KewnessKurtosisPoints[1:100,1], X2 =KewnessKurtosisPoints[1:100,2] , X3 =KewnessKurtosisPoints[101:200,1] , X4 =KewnessKurtosisPoints[101:200,2] , Y = BayesLinearPerformance )
model = lm( Y~ X1+I(X1^2) + X2 + I(X2^2) +X3+I(X3^2)+X4+I(X4^2),data = DataforRegression)
summary(model)

###### History Match   ######
xstar <- DP_RescaleZeroOneToab( rand(100000,1 ) , -10,10)
PriorNonImplausibleSet <- dnorm(xstar , E_mu , sqrt(V_mu) )
PriorNonImplausibleSet <- xstar[PriorNonImplausibleSet> 0.000001]

LB <- -3*(sqrt(var(PriorNonImplausibleSet) ) + sigma)
UB <- 3*(sqrt(var(PriorNonImplausibleSet) ) + sigma)


xstar <- seq( LB , UB , (UB-LB)/201 )
FStruct <- t(apply(as.matrix(PriorNonImplausibleSet) , 1 , function(X){ dnorm(xstar , mean = X , sd = sigma) } ))

Xred <- t(apply(FStruct ,1,  function(X){as.matrix(c(xstar[sort(which(abs(X)>=0.01))[1]],xstar[sort(which(abs(X)>=0.01))[length(which(abs(X)>=0.01))]]) ) } ))
Xred <- t(apply(Xred , 1 , function(X){seq(from = X[1] ,to =  X[2] ,by =  (abs( X[1] - X[2] )/201) )[1:201] }))

FStruct <- apply(cbind(PriorNonImplausibleSet ,Xred ) , 1 , function(X){
  return(  pnorm(X[2:length(X)], X[1] , sigma) )
})




