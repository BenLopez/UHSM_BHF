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


SampleTdis <- function(n =1 , mu =0 , V =1  , df= 5){
  output <-t(rmvt(n = n , sigma = as.matrix(((df - 2)/(df))*(V)) , df = df , delta = as.matrix(mu) ))
  return(as.vector(output))
}
CalculateConjugatePosterior <- function(E_mu , V_mu , X , n , sigma){
  
  x_bar <- mean(X)
  mu_X <- ((n/sigma)*x_bar + (1/V_mu)*E_mu)/((n/sigma) +(1/V_mu)) 
  V_X <-((n/sigma) +(1/V_mu))^(-1)
  return( data.frame( mu_X = mu_X , V_X = V_X ) )
}



{
E_mu <- 0
V_mu <- 1
sigma <- 5
n <- 20

mu <- SampleTdis( n = 1 , mu = E_mu , V  = V_mu )
# mu <- rnorm( n = 1 , mean = E_mu , sd = sqrt(V_mu) )

# X <- rnorm( n = n , mean = mu , sd = sqrt(sigma) )
X <- SampleTdis( n = n , mu = mu , V = sigma )
X_bar <- mean(X)

PosteriorMoments <- CalculateConjugatePosterior(E_mu , V_mu , X , n , sigma)

x <- seq(-5 , 5 , 0.01)
plot(x ,  dnorm(x = x , mean = PosteriorMoments$mu_X , sd = sqrt(PosteriorMoments$V_X))  , type='l', col = 'red', xlab = 'x' , ylab = 'Density')
lines(x ,   dnorm(x = x , mean = E_mu , sd = sqrt(V_mu) ) , type = 'l' )
title('Prior and Posterior')

{
numbersimulations <- 1000000
MulitpleMumuSamples <- as.matrix(SampleTdis( n = numbersimulations , mu = E_mu , V = V_mu ))
# MulitpleMumuSamples <- as.matrix(rnorm( n = numbersimulations , mean = E_mu , sd = sqrt(V_mu) ))


# MulitpleDataSamples <- rowMeans(t(apply(MulitpleMumuSamples , 1 , function(X){rnorm(n , X , sqrt(sigma) ) } )))}
MulitpleDataSamples <- rowMeans(t(apply(MulitpleMumuSamples , 1 , function(X){SampleTdis(n , X , sigma ) } )))}


hist(MulitpleMumuSamples[abs(MulitpleDataSamples - X_bar) < 0.01 ], breaks = 100, freq = F , add = T )
  
}
