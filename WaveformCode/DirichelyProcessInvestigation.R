
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
FilestoProcess = DP_ChooseECGstoProcess()

source('FM_CreateRhythumPriors.R')

PNIS <- PriorNonImplausibleSetRegularyIreRegular
F_x <- F_x_ReIre
f_x <- f_x_ReIre
  
{

{
  c <- 10000
  l <- 0.005
  n <- 500 - 19
  G_0 <- function(n) { FM_SampleGMM( X = PNIS[index,] , n) } 
}
  
  x11()  
  index <-1
  plot(x  , F_x[index,] , type ='l' , ylab = 'Culmulative Probability' , xlab = 'RR' , xlim = c(0,2) )
  title( 'Samples of CDFs for z')
  lines(x , F_x[index,]+ 3*sqrt( (F_x[index,]*(1-F_x[index,]))/(c + 1) + (F_x[index,]*(1-F_x[index,]))/(n) + 0.005^2) , type ='l'  )
  lines(x , F_x[index,]- 3*sqrt( (F_x[index,]*(1-F_x[index,]))/(c + 1) + (F_x[index,]*(1-F_x[index,]))/(n) + 0.005^2) , type ='l'  )
  

for(i in 1:100){
y <- G_0(n)
theta <- FM_SampleDP(c , l , n , G_0)
  
A <- sqrt(var(theta)/(var(theta) + l^2))
theta <- A*(theta + rnorm(n = n  , mean = 0 , l )) + (1 - A)*mean(theta)
#lines(x , FM_CalculateCDFS(y , xx = x) , type = 'l' , col =rgb(1,0,0,alpha = 0.05))
lines(x , FM_CalculateCDFS(theta , xx = x) , type = 'l' , col =rgb(0,0,1,alpha = 0.05))
}
}

x11()
par(mfrow = c(2 , 2))

{
 
{
  c <- 10000
  l <- 0.005
  G_0 <- function(n) { FM_SampleGMM( X = PNIS[index,] , n)}
  n <- 5000 - 19
  
}
  index <- 1
  plot(x   , f_x[index,] , type ='l' , ylim = c(0,max(f_x[index,])+0.5 ) , xlab = 'RR' , ylab = 'Density', xlim = c(0.25,2))
  title(TeX(paste0('Samples of DP: $\\alpha = ' , c  , '$')))
  
  
musigmastruct <- matrix(0 , 100 , 4)

for(i in 1:10){
  y <- G_0(n)
  theta <- FM_SampleDP(c , l , n , G_0)
  
  kde1 <- kde(theta)
  kde2 <- kde(y)
  
  musigmastruct[i , 3] <- mean(theta)
  musigmastruct[i , 4] <- var(theta)
  
  lines(x , predict( object = kde1 , x = x- 0.005 ) , type = 'l' , col =rgb(0,0,1,alpha = 0.1))
  #lines(x - 0.005 , predict( object = kde2 , x = x- 0.005 ) , type = 'l' , col =rgb(1,0,0,alpha = 0.1))
}
}


