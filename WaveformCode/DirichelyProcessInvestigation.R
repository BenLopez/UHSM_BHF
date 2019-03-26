
{


{
  c <- 10000
  l <- 0.05
  G_0 <- function(n) { FM_SampleGMM( X = PriorNonImplausibleSetRegularyIreRegular[index,] , n)}
  n <- 500
  
}
  
  x11()  
  index <-1
  plot(x , F_x_ReIre[index,] , type ='l' , ylab = 'RR' , xlab = 'Culmulative Probability' )
  title( 'Samples of CDFs for z')
  lines(x , F_x_ReIre[index,]+ 2*sqrt( (F_x_ReIre[index,]*(1-F_x_ReIre[index,]))/(c +1) + (F_x_ReIre[index,]*(1-F_x_ReIre[index,]))/(500) + 0.005^2) , type ='l'  )
  lines(x , F_x_ReIre[index,]- 2*sqrt( (F_x_ReIre[index,]*(1-F_x_ReIre[index,]))/(c +1) + (F_x_ReIre[index,]*(1-F_x_ReIre[index,]))/(500) + 0.005^2) , type ='l'  )
  

for(i in 1:10){
theta <- FM_SampleDP(c , l , n , G_0)
  
A <- sqrt(var(theta)/(var(theta) + l^2))
theta <- A*(theta + rnorm(n = n  , mean = 0 , l )) + (1 - A)*mean(theta)

lines(x , FM_CalculateCDFS(theta , xx = x) , type = 'l' , col =rgb(0,0,1,alpha = 0.05))
lines(x , FM_CalculateCDFS(y , xx = x) , type = 'l' , col =rgb(1,0,0,alpha = 0.05))
}
}

{index <- 1
plot(x  , f_x_ReIre[index,] , type ='l' , ylim = c(0,3) , xlab = 'RR' , ylab = 'Probability Density Function')
title('Samples of Density Functions')

{
  c <- 100000
  l <- 0.01
  G_0 <- function(n) { FM_SampleGMM( X = PriorNonImplausibleSetRegularyIreRegular[index,] , n)}
  n <- 500
  
}

musigmastruct <- matrix(0 , 100 , 4)

for(i in 1:10){
  y <- G_0(n)
  theta <- FM_SampleDP(c , l , n , G_0)
  
  kde1 <- kde(theta)
  kde2 <- kde(y)
  
  musigmastruct[i , 3] <- mean(theta)
  musigmastruct[i , 4] <- var(theta)
  
  lines(kde1$eval.points , kde1$estimate , type = 'l' , col =rgb(0,0,1,alpha = 0.1))
  lines(kde2$eval.points , kde2$estimate , type = 'l' , col =rgb(1,0,0,alpha = 0.1))
}
}


