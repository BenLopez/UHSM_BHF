# Multiple output    
{X <- matrix(c(1:500) ,500 , 1)
l <- 10
p <- 1.45
KXX <- CF_ExponentialFamily(X , X , l , p)
alpha = 0.5

mu1 <-  colMeans( DP_RemoveNaRows(DataBase[[3]]) ) 
mu2 <-  mu1
Sigma1 <-  cov( DP_RemoveNaRows(DataBase[[3]]) ) 
Sigma2 <-  1.1*cov( DP_RemoveNaRows(DataBase[[3]]) ) 

SampleGP <- t(apply(BE_SampleSeparableMVGP(KXX , Sigma2)   , 1 , function(X){X + mu2}))

Probabilities <- CD_CalculateActualUpdatedProbabilitiesMOGP(SampleGP , alpha , KXX , invD = 0 , mu1 , mu2 , Sigma1 , Sigma2 )

#SampleGP <- DataBaseMaster$AFPatientsDatabase[2 , 10001:10500 , 1:11]
test2 <- CD_CalculateIndividualDensitiesMOPG( SampleGP , mu1 , mu2 , Sigma1 , Sigma2 )
test2 <- alpha*test2[,1] / (alpha*test2[,1] + (1-alpha)*test2[,2])

Indtest <- CD_CalulateIndenpendentEstimatedProbabilitiesMO( SampleGP ,alpha = alpha, mu1 , mu2 , Sigma1 , Sigma2, weight = 0 )

plot(test2, type = 'l' , ylim = c(0,1))
abline(alpha , 0)
lines(Probabilities , type = 'l', col = 'red')
lines(Indtest , type = 'l', col = 'blue')
Specification1 <- CD_CreateDefaultSpecification1D( l , p , mu1 , Sigma1 ) 
Specification2 <- CD_CreateDefaultSpecification1D( l , p , mu2 , Sigma2 ) 
}


numberofsamples <- 100
SetofSamples <- array(0 , c( Specification1$n , length(Specification2$mu) , round(alpha*numberofsamples) + round((1-alpha)*numberofsamples)) )
SetofSamples[,,1:round(alpha*numberofsamples)] <- CD_SampleDataFromSpecification1DMO( Specification = Specification1 ,  numberofsamples = round(alpha*numberofsamples) )
SetofSamples[,,(round(alpha*numberofsamples) + 1):size(SetofSamples)[3]] <- CD_SampleDataFromSpecification1DMO( Specification = Specification2 ,  numberofsamples = round((1-alpha)*numberofsamples) ) 

TestProbabilities <- CD_CalculateUpdatedProbabilitiesFromSpecificationsMO(SetofSamples , Specification1 , Specification2 , alpha)

SOS <- CD_CalculateSecondOrderSpecificationMO( Specification1 , Specification2 , alpha = alpha  , numberofsamples = 1000 , numberinupdate = 1)

SampleGP <- t(apply(BE_SampleSeparableMVGP(KXX , Sigma1)   , 1 , function(X){X + mu1}))
Probabilities <- CD_CalculateActualUpdatedProbabilitiesMOGP(SampleGP , alpha , KXX , invD = 0 , mu1 , mu2 , Sigma1 , Sigma2 )
AdjustedBeliefPrevision <- CD_CalulateAdjustedBeliefsForPrevision1DMO( SOS , SampleGP , Specification1 , Specification2,  alpha = alpha )

{
  x11(20,14)
  p3 <- ggplot(  ) +
    geom_line(data = data.frame(time = c( 1:501 ) , Probabilities = Probabilities[,1] ) , aes(time , Probabilities) ,  color = rgb(0,0,1 , alpha = 0.5)) +
    ggtitle('Adjusted Beliefs for Actual Probabilities')
  tmp1 <- AdjustedBeliefPrevision[,1] + 3*sqrt(AdjustedBeliefPrevision[,2]) 
  tmp2 <- AdjustedBeliefPrevision[,1] - 3*sqrt(AdjustedBeliefPrevision[,2]) 
  tmp1[tmp1 > 1] <- 1
  tmp2[tmp2 < 0] <- 0
  p3 <- p3 +
    geom_line(data = data.frame( time = c( 1:501 ) , Probabilities = AdjustedBeliefPrevision[,1] ) , aes(time , Probabilities) ,  color = rgb(1,0,0 , alpha = 0.5))+
    geom_line(data = data.frame( time = c( 1:501 ) , Probabilities =  tmp1 ) , aes(time , Probabilities) ,  color = rgb(0,0,0 , alpha = 0.5))+
    geom_line(data = data.frame( time = c( 1:501 ) , Probabilities =  tmp2 ) , aes(time , Probabilities) ,  color = rgb(0,0,0 , alpha = 0.5))
  print(p3)
}




Probabilities2 <- Probabilities
for(i in 2:501){
  Probabilities2[i,] <- dnorm(x = mean(mahal[1:(i-1),1]) , mean = 11 , sd = sqrt(var(mahal[1:(i-1),1])/(i-1))) / (dnorm(x = mean(mahal[1:(i-1),1]) , mean = 11 , sd = sqrt(var(mahal[1:(i-1),1])/(i-1))) + dnorm(x = mean(mahal[1:(i-1),2]) , mean = 11 , sd = sqrt(var(mahal[1:(i-1),2])/(i-1))))
}

plot(Probabilities2 , Probabilities , pch =16 , col = rgb(1 , 0, 0 , alpha = 0.1))
