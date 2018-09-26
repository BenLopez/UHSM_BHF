Mcluststruct <- LocalDistributionStruct[[1]] 

#Mcluststruct$parameters$mean <- Mcluststruct$parameters$mean[1:2 , ]
#Mcluststruct$parameters$variance$sigma <- Mcluststruct$parameters$variance$sigma[1:2 ,1:2, ]

Trainingset <- BC_SampleGMM( LocalDistributionStruct[[1]] , 10000 )
Trainingset2 <- BC_SampleGMM( LocalDistributionStruct[[2]] , 10000 )

Validationset <- rbind( BC_SampleGMM( LocalDistributionStruct[[1]] , 1000 ) ,  BC_SampleGMM( LocalDistributionStruct[[2]] , 1000 ))

H <- 0.1*cov(Trainingset)
tildeH <- 0.075*cov(Trainingset2)

blah <- matrix(0 , dim(Validationset)[1] , 1)
for(i in 1:dim(Validationset)[1]){
  blah[i,] <-  (sum(mahalanobis(Trainingset , center  = Validationset[i,] , cov = H) < 28)/dim(Trainingset)[1]) / ((sum(mahalanobis(Trainingset , center  = Validationset[i,] , cov = H) < 28)/dim(Trainingset)[1])+(sum(mahalanobis(Trainingset2 , center  = Validationset[i,] , cov = tildeH) < 28)/dim(Trainingset2)[1]))
}
plot(blah)
title(mean(blah[!is.na(blah)]))

0.999*(( 2*pi*det(H) )^(-0.5))

x <- t(as.matrix(Validationset[3,]))
       
KDE_DensityEstimate(  X = Trainingset 
                    , Xstar =  x
                    , H = H) / BC_PredictGMMDensity(Mcluststruct , x)
