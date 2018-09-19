{

MclustDistributionStruct <- LocalDistributionStruct[[1]]

MclustDistributionStruct$parameters$pro <- c(0.2 , 0.7 ,0.1)
MclustDistributionStruct$parameters$mean <- matrix(0, 2 , 3)
MclustDistributionStruct$parameters$mean[ , 1] <- c(2 , 2)
MclustDistributionStruct$parameters$mean[ , 2] <- c(10 , 10)
MclustDistributionStruct$parameters$mean[ , 3] <- c(20 , 20)

MclustDistributionStruct$parameters$variance$sigma <- array(0 , c(2 , 2 , 3))
MclustDistributionStruct$parameters$variance$sigma[ , , 1] <- diag(c(1 , 1)) 
MclustDistributionStruct$parameters$variance$sigma[ , , 2] <- diag(c(1 , 1)) 
MclustDistributionStruct$parameters$variance$sigma[ , , 3] <- diag(c(1 , 1)) 

SampleofPoints <- BC_SampleGMM(MclustDistributionStruct = MclustDistributionStruct , numberofsamples =  1000)
pairs(SampleofPoints , col = rgb(0,0,1 , alpha = 0.2))
}
