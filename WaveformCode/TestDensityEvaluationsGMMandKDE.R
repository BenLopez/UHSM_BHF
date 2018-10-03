#Mcluststruct$parameters$mean <- Mcluststruct$parameters$mean[1:2 , ]
#Mcluststruct$parameters$variance$sigma <- Mcluststruct$parameters$variance$sigma[1:2 ,1:2, ]


alpha = c(0.5 , 0.5)

numberoftrainingsamples <- 10000
numberofvalidationsamples <- 10000

Trainingset <- BC_SampleGMM( LocalDistributionStruct[[1]] , numberoftrainingsamples )
Trainingset2 <- BC_SampleGMM( LocalDistributionStruct[[2]] , numberoftrainingsamples )

Validationset <- rbind( BC_SampleGMM(LocalDistributionStruct[[1]] , round(alpha[1]* numberofvalidationsamples) ) ,  BC_SampleGMM( LocalDistributionStruct[[2]] , round(alpha[2]*numberofvalidationsamples) ))


KDE_CalulatePuesdoProd <- function( Trainingset , x , H , thresh = 28 ){
  return(sum(mahalanobis(Trainingset , center  = x , cov = H) < 28))/length(Trainingset)
}
KDE_CalulatePuesdoPosteriorProb <- function(Trainingset ,Trainingset2 , alpha = c(0.5,0.5) , x , H , H2 , thresh = 28){
  return((alpha[1]*KDE_CalulatePuesdoProd(Trainingset , x , H , thresh )) / (alpha[1]*KDE_CalulatePuesdoProd(Trainingset , x , H , thresh ) + alpha[2]*KDE_CalulatePuesdoProd(Trainingset2 , x , H2 , thresh ) ))
}
CalculateGMMPosteriorProb <- function(LocalDistributionStruct , x , alpha = c(0.5,0.5)){
  return( (alpha[1]*BC_PredictGMMDensity(LocalDistributionStruct[[1]] , x)) / ((alpha[1]*BC_PredictGMMDensity(LocalDistributionStruct[[1]] , x)) + (alpha[2]*BC_PredictGMMDensity(LocalDistributionStruct[[2]] , x))) )
} 



Simulator <- function(X){
  H <- X*cov( rbind(Trainingset , Trainingset2) )
  H2 <-H
  
  PuesdoPosteriorProbabilities <- apply(Validationset , 1 , function(X){KDE_CalulatePuesdoPosteriorProb(Trainingset , Trainingset2 , alpha , t( as.matrix(X) ) , H , H2 , thresh = 28)} )
  
  ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(
    GlobalProbCalibrationStruct = DP_RemoveNaRows(BC_CreateProbCalibrationStruct( as.matrix(PuesdoPosteriorProbabilities) 
                                                                                  ,  alpha 
                                                                                  , numberofvalidationsamples )) 
    , BinWidth = 0.05))
  return(ProbabiliticCalibrationOutput)
  
}

ImMeasure <- function( ProbabiliticCalibrationOutput ){
  return(abs( mean((ProbabiliticCalibrationOutput$x - ProbabiliticCalibrationOutput$y)/ProbabiliticCalibrationOutput$sd)) )
}


# Define some parameters and settings for emulation and history matching.
PriorRange = c(0,0.5)

HistoryMatchSettings <- BE_CreateDefaultHistoryMatchClass()
HistoryMatchSettings$Im_Thresh <- 0.3

EmulatorSettings <- BE_CreateDefaultEmulationClass()
EmulatorSettings$w <- function(X){
  return(0.01)
}


Chi_star <- BE_HistoryMatch(TrainingSet , TrainingSet2, EmulatorSettings = EmulatorSettings  , HistoryMatchSettings = HistoryMatchSettings , PriorRange = PriorRange )
EmulatorSettings = Chi_star$EmulatorSettings
chi_star = Chi_star$chi_star

PosteriorProbabilities <- apply(Validationset , 1 , function(X){CalculateGMMPosteriorProb(LocalDistributionStruct , t( as.matrix(X) ))} )

plot( PosteriorProbabilities )
title(mean( PosteriorProbabilities ))

ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(
  GlobalProbCalibrationStruct = DP_RemoveNaRows(BC_CreateProbCalibrationStruct( as.matrix(PosteriorProbabilities) 
                                                                                ,  alpha 
                                                                                , numberofvalidationsamples )) 
  , BinWidth = 0.05))



BC_PlotCreateProbabilityCalibrationPlot(ProbabiliticCalibrationOutput) + ggtitle('Local Calibration Probabilities')

mean((ProbabiliticCalibrationOutput$x - ProbabiliticCalibrationOutput$y)/ProbabiliticCalibrationOutput$sd)

H <- 0.15*cov( Trainingset )
H2 <-H
#H2 <- 0.075*cov( Trainingset2 )

PuesdoPosteriorProbabilities <- apply(Validationset , 1 , function(X){KDE_CalulatePuesdoPosteriorProb(Trainingset , Trainingset2 , alpha , t( as.matrix(X) ) , H , H2 , thresh = 28)} )

plot( PuesdoPosteriorProbabilities )
title( mean( PuesdoPosteriorProbabilities[!is.na(PuesdoPosteriorProbabilities)] ) )


ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(
  GlobalProbCalibrationStruct = DP_RemoveNaRows(BC_CreateProbCalibrationStruct( as.matrix(PuesdoPosteriorProbabilities) 
                                                                ,  alpha 
                                                                , numberofvalidationsamples )) 
  , BinWidth = 0.05))



BC_PlotCreateProbabilityCalibrationPlot(ProbabiliticCalibrationOutput) + ggtitle('Local Calibration Probabilities')

mean((ProbabiliticCalibrationOutput$x - ProbabiliticCalibrationOutput$y)/ProbabiliticCalibrationOutput$sd)


######
















E_X <- alpha[1]
V_X <- alpha[1]*(1 - alpha[1])
E_D <- mean(ProbabiliticCalibrationOutput$y)
V_D <- var(ProbabiliticCalibrationOutput$y)

corr_DX <- cor( ProbabiliticCalibrationOutput$y ,  ProbabiliticCalibrationOutput$x)
cov_DX  <- corr_DX*sqrt( V_X*V_D )

E_D_X <- E_X + (cov_DX/V_D)*(ProbabiliticCalibrationOutput$y - E_D)
V_D_X <- V_X - (cov_DX/V_D)*cov_DX
 


ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(
  GlobalProbCalibrationStruct = DP_RemoveNaRows(BC_CreateProbCalibrationStruct( as.matrix(E_D_X) 
                                                                                ,  alpha 
                                                                                , numberofvalidationsamples )) 
  , BinWidth = 0.05))



BC_PlotCreateProbabilityCalibrationPlot(ProbabiliticCalibrationOutput) + ggtitle('Local Calibration Probabilities')
mean((ProbabiliticCalibrationOutput$x - ProbabiliticCalibrationOutput$y)/ProbabiliticCalibrationOutput$sd)
