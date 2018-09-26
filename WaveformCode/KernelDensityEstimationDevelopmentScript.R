
# Define dataset
Trainingset <- BC_SampleGMM(LocalDistributionStruct[[1]] , 10000)

Validationset <- BC_SampleGMM(LocalDistributionStruct[[1]] , 100000)



# Define simulator for history match.
Simulator <- function( X ){
  H <- X*diag(sqrt(apply( Trainingset , 2 , var)))
  MclustDistributionStruct <- KDE_CreateMclustClassFromSample( Trainingset , H = H )
  fx <- KDE_CalulateImplausabilityVectorforKDE(Trainingset , MclustDistributionStruct , numberofsamples =  1000 , numberofrepitions = 10)
  return(fx)
}             

# Define implausability for history match.
ImMeasure <- function(y ){
  return( KDE_MaxIm( y , SecondOrderImSpecifictaion ) )
}



# Define some parameters and settings for emulation and history matching.
PriorRange = c(0,0.1)

HistoryMatchSettings <- BE_CreateDefaultHistoryMatchClass()
HistoryMatchSettings$Im_Thresh <- 6

EmulatorSettings <- BE_CreateDefaultEmulationClass()
EmulatorSettings$w <- function(X){
  return(0.1)
}

chi_star <- KDE_HistoryMatchBandWidth( Trainingset 
                                      , Validationset 
                                      , numberofsamples = 1000 
                                      , EmulatorSettings = EmulatorSettings 
                                      , HistoryMatchSettings = HistoryMatchSettings
                                      , PriorRange = PriorRange)


# Validation plots
h_hat <- min(chi_star)
H = h_hat*diag(sqrt(apply( Trainingset , 2 , var)))
MclustDistributionStruct <- KDE_CreateMclustClassFromSample(X = Trainingset , H = H )

BC_PlotPairsFromTwoVariables( DP_SampleRowsFromMatrix(Validationset , 1000) , BC_SampleGMM(MclustDistributionStruct ,  1000) )
BC_PlotPairs( DP_SampleRowsFromMatrix(Validationset , 1000 ) )
BC_PlotPairs( BC_SampleGMM(MclustDistributionStruct , 1000 ) )

BC_PlotCompareTwoHists( DP_SampleRowsFromMatrix(Validationset , 1000) , BC_SampleGMM(MclustDistributionStruct ,  1000) )






#X2 <- BC_SampleGMM(MclustDistributionStruct , 1000)
#BC_PlotPairsFromTwoVariables( X2 , X[1:1000 , ] )
#BC_PlotPairs( X2 )
#BC_PlotPairs( X[1:1000 , ] )


