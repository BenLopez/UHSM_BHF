
Trainingset <- BC_SampleGMM(LocalDistributionStruct[[1]] , 10000)
Validationset <- BC_SampleGMM(LocalDistributionStruct[[1]] , 10000)
SecondOrderImSpecifictaion <- KDE_CreateSecondOrderSpecificationImp(Trainingset = Trainingset , Validationset = Validationset  , numberofsamples = 1000)  


SampleFromNonImpH <- runif(100 , min = 0.0004 , max = 0.0006)

Imatrix <- matrix(0 , length(SampleFromNonImpH) , 1)
for(i in 1:length(SampleFromNonImpH)){
  H = SampleFromNonImpH[i]*diag(sqrt(apply( Trainingset , 2 , var)))
  MclustDistributionStruct <- KDE_CreateMclustClassFromSample(X = Trainingset , H = H )
  Im <- KDE_CalulateImplausabilityVectorforKDE(Trainingset , MclustDistributionStruct , numberofsamples =  1000 , numberofrepitions = 10)
  
  Imatrix[i,] <- KDE_MaxIm( Im , SecondOrderImSpecifictaion )
  DP_WaitBar(i/length(SampleFromNonImpH))
}

plot(SampleFromNonImpH , Imatrix)


H = 0.00047*diag(sqrt(apply( Trainingset , 2 , var)))
MclustDistributionStruct <- KDE_CreateMclustClassFromSample(X = Trainingset , H = H )
Im <- KDE_MaxIm(KDE_CalulateImplausabilityVectorforKDE(Trainingset , MclustDistributionStruct , numberofsamples =  1000 , numberofrepitions = 1) , SecondOrderImSpecifictaion = SecondOrderImSpecifictaion)

BC_PlotPairsFromTwoVariables( DP_SampleRowsFromMatrix(Validationset , 1000) , BC_SampleGMM(MclustDistributionStruct ,  1000) )
BC_PlotPairs( DP_SampleRowsFromMatrix(Validationset , 1000 ) )
BC_PlotPairs( BC_SampleGMM(MclustDistributionStruct , 1000 ) )

BC_PlotCompareTwoHists( DP_SampleRowsFromMatrix(Validationset , 1000) , BC_SampleGMM(MclustDistributionStruct ,  1000) )


Simulator <- function( X){
  H <- X*diag(sqrt(apply( Trainingset , 2 , var)))
  MclustDistributionStruct <- KDE_CreateMclustClassFromSample( Trainingset , H = H )
  fx <- KDE_CalulateImplausabilityVectorforKDE(Trainingset , MclustDistributionStruct , numberofsamples =  1000 , numberofrepitions = 10)
  return(fx)
}             

ImMeasure <- function(y ){
  return( KDE_MaxIm( y , SecondOrderImSpecifictaion ) )
}

PriorRange = c(0,0.1)
HistoryMatchSettings <- BE_CreateDefaultHistoryMatchClass()
EmulatorSettings <- BE_CreateDefaultEmulationClass()
EmulatorSettings$CorrelationLength <- function(X , n){
  return(2*(DP_FindMeanMindistances(X)))
}
HistoryMatchSettings$KDE <- 0
EmulatorSettings$w <- function(X){
  return(0.05)
}
HistoryMatchSettings$Im_Thresh <- 6
SizechiStariMinus1 <- 1
WaveCounter <- 1
chi_star <- matrix(0 , 1 , 1)
while( length(chi_star) >=  SizechiStariMinus1 ){

SizechiStariMinus1 <- length(chi_star)

tmp <- BE_CalulateNonImplausibleSets(EmulatorSettings = EmulatorSettings , 
                                     HistoryMatchSettings = HistoryMatchSettings , 
                                     PriorRange = PriorRange,
                                     chi_star = chi_star)
chi_star <- tmp$chi_star
EmulatorSettings <- tmp$EmulatorSettings


PriorRange <- DP_CalculateLimits(chi_star)
print(paste0('Number of Non-Implausible Points = ', length(chi_star)))
print(paste0('Wave ' , WaveCounter , 'completed.'))
abline(v = PriorRange[1])
abline(v = PriorRange[2])

if(is.na(chi_star[1])){
  break
}

WaveCounter <- WaveCounter + 1

}
HistoryMatchSettings$AddPointsToEmulator = 1
chi_star <- sample( chi_star , SizechiStariMinus1 , replace = TRUE)

SizechiStariMinus1 <- 10
while(var(chi_star) <=  SizechiStariMinus1){
  
  SizechiStariMinus1 <- var(chi_star)
  
  tmp <- BE_CalulateNonImplausibleSets(EmulatorSettings = EmulatorSettings , 
                                       HistoryMatchSettings = HistoryMatchSettings , 
                                       PriorRange = PriorRange,
                                       chi_star = chi_star)
  chi_star <- tmp$chi_star
  EmulatorSettings <- tmp$EmulatorSettings
  
  PriorRange <- DP_CalculateLimits(chi_star)
  print(paste0('Number of Non-Implausible Points = ', length(chi_star)))
  print(paste0('Wave ' , WaveCounter , 'completed.'))
  WaveCounter <- WaveCounter + 1
  abline(v = PriorRange[1])
  abline(v = PriorRange[2])
  
  
  if(is.na(chi_star[1])){
    break
  }
  
}


#X2 <- BC_SampleGMM(MclustDistributionStruct , 1000)
#BC_PlotPairsFromTwoVariables( X2 , X[1:1000 , ] )
#BC_PlotPairs( X2 )
#BC_PlotPairs( X[1:1000 , ] )


