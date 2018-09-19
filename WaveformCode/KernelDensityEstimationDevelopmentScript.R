
Trainingset <- BC_SampleGMM(LocalDistributionStruct[[1]] , 100000)
Validationset <- BC_SampleGMM(LocalDistributionStruct[[1]] , 100000)
SecondOrderImSpecifictaion <- KDE_CreateSecondOrderSpecificationImp(Trainingset = Trainingset , Validationset = Validationset  , numberofsamples = 1000)  



SampleFromNonImpH <- runif(1000 , min = 0.00000001 , max = 0.01)

Imatrix <- matrix(0 , length(SampleFromNonImpH) , 1)
for(i in 1:length(SampleFromNonImpH)){
  H = SampleFromNonImpH[i]*diag(sqrt(apply( Trainingset , 2 , var)))
  MclustDistributionStruct <- KDE_CreateMclustClassFromSample(X = Trainingset , H = H )
  Im <- KDE_CalulateImplausabilityVectorforKDE(Trainingset , MclustDistributionStruct , numberofsamples =  1000 , numberofrepitions = 1)
  
  Imatrix[i,] <- KDE_MaxIm( Im , SecondOrderImSpecifictaion )
  DP_WaitBar(i/length(SampleFromNonImpH))
}

#X2 <- BC_SampleGMM(MclustDistributionStruct , 1000)
#BC_PlotPairsFromTwoVariables( X2 , X[1:1000 , ] )
#BC_PlotPairs( X2 )
#BC_PlotPairs( X[1:1000 , ] )



