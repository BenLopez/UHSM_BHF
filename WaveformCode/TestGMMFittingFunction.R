
DataBase2 <- DataBase
DataBase2[[1]] <- BC_SampleGMM( LocalDistributionStruct[[1]] , 100000 )
DataBase2[[2]] <- BC_SampleGMM( LocalDistributionStruct[[2]] , 100000 )

LocalDistributionStruct2 <- BC_EstimateLocalDensitiesGMM(  DataBase2 , numberofcomponents = length(LocalDistributionStruct[[1]]$parameters$pro) )

BC_plotValidateDensityEstimationMarginalHistograms(DataBase2[[1]] , DataBase2[[2]]  ,  BC_SampleGMM(LocalDistributionStruct2[[1]], 100000) , BC_SampleGMM(LocalDistributionStruct2[[2]], 100000) )
