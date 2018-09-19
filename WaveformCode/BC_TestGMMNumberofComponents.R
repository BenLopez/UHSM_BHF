
NumSamples <- 100000
SampleIndexs = sample(1:length(DataBase[[1]][ , variabletoview]) , min(NumSamples , length(DataBase[[1]][ , variabletoview])))
SampleIndexs2 = sample(1:length(DataBase[[2]][ , variabletoview]) , min(NumSamples , length(DataBase[[2]][ , variabletoview])))


BC_plotValidateDensityEstimationMarginalHistograms(DataBase[[1]] , DataBase[[2]] , DataBase[[1]][SampleIndexs1 , ] ,  DataBase[[2]][SampleIndexs2 , ] )


DataBase2 <- DataBase
DataBase2[[1]] <- DataBase[[1]][SampleIndexs1 , ]  
DataBase2[[2]] <- DataBase[[2]][SampleIndexs2 , ]   

NumberofComponents <- c(5 , 10 , 20 , 30 , 40 , 50 , 100 , 200 , 400)

for(kk in 8:length(NumberofComponents)){
  
  LocalDistributionStruct <-  BC_EstimateLocalDensitiesGMM( DataBase2 , numberofcomponents = NumberofComponents[kk] )
  save( LocalDistributionStruct , file = paste0(precomputedfolderpath, '\\old\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]]  , BCOptions[[4]], NumberofComponents[kk] ,  'DistributionSummaries'  ,'.RData')) )

}


for(kk in 1:length(NumberofComponents)){
  load( paste0(precomputedfolderpath, '\\old\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]]  , BCOptions[[4]], NumberofComponents[kk] ,  'DistributionSummaries'  ,'.RData')) )
  BC_plotValidateDensityEstimationMarginalHistograms( A = DataBase2[[1]] , B = DataBase2[[2]]  ,  C = BC_SampleGMM(MclustDistributionStruct = LocalDistributionStruct[[1]], numberofsamples = 100000) , D = BC_SampleGMM(LocalDistributionStruct[[2]], 100000) )
}  
