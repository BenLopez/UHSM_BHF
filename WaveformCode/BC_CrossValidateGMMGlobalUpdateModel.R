DataSetPriorProbabilities <- BC_ExtractValidationPriors(Priorprobabilities , DataBaseMaster , DataBase)


ArchiveDataBase <-  DataBase
CrossValidatedProbability <- matrix(0,dim(DataBase[[3]])[1] , 1)
for(i in 1:dim(DataBase[[3]])[1]){

DataBase[[3]] <- DataBase[[3]][-i , ]
GlobalDistributionStructtmp  <-  BC_EstimateGlobalDensitiesGMM( DataBase = DataBase , numberofcomponents = BCParameters$GlobalNumberComponentsforGMM)
CrossValidatedProbability[i,] <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix(ArchiveDataBase[[3]][i,])) , GlobalDistributionStruct = GlobalDistributionStructtmp , Priorprobabilities = DataSetPriorProbabilities)$A
DataBase <- ArchiveDataBase

DP_WaitBar(i/dim(DataBase[[3]])[1])
}

hist( CrossValidatedProbability , col=rgb(1,1,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(0 , tmp$breaks, 1) )
