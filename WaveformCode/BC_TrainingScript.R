# Load or create master data.
if( BCOptions[[1]] == 'DistributionSummaries' ){
  if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath, paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData')) ){
    print(paste0('Pre computed database master found. If you do not want to use this delete ', paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData'), 'from'  , precomputedfolderpath , '.'))
    load(file = paste0(precomputedfolderpath , '\\' , paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData')))  
  }else{ 
    print('Pre computed database master not found. Processing.')
    DataBaseMaster <- BC_LoadInAllDistributionSummaries(path , listAllPatients , PatIndex2017)
    print('Saving file.')
    save(DataBaseMaster , file = paste0(precomputedfolderpath , '\\' , paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData')))
    print('File Saved.')
  }
}

if( BCOptions[[1]] == 'CDFs' ){
  if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath, paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData')) ){
    print(paste0('Pre computed database master found. If you do not want to use this delete ', paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData'), 'from'  , precomputedfolderpath , '.'))
    load(file = paste0(precomputedfolderpath , '\\' , paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData')))  
  }else{ 
    print('Pre computed database master not found. Processing.')
    DataBaseMaster <- BC_LoadInAllCDFS(path , listAllPatients , PatIndex2017)
    print('Saving file.')
    save(DataBaseMaster , file = paste0(precomputedfolderpath , '\\' , paste0(BCOptions[[1]] , 'DataBaseMaster'  ,'.RData')))
    print('File Saved.')
  }
}

# Load or create database.
if(BCOptions[[2]] == 'AFClassifier' ){
  if(DP_CheckfileinPrecomputedfolder(precomputedfolderpath, paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData')
  )){
    print(paste0('Pre computed database master found. If you do not want to use this delete', paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData'), ' from ' , precomputedfolderpath , '.'))
    load(paste0(precomputedfolderpath, '\\' , paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData')))
  }else{
    print('Pre computed database master not found. Processing.')
    DataBase <- BC_CreateAFandNOAFDataStructure(DataBaseMaster =  DataBaseMaster,PatIndex2017 = PatIndex2017)
for(i in 1:length(DataBase)){DataBase[[i]] <- DataBase[[i]][ , 2:31] }
  save( DataBase , file = paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData')) )
}
#source('BC_ValidationPlots.R')  
}

# Load or create database.
if(BCOptions[[2]] == 'AFPreAFClassifier' ){
  if(DP_CheckfileinPrecomputedfolder(precomputedfolderpath, paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData')
  )){
    print(paste0('Pre computed database master found. If you do not want to use this delete', paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData'), ' from ' , precomputedfolderpath , '.'))
    load(paste0(precomputedfolderpath, '\\' , paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData')))
  }else{
    print('Pre computed database master not found. Processing.')
    DataBase <- BC_CreateAFPreAFandNOAFDataStructure( DataBaseMaster =  DataBaseMaster,PatIndex2017 = PatIndex2017)
    save( DataBase , file = paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,  'DataBase'  ,'.RData')) )
  }
  #source('BC_ValidationPlots.R')  
}

# Estimate population densities.
if(BCOptions[[3]] == 'MVN'){
  GlobalDistributionStruct <-  BC_EstimateGlobalDensitiesMVN(DataBase)
  GlobalSecondOrderStruct  <-  GlobalDistributionStruct
  }
if(BCOptions[[3]] == 'GMM'){
  if(DP_CheckfileinPrecomputedfolder(precomputedfolderpath, paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]]  ,   'DistributionSummaries'  ,'.RData')
  )){
    print(paste0('Pre computed database master found. If you do not want to use this delete',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]]  ,   'DistributionSummaries'  ,'.RData') ,' from ' , precomputedfolderpath , '.'))
    load(paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]]  ,   'DistributionSummaries'  ,'.RData')))
  }
  GlobalDistributionStruct  <-  BC_EstimateGlobalDensitiesGMM( DataBase = DataBase , numberofcomponents = BCParameters$GlobalNumberComponentsforGMM)
  GlobalSecondOrderStruct   <-  BC_EstimateGlobalDensitiesMVN( DataBase )
  save( GlobalDistributionStruct , file = paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]]  ,   'DistributionSummaries'  ,'.RData')) )
}  

if(BCOptions[[4]] == 'MVN'){
  LocalDistributionStruct <-  BC_EstimateLocalDensitiesMVN(DataBase)
  LocalSecondOrderStruct <- LocalDistributionStruct
}
if(BCOptions[[4]] == 'GMM'){
  if(DP_CheckfileinPrecomputedfolder(precomputedfolderpath, paste0(BCOptions[[1]] ,BCOptions[[2]] , 'MVN'  , BCOptions[[4]] ,  'DistributionSummaries'  ,'.RData')
  )){
    print(paste0('Pre computed database master found. If you do not want to use this delete',paste0(BCOptions[[1]] ,BCOptions[[2]] ,'MVN'  , BCOptions[[4]],  'DistributionSummaries'  ,'.RData') ,' from ' , precomputedfolderpath , '.'))
    load(paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,'MVN' , BCOptions[[4]] ,   'DistributionSummaries'  ,'.RData')))
    LocalSecondOrderStruct <- BC_EstimateLocalDensitiesMVN(DataBase)
  }else
    LocalDistributionStruct <-  BC_EstimateLocalDensitiesGMM( DataBase , numberofcomponents = BCParameters$NumberComponentsforGMM )
    LocalSecondOrderStruct <- BC_EstimateLocalDensitiesMVN( DataBase )
  save( LocalDistributionStruct , file = paste0(precomputedfolderpath, '\\',paste0(BCOptions[[1]] ,BCOptions[[2]] ,BCOptions[[3]]  , BCOptions[[4]],  'DistributionSummaries'  ,'.RData')) )
}

ImSecondOrderStruct <- BC_CalulateImplausabiltySecondOrderStatistics(DataBase = DataBase , SecondOrderStruct = LocalSecondOrderStruct )
GlobalImSecondOrderStruct <- BC_GlobalCalulateImplausabiltySecondOrderStatistics(DataBase = DataBase , SecondOrderStruct = GlobalSecondOrderStruct )

DataSetPriorProbabilities <- Priorprobabilities
DataSetPriorProbabilities$A <- length(DataBaseMaster$AFPatinetsNames)/(length(DataBaseMaster$NAFPatinetsNames)+length(DataBaseMaster$AFPatinetsNames))
DataSetPriorProbabilities$`A^c` <- 1 - DataSetPriorProbabilities$A

if( BCOptions$GlobalUpdate == 'Yes' ){
  {
    GlobalUpdateDiagnostics <- matrix(0 , dim(DataBase[[3]])[1] + dim(DataBase[[4]])[1] , 1 )
    counter <-1
    
if(BCOptions$DensityEstimationGlobal == 'MVN'){
    for(i in 1:dim(DataBase[[3]])[1]){
      GlobalUpdateDiagnostics[counter,] <- BC_GlobalBayesianBeliefUpdateMVN(M = DataBase[[3]][i,] , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = Priorprobabilities)$A
      counter <- counter +1
    }
    for(i in 1:dim(DataBase[[4]])[1]){
      GlobalUpdateDiagnostics[counter,] <- BC_GlobalBayesianBeliefUpdateMVN(M = DataBase[[4]][i,] , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = Priorprobabilities)$A
      counter <- counter +1
    }
}
if(BCOptions$DensityEstimationGlobal == 'GMM'){
      for(i in 1:dim(DataBase[[3]])[1]){
        if(sum(is.na(t(as.matrix(DataBase[[3]][i,]))))>0){next}
        GlobalUpdateDiagnostics[counter,] <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix(DataBase[[3]][i,])) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = Priorprobabilities)$A
        counter <- counter +1
      }
      for(i in 1:dim(DataBase[[4]])[1]){
        if(sum(is.na(t(as.matrix(DataBase[[4]][i,]))))>0){next}
        GlobalUpdateDiagnostics[counter,] <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix(DataBase[[4]][i,])) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = Priorprobabilities)$A
        counter <- counter +1
      }
}
x11()
        tmp <- hist(GlobalUpdateDiagnostics[52:dim(GlobalUpdateDiagnostics)[1], ] , col=rgb(0,0,1,alpha = 0.5) ,
                main= paste0('Global Update Diagnostics' , ' Histogram') , 
                xlab = 'Updated Probabilities' , 
                freq = FALSE , 
                ylim = c(0,100) , breaks = 50 )
    hist(GlobalUpdateDiagnostics[1:52, ], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(0 , tmp$breaks, 1))
    hist(0.2 + 0*GlobalUpdateDiagnostics, col=rgb(0,1,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(0 , tmp$breaks, 1))
  }
}
