GMMmodelbyPatient <- list()

for(i in 1:dim( DataBaseMaster$NAFPatientsDatabase)[1]){
TempData1 <- DataBaseMaster$NAFPatientsDatabase[i , DataBaseMaster$NAFPatientsDatabase[i ,  ,1] !=0  ,1:11]
if(dim(TempData1)[1] < 5000){next}
DistributionSummaries <- list(1 , 1 ,1 , 1, 1, 1, 1, 1, 1, 1, 1 , 1)
for(j in 1:size(TempData1)[2]){
  DistributionSummaries[[j]] <- TempData1[ , j]  
}
TempData1 <- BC_EstimateGlobalandLocalParametersDisSum(DistributionSummaries )
GMMmodelbyPatient[[i]] <- densityMclust( DP_RemoveNaRows( TempData1$W ) , G = 100 )
DP_WaitBar(i/dim(DataBaseMaster$NAFPatientsDatabase)[1])
}

NN <- matrix(0,dim( DataBaseMaster$NAFPatientsDatabase)[1] , 1)
for(i in c(1:dim( DataBaseMaster$NAFPatientsDatabase)[1])){
  TempData1 <- DataBaseMaster$NAFPatientsDatabase[i , DataBaseMaster$NAFPatientsDatabase[i ,  ,1] !=0  ,1:11]
  NN[i,] <- dim(TempData1)[1]  
}

NN <- NN[NN >= 5000]

GMMmodelbyPatient <- DP_RemoveEmptyElementsfromlist(GMMmodelbyPatient)


BC_PlotPairs( DataBase[[2]][sample(1:dim( DataBase[[2]])[1] , 1000) , ] )
BC_PlotPairs(BC_SampleMGMM( GMMmodelbyPatient  , numberofsamples = 1000  , alpha = NN/sum(NN) ))


BC_PlotPairsFromTwoVariables( DataBase[[2]][sample(1:dim( DataBase[[2]])[1] , 1000) , ] , BC_SampleMGMM( GMMmodelbyPatient  , numberofsamples = 1000  , alpha = NN/sum(NN) ) )
BC_PlotPairsFromTwoVariables( DataBase[[2]][sample(1:dim( DataBase[[2]])[1] , 1000) , ] , BC_SampleGMM( LocalDistributionStruct[[2]]  , numberofsamples = 1000  ) )


BC_PlotCompareThreeHists(X = DataBase[[2]][sample(1:dim( DataBase[[2]])[1] , 10000) , ] ,
                         Y= BC_SampleMGMM( GMMmodelbyPatient  , numberofsamples = 10000  , alpha = NN/sum(NN) ), 
                         Z=  BC_SampleGMM( LocalDistributionStruct[[2]]  , numberofsamples = 10000  ) )




HisImSecondOrderSpecifiction <- BC_CreateSecondOrderSpecificationHistImp( Trainingset = DataBase[[2]] ,
                                                                          Validationset =  DataBase[[2]] ,
                                                                          numberofsamples = 10000)


(KDE_CalulateHistImplausability( DataBase[[2]] , BC_SampleMGMM( GMMmodelbyPatient  , numberofsamples = 10000  , alpha = NN/sum(NN) ) ) - HisImSecondOrderSpecifiction$mu)/sqrt(diag(HisImSecondOrderSpecifiction$Sigma))
(KDE_CalulateHistImplausability( DataBase[[2]] , BC_SampleGMM( LocalDistributionStruct[[2]]  , numberofsamples = 10000  ) )- HisImSecondOrderSpecifiction$mu)/sqrt(diag(HisImSecondOrderSpecifiction$Sigma))



BC_PlotCompareThreeHists <- function(X , Y , Z){
  x11(20 , 14)
  par( mfrow = c(2 , ceiling(dim(X)[2]/2)) )    
  
  for(variabletoview in 1:dim(X)[2]){
    
    tmphist <- hist(rbind(X[ , variabletoview] , Y[ , variabletoview] , Z[, variabletoview]) , breaks = 30 , plot = FALSE) 
    hist(X[ , variabletoview]
         , col=rgb(1,0,0,alpha =0.5) 
         , freq = FALSE 
         , breaks = tmphist$breaks
         , main = paste0('Variable ' , variabletoview)
         , xlabel = paste0('Variable = ' , variabletoview))
    hist(Y[ , variabletoview]
         , col = rgb(0,0,1,alpha =0.5) 
         , freq = FALSE 
         , breaks = tmphist$breaks
         , add = T)
    hist(Z[ , variabletoview]
         , col = rgb(0,1,0,alpha =0.5) 
         , freq = FALSE 
         , breaks = tmphist$breaks
         , add = T)
  }
}  

  

