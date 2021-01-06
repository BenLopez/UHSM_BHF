{
  c <- 5000
  l <- 0.005
  n <- 501
  numberofsamples <- 100000
  
  {
    #numberofsamples <- 100000
    #PriorNonImplausibleSetRegular <- matrix(0 , numberofsamples , 10)
    
    #PriorNonImplausibleSetRegular[, 5] <- runif(numberofsamples , 0.4 , 1)
    #PriorNonImplausibleSetRegular[, 4] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] - 0.15*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] - 0.05*PriorNonImplausibleSetRegular[, 5])
    #PriorNonImplausibleSetRegular[, 6] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] + 0.05*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] + 0.15*PriorNonImplausibleSetRegular[, 5] )
    
    #PriorNonImplausibleSetRegular[, 1] <- runif(numberofsamples , 0 , 0.1)
    #PriorNonImplausibleSetRegular[, 3] <- runif(numberofsamples , 0 , 0.1)
    #PriorNonImplausibleSetRegular[, 2] <- 1 - (PriorNonImplausibleSetRegular[, 1] + PriorNonImplausibleSetRegular[, 3]) 
    
    #PriorNonImplausibleSetRegular[, 7] <- runif(numberofsamples , 0.004 , (0.06*PriorNonImplausibleSetRegular[, 4])  )
    #PriorNonImplausibleSetRegular[, 8] <- runif(numberofsamples , 0.004 , (0.06*PriorNonImplausibleSetRegular[, 5]) )
    #PriorNonImplausibleSetRegular[, 9] <- runif(numberofsamples ,  apply(as.matrix(PriorNonImplausibleSetRegular[, 7]- 0.15*PriorNonImplausibleSetRegular[, 7]) , 1 , function(X){max(0,X)} ) , 1.15*PriorNonImplausibleSetRegular[, 7] )
    
    #PriorNonImplausibleSetRegular[, 10]<- runif(numberofsamples , 0 , 1)
  
    
    PriorNonImplausibleSetRegular <- matrix(0 , numberofsamples , 10)
    
    PriorNonImplausibleSetRegular[, 5] <- runif(numberofsamples , 0.4 , 1)
    PriorNonImplausibleSetRegular[, 4] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] - 0.25*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] - 0.05*PriorNonImplausibleSetRegular[, 5])
    PriorNonImplausibleSetRegular[, 6] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] + 0.05*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] + 0.25*PriorNonImplausibleSetRegular[, 5] )
    
    PriorNonImplausibleSetRegular[, 1] <- runif(numberofsamples , 0 , 0.25)
    PriorNonImplausibleSetRegular[, 3] <- runif(numberofsamples , 0 , 0.25)
    PriorNonImplausibleSetRegular[, 2] <- 1 - (PriorNonImplausibleSetRegular[, 1] + PriorNonImplausibleSetRegular[, 3]) 
    
    PriorNonImplausibleSetRegular[, 7] <- runif(numberofsamples , 0.004 , (0.06*PriorNonImplausibleSetRegular[, 4])  )
    PriorNonImplausibleSetRegular[, 8] <- runif(numberofsamples , 0.004 , (0.06*PriorNonImplausibleSetRegular[, 5]) )
    PriorNonImplausibleSetRegular[, 9] <- runif(numberofsamples ,  apply(as.matrix(PriorNonImplausibleSetRegular[, 7]- 0.15*PriorNonImplausibleSetRegular[, 7]) , 1 , function(X){max(0,X)} ) , 1.15*PriorNonImplausibleSetRegular[, 7] )
    
    PriorNonImplausibleSetRegular[, 10]<- runif(numberofsamples , 0 , 1)
    
    if(UseAnnotatedData == 1){
      
      PriorNonImplausibleSetRegular <- unique(rbind(PriorNonImplausibleSetRegular , FM_GetNonImplausibleSetsFromlabelledDataset()))
    
    }
      
  }
  # SPecify function of measurement error and model discrepancy.
    MDFunction <- function(F_x , X,gamma = 0.005^2){
      output <- sqrt((F_x*(1-F_x))/c + (F_x*(1-F_x))/n + gamma )
      return(output)
    }
    
    #MDFunction <- function(f_x , X){
    #  return( 0.3  )
    #}
    
    x <- seq(0 , 2, 0.01)
    f_xreg <- t(apply(PriorNonImplausibleSetRegular , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
    F_xreg <- t(apply(PriorNonImplausibleSetRegular , 1 , function(X){FM_EvalulateCDFEstimate(x , X) }))
    
    xreg <- t(apply(F_xreg , 1 , function(X){FM_ExtractActiveOutputs(X , x)} ))
    f_xreg <- t(apply(cbind(PriorNonImplausibleSetRegular , xreg) , 1 , function(X){FM_EvaluateDenistyEstimate(  X[11:211] , X[1:10]) }))
    F_xreg <- t(apply(cbind(PriorNonImplausibleSetRegular , xreg) , 1 , function(X){FM_EvalulateCDFEstimate(X[11:211] , X[1:10] ) }))
    
    MD_Reg <- matrix(0 , dim(F_xreg)[1] , length(x))
    for( i in 1:dim(F_xreg)[1]){
      MD_Reg[i , ] <- MDFunction( F_xreg[i , ] , PriorNonImplausibleSetRegular[i , ] )
    }
  }
  
  # Prior specification regularly-irregular
  {
    
    {
    {  
    {PriorNonImplausibleSetRegularyIreRegular <- matrix(0 , numberofsamples , 10)
    
    #PriorNonImplausibleSetRegularyIreRegular[, 5] <- runif(numberofsamples , 0.1 , 1.5)
    #PriorNonImplausibleSetRegularyIreRegular[, 4] <- runif(numberofsamples , 0.05 ,  apply(as.matrix(0.9*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 , function(X){max(X , 0.1)}) )
    #PriorNonImplausibleSetRegularyIreRegular[, 6] <- runif(numberofsamples , apply(as.matrix(PriorNonImplausibleSetRegularyIreRegular[, 5] + 0.1*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 ,function(X){max(0.1 , X)} ) , 2)
      
    PriorNonImplausibleSetRegularyIreRegular[, 5] <- runif(numberofsamples , 0.3 , 1)
    PriorNonImplausibleSetRegularyIreRegular[, 4] <- runif(numberofsamples , 0.15 ,  apply(as.matrix(0.9*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 , function(X){max(X , 0.26)}) )
    PriorNonImplausibleSetRegularyIreRegular[, 6] <- runif(numberofsamples , apply(as.matrix(PriorNonImplausibleSetRegularyIreRegular[, 5] + 0.1*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 ,function(X){max(0.1 , X)} ) , 2)
    
    #PriorNonImplausibleSetRegularyIreRegular[, 2] <- runif(numberofsamples , 0.4 , 1)
    PriorNonImplausibleSetRegularyIreRegular[, 2] <- runif(numberofsamples , 0.4 , 1)
    PriorNonImplausibleSetRegularyIreRegular[, 1] <- runif(numberofsamples , rep(0 , numberofsamples ) ,  apply( cbind(PriorNonImplausibleSetRegularyIreRegular[, 2] , (1 - PriorNonImplausibleSetRegularyIreRegular[,2])) , 1 , min) )
    PriorNonImplausibleSetRegularyIreRegular[, 3] <- 1 - PriorNonImplausibleSetRegularyIreRegular[, 1] - PriorNonImplausibleSetRegularyIreRegular[, 2]
    
    #PriorNonImplausibleSetRegularyIreRegular[, 8] <- runif(numberofsamples , 0.05 , 0.2 )
     PriorNonImplausibleSetRegularyIreRegular[, 8] <- runif(numberofsamples , 0.03 , 0.2 )
    #PriorNonImplausibleSetRegularyIreRegular[, 7] <- runif(numberofsamples , 0.8*PriorNonImplausibleSetRegularyIreRegular[, 8] , PriorNonImplausibleSetRegularyIreRegular[, 8] )
     PriorNonImplausibleSetRegularyIreRegular[, 7] <- runif(numberofsamples , 0.2*PriorNonImplausibleSetRegularyIreRegular[, 8] , PriorNonImplausibleSetRegularyIreRegular[, 8] )
    #PriorNonImplausibleSetRegularyIreRegular[, 9] <- runif(numberofsamples , 0.8*PriorNonImplausibleSetRegularyIreRegular[, 8] , 0.2 )
     PriorNonImplausibleSetRegularyIreRegular[, 9] <- runif(numberofsamples , 0.2*PriorNonImplausibleSetRegularyIreRegular[, 8] , PriorNonImplausibleSetRegularyIreRegular[, 8] )
     PriorNonImplausibleSetRegularyIreRegular[, 10]<- runif(numberofsamples , 0 , 1)
    #BC_PlotPairsFromTwoVariables(PriorNonImplausibleSetRegularyIreRegular[1:500,7:9] , PriorNonImplausibleSetRegular[1:500,7:9]  , alpha = 0.1)  
    }
    x <- seq(0 , 2, 0.01)
    f_x_ReIre <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
    #F_x_ReIre <- t(apply(f_x_ReIre , 1 , function(X){cumsum( c(0,diff(x))*(X/sum(c(0,diff(x))*(X))) ) }))
    F_x_ReIre <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvalulateCDFEstimate(x , X) }))
    
    PeakHeights <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvaluateDenistyEstimate(X[4:6], X) }))
    
    # max height must be central peak
    Valid1 = (PeakHeights[,1] < PeakHeights[,2])*(PeakHeights[,3] < PeakHeights[,2]) == 1
    
    # No daylight between peaks
    #MinBetweenPeaks = apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){ min(c( min(FM_EvaluateDenistyEstimate(seq(X[4] , X[5] , 0.01), X)/max(FM_EvaluateDenistyEstimate(x, X)) ) , min(FM_EvaluateDenistyEstimate(seq(X[5] , X[6] , 0.01), X)/max(FM_EvaluateDenistyEstimate(x, X)) ))) } ) 
    MinBetweenPeaks = apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){ min(c( min(FM_EvaluateDenistyEstimate(seq(X[4] , X[5] , 0.01), X) ) , min(FM_EvaluateDenistyEstimate(seq(X[5] , X[6] , 0.01), X) ))) } ) 
    Valid2 = MinBetweenPeaks > 0.2
    Valid = (Valid1*Valid2) == 1
    
    PriorNonImplausibleSetRegularyIreRegular <- PriorNonImplausibleSetRegularyIreRegular[Valid , ]
    f_x_ReIre <- f_x_ReIre[Valid , ]
    F_x_ReIre <- F_x_ReIre[Valid , ]
    
    if(UseAnnotatedData == 1){
      load(paste0(precomputedfolderpath , '\\AFNonImplausibleSets.RData'))
      PriorNonImplausibleSetRegularyIreRegular <- rbind(PriorNonImplausibleSetRegularyIreRegular,AFNonImplausibleSets)
    }
    f_x_ReIre <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
    F_x_ReIre <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvalulateCDFEstimate(x , X) }))
  
    xIrIreg <- t(apply(F_x_ReIre , 1 , function(X){FM_ExtractActiveOutputs(X , x)} ))
    f_x_ReIre <- t(apply(cbind(PriorNonImplausibleSetRegularyIreRegular , xIrIreg) , 1 , function(X){FM_EvaluateDenistyEstimate(  X[11:211] , X[1:10]) }))
    F_x_ReIre <- t(apply(cbind(PriorNonImplausibleSetRegularyIreRegular , xIrIreg) , 1 , function(X){FM_EvalulateCDFEstimate(X[11:211] , X[1:10] ) }))

    MD_ReIre <- matrix(0 , dim(PriorNonImplausibleSetRegularyIreRegular)[1] , length(x))
    for( i in 1:dim(PriorNonImplausibleSetRegularyIreRegular)[1]){
      MD_ReIre[i , ] <- MDFunction(F_x_ReIre[i , ] , PriorNonImplausibleSetRegularyIreRegular[i ,])
    }
  }
  
  BC_PlotPairsFromTwoVariables(PriorNonImplausibleSetRegularyIreRegular[sample(1:dim(PriorNonImplausibleSetRegularyIreRegular)[1] , 1000 ),] ,
                               PriorNonImplausibleSetRegular[sample(100000:dim(PriorNonImplausibleSetRegular)[1] , 1000 ),]  , alpha = 0.1)  
  }
  
  # Calculate approximate correlation matrix.
  { #SampleMatrix  <-  FM_SampleRealisationsSet( PriorNonImplausibleSet = PriorNonImplausibleSetRegularyIreRegular , N = n , x = xIrIreg  )
    #ME_matrix     <-  SampleMatrix - F_x_ReIre
    #Cov_sdhat1    <-  cov( ME_matrix )
    #Corr_sdhat1   <-  cov2cor( DP_AddNugget(Cov_sdhat1 , 0.00000045) )
    Corr_sdhat1   <- 1
    
    #SampleMatrix <- FM_SampleRealisationsSet( PriorNonImplausibleSet =  PriorNonImplausibleSetRegular[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] , N = n , x =  xreg   )
    #ME_matrix <- SampleMatrix - F_xreg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],]
    #Cov_sdhat2 <- cov( ME_matrix )
    #Corr_sdhat2 <- cov2cor( DP_AddNugget(Cov_sdhat2 , 0.00000045) )
    Corr_sdhat2 <- 1
  }
  
  { 
    #ImplausabilityMatrix <-  FM_CalulateImForGroundTruth(x = xIrIreg  , 
    #                                                     F_x = F_x_ReIre , 
    #                                                     f_x = f_x_ReIre ,  
    #                                                     PriorNonImplausibleSet= PriorNonImplausibleSetRegularyIreRegular , 
    #                                                     MD = MD_ReIre , 
    #                                                     Corr_sdhat = Corr_sdhat1 , 
    #                                                     N = n  )
    #ImplausabilityMatrix1 <- FM_CalulateImForGroundTruth(x = xreg  ,
    #                                                     F_x = F_xreg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] ,
    #                                                     f_x = f_xreg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] ,
    #                                                     PriorNonImplausibleSet= PriorNonImplausibleSetRegular[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] ,
    #                                                     MD = MD_Reg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],],
    #                                                     Corr_sdhat =  Corr_sdhat2 ,
    #                                                     N = n  )
    
    #{ 
    #  ImThreshold1 <- apply( ImplausabilityMatrix  , 2 , function(X){ quantile(X , 0.99) } )  # 1.05 for 5% model discrepancy
    #  ImThreshold2 <- apply( ImplausabilityMatrix1  , 2 , function( X ){ quantile(X , 0.99) } ) 
    #}
    
  }
  
}
{
  PriorNonImplausibleSetTotal <- BE_SampleLHSinab( a = apply(rbind(PriorNonImplausibleSetRegular , PriorNonImplausibleSetRegularyIreRegular ) , 2 , min) , b = apply(rbind(PriorNonImplausibleSetRegular , PriorNonImplausibleSetRegularyIreRegular ) , 2 , max) , numbersamples = 1000000 )
  PriorNonImplausibleSetTotal <- PriorNonImplausibleSetTotal[(PriorNonImplausibleSetTotal[,4] < PriorNonImplausibleSetTotal[,5]) ,  ]
  PriorNonImplausibleSetTotal <- PriorNonImplausibleSetTotal[(PriorNonImplausibleSetTotal[,5] < PriorNonImplausibleSetTotal[,6]) ,  ]
  PriorNonImplausibleSetTotal <- PriorNonImplausibleSetTotal[(PriorNonImplausibleSetTotal[,1] < PriorNonImplausibleSetTotal[,2]),  ]
  PriorNonImplausibleSetTotal <- PriorNonImplausibleSetTotal[(PriorNonImplausibleSetTotal[,2] > PriorNonImplausibleSetTotal[,3]),  ]
  PriorNonImplausibleSetTotal[,1:3] <- PriorNonImplausibleSetTotal[ , 1:3] / apply(PriorNonImplausibleSetTotal[ , 1:3] , 1 , sum )
  
  BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[1:500,1:3] , PriorNonImplausibleSetRegular[1:500,1:3] , PriorNonImplausibleSetTotal[1:5000 ,1:3] )  
  x <- seq(0 , 2, 0.01)
  f_total <- t(apply(PriorNonImplausibleSetTotal , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
  F_total <- t(apply(PriorNonImplausibleSetTotal , 1 , function(X){FM_EvalulateCDFEstimate(x , X) }))
  
  xTotal<- t(apply(F_total , 1 , function(X){FM_ExtractActiveOutputs(X , x)} ))
  f_total <- t(apply(cbind(PriorNonImplausibleSetTotal , xTotal) , 1 , function(X){FM_EvaluateDenistyEstimate(  X[11:211] , X[1:10]) }))
  F_total <- t(apply(cbind(PriorNonImplausibleSetTotal , xTotal) , 1 , function(X){FM_EvalulateCDFEstimate(X[11:211] , X[1:10] ) }))
  
  MD_Total <- matrix(0 , dim(F_total)[1] , length(x))
  for( i in 1:dim(F_total)[1] ){
    MD_Total[i , ] <- MDFunction( F_total[i , ] , PriorNonImplausibleSetTotal[i , ] )
  }
  
}


x11()
p1 <- ggplot(data = data.frame(x=xreg[3,] ,  y = f_xreg[3,]), aes(x , y) ) + geom_line(col = 'blue')+
      geom_line(data =data.frame(x=xIrIreg[6,] ,  y = FM_EvaluateDenistyEstimate(xIrIreg[6,]  , c(PriorNonImplausibleSetRegularyIreRegular[6 , 1:9] , 1 ) )) , aes(x , y), col = 'red') + xlab(TeX('r') ) + ylab(TeX('Density') ) + 
      ggtitle('Regular and Irregularly-irregular Heart-rhythm Distributions')
print(p1)