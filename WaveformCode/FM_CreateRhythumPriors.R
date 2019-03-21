{
  {
    numberofsamples <- 100000
    PriorNonImplausibleSetRegular <- matrix(0 , numberofsamples , 10)
    
    PriorNonImplausibleSetRegular[, 5] <- runif(numberofsamples , 0.6 , 1)
    PriorNonImplausibleSetRegular[, 4] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] - 0.15*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] - 0.05*PriorNonImplausibleSetRegular[, 5])
    PriorNonImplausibleSetRegular[, 6] <- runif(numberofsamples , PriorNonImplausibleSetRegular[, 5] + 0.05*PriorNonImplausibleSetRegular[, 5] , PriorNonImplausibleSetRegular[, 5] + 0.15*PriorNonImplausibleSetRegular[, 5] )
    
    PriorNonImplausibleSetRegular[, 1] <- runif(numberofsamples , 0 , 0.1)
    PriorNonImplausibleSetRegular[, 3] <- runif(numberofsamples , 0 , 0.1)
    PriorNonImplausibleSetRegular[, 2] <- 1 - (PriorNonImplausibleSetRegular[, 1] + PriorNonImplausibleSetRegular[, 3]) 
    
    PriorNonImplausibleSetRegular[, 7] <- runif(numberofsamples , 0.004 , (0.06*PriorNonImplausibleSetRegular[, 4])  )
    PriorNonImplausibleSetRegular[, 8] <- runif(numberofsamples , 0.001 , (0.06*PriorNonImplausibleSetRegular[, 5]) )
    PriorNonImplausibleSetRegular[, 9] <- runif(numberofsamples ,  apply(as.matrix(PriorNonImplausibleSetRegular[, 7]- 0.15*PriorNonImplausibleSetRegular[, 7]) , 1 , function(X){max(0,X)} ) , PriorNonImplausibleSetRegular[, 7] + 0.15*PriorNonImplausibleSetRegular[, 7] )
    
    PriorNonImplausibleSetRegular[, 10]<- runif(numberofsamples , 0 , 1)
    
    MDFunction <- function(F_x , X){
      output <- sqrt( (F_x*(1-F_x))/481 + 0.000000001 )
      output[output <= 0.0001] <- quantile(output[output > 0] , 0.05, na.rm = T)
      return(output)
    }
    
    #MDFunction <- function(f_x , X){
    #  return( 0.3  )
    #}
    
    x <- seq(0 , 2, 0.01)
    f_xreg <- t(apply(PriorNonImplausibleSetRegular , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
    F_xreg <- t(apply(f_xreg , 1 , function(X){cumsum( c(0,diff(x))*(X/sum(c(0,diff(x))*(X))) ) }))
    
    MD_Reg <- matrix(0 , numberofsamples , length(x))
    for( i in 1:numberofsamples){
      MD_Reg[i , ] <- MDFunction( F_xreg[i , ] , PriorNonImplausibleSetRegular[i , ] )
    }
  }
  
  # Prior specification regularly-irregular
  {
    PriorNonImplausibleSetRegularyIreRegular <- matrix(0 , numberofsamples , 10)
    
    PriorNonImplausibleSetRegularyIreRegular[, 5] <- runif(numberofsamples , 0.3 , 1)
    PriorNonImplausibleSetRegularyIreRegular[, 4] <- runif(numberofsamples , 0.25 ,  apply(as.matrix(0.9*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 , function(X){max(X , 0.26)}) )
    PriorNonImplausibleSetRegularyIreRegular[, 6] <- runif(numberofsamples , apply(as.matrix(PriorNonImplausibleSetRegularyIreRegular[, 5] + 0.1*PriorNonImplausibleSetRegularyIreRegular[, 5]) , 1 ,function(X){max(0.1 , X)} ) , 2)
    
    PriorNonImplausibleSetRegularyIreRegular[, 2] <- runif(numberofsamples , 0.4 , 0.8)
    PriorNonImplausibleSetRegularyIreRegular[, 1] <- runif(numberofsamples , rep(0 , numberofsamples ) ,  apply( cbind(PriorNonImplausibleSetRegularyIreRegular[, 2] , (1 - PriorNonImplausibleSetRegularyIreRegular[,2])) , 1 , min) )
    PriorNonImplausibleSetRegularyIreRegular[, 3] <- 1 - PriorNonImplausibleSetRegularyIreRegular[, 1] - PriorNonImplausibleSetRegularyIreRegular[, 2]
    
    PriorNonImplausibleSetRegularyIreRegular[, 8] <- runif(numberofsamples , 0.05 , 0.2 )
    PriorNonImplausibleSetRegularyIreRegular[, 7] <- runif(numberofsamples , 0.8*PriorNonImplausibleSetRegularyIreRegular[, 8] , PriorNonImplausibleSetRegularyIreRegular[, 8] )
    PriorNonImplausibleSetRegularyIreRegular[, 9] <- runif(numberofsamples , 0.8*PriorNonImplausibleSetRegularyIreRegular[, 8] , PriorNonImplausibleSetRegularyIreRegular[, 8] )
    PriorNonImplausibleSetRegularyIreRegular[, 10] <- 1
    
    f_x_ReIre <- t(apply(PriorNonImplausibleSetRegularyIreRegular , 1 , function(X){FM_EvaluateDenistyEstimate(x , X) }))
    
    F_x_ReIre <- t(apply(f_x_ReIre , 1 , function(X){cumsum( c(0,diff(x))*(X/sum(c(0,diff(x))*(X))) ) }))
    
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
    
    MDFunction2 <- function(F_x , X){
      return(MDFunction(F_x , X))
    }
    
    MD_ReIre <- matrix(0 , dim(PriorNonImplausibleSetRegularyIreRegular)[1] , length(x))
    for( i in 1:dim(PriorNonImplausibleSetRegularyIreRegular)[1]){
      MD_ReIre[i , ] <- MDFunction2(F_x_ReIre[i , ] , PriorNonImplausibleSetRegularyIreRegular[i ,])
    }
  }
  
  x <- x + median(diff(x)/2)
  
  # Calculate approximate correlation matrix.
  { SampleMatrix  <-  FM_SampleRealisationsSet( PriorNonImplausibleSetRegularyIreRegular , N = (500 - 19) , x  )
    ME_matrix     <-  SampleMatrix - F_x_ReIre
    Cov_sdhat1    <-  cov( ME_matrix )
    Corr_sdhat1   <-  cov2cor( DP_AddNugget(Cov_sdhat1 , 0.0045) )
    
    
    SampleMatrix <- FM_SampleRealisationsSet( PriorNonImplausibleSet =  PriorNonImplausibleSetRegular[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] , N = (500 - 19) , x =  x   )
    ME_matrix <- SampleMatrix - F_xreg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],]
    Cov_sdhat2 <- cov( ME_matrix )
    Corr_sdhat2 <- cov2cor( DP_AddNugget(Cov_sdhat2 , 0.0045) )
  }
  
  { 
    ImplausabilityMatrix <- FM_CalulateImForGroundTruth(x = x  , F_x = F_x_ReIre , f_x = f_x_ReIre ,  PriorNonImplausibleSet= PriorNonImplausibleSetRegularyIreRegular , MD = MD_ReIre , Corr_sdhat1 , N = (500 - 19)  )
    ImplausabilityMatrix1 <- FM_CalulateImForGroundTruth(x = x  , F_x = F_xreg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] , f_x = f_xreg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] , PriorNonImplausibleSet= PriorNonImplausibleSetRegular[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],] , MD = MD_Reg[1:dim(PriorNonImplausibleSetRegularyIreRegular)[1],], Corr_sdhat2 , N = (500 - 19)  )
    
    { 
      ImThreshold1 <- 1.05*apply( ImplausabilityMatrix  , 2 , function(X){ quantile(X , 0.99) } )  # 1.05 for 5% model discrepancy
      ImThreshold2 <- 1.05*apply( ImplausabilityMatrix1  , 2 , function( X ){ quantile(X , 0.99) } ) 
    }
    
  }
  
}
