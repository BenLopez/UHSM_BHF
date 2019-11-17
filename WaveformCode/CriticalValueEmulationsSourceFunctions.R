CTEm_CalculateImThreshold <- function(x, X , F_X, MD_X ,f_X, N=500 , numberofrepeats = 100 , c = 5000 , l = 0.005 , p = 0.99){
  
  ImplausabilityMatrix <- matrix(0 , numberofrepeats , 4)
  G_0 <- function(N){ FM_SampleGMM( X = X , N) }
  
  for(ii in 1:numberofrepeats){
    RRtimes <- FM_SampleDP(c , l , N , G_0)
    y <- FM_CalculateCDFS( RRtimes = RRtimes , xx = x )
    
    NonZeroLogical <- matrix(T , length(RRtimes) , 1)
    
    #RPeakKDEEstmate <- kde( RRtimes )
    #z <- predict(RPeakKDEEstmate , x = x )
    
    ImplausabilityMatrix[ii , 1] <- mean( abs(y[NonZeroLogical] - F_X[NonZeroLogical]) / (MD_X[NonZeroLogical] + 0.0045) , na.rm = T)
    ImplausabilityMatrix[ii , 2] <- max( abs(y[NonZeroLogical] - F_X[NonZeroLogical]) / (MD_X[NonZeroLogical]+ 0.0045) , na.rm = T)
    #ImplausabilityMatrix[ii , 3] <- (( (z[NonZeroLogical] - f_X[NonZeroLogical])/f_X[NonZeroLogical] ) )[which.max(z[NonZeroLogical])]
    ImplausabilityMatrix[ii , 3] <- 0
    #ImplausabilityMatrix[ii , 4] <- sum(log(FM_EvaluateDenistyEstimate(RRtimes ,X )))
    ImplausabilityMatrix[ii , 4] <- 0
  }
  
  return(apply(ImplausabilityMatrix , 2 , function(X){quantile(X , p)}) )
  
}
CTEm_CreateTrainingSet <- function(x,PriorNonImplausibleSetTotal,F_total,f_total,MD_Total  ,numberofsimulations = 1100, c = 5000 , l = 0.005 , N = 500-19 ,numberofrepeats = 50000){
  outputstruct <- matrix(0 , numberofsimulations , 4)
  
  for(ii in 1:numberofsimulations){
    X <- PriorNonImplausibleSetTotal[ ii , ]
    F_X = F_total[ ii , ]
    MD_X = MD_Total[ ii , ]
    f_X = f_total[ ii , ]
    xx = as.matrix(x[ii,])
    outputstruct[ii,] = CTEm_CalculateImThreshold( x = xx ,
                                                   X =  X ,
                                                   F_X =  F_X , 
                                                   MD_X = MD_X ,
                                                   f_X =  f_X ,
                                                   N = N ,
                                                   numberofrepeats = numberofrepeats)
    DP_WaitBar(ii/numberofsimulations)
  }  
  
  return(outputstruct)
}
CTEm_CreateTrainingSetCDF <- function(xx,x,PriorNonImplausibleSetTotal,F_total,f_total,MD_Total  ,numberofsimulations = 1100, c = 5000 , l = 0.005 , N = 500-19 ,numberofrepeats = 50000){
  outputstruct <- array(0 , c(numberofsimulations , length(xx) , 2) )
  
  for(ii in 1:numberofsimulations){
    X <- PriorNonImplausibleSetTotal[ ii , ]
    F_X = F_total[ ii , ]
    MD_X = MD_Total[ ii , ]
    f_X = f_total[ ii , ]
    xxx = as.matrix(x[ii,])
    outputstruct[ii,,] = CTEm_CalculateImThresholdCDF(xx,
                                                      x = xxx,
                                                      X ,
                                                      F_X,
                                                      MD_X ,
                                                      f_X, 
                                                      N = N , numberofrepeats = numberofrepeats)
    DP_WaitBar(ii/numberofsimulations)
  }  
  
  return(outputstruct)
}
CTEm_CalculateImThresholdCDF <- function(xx,x,X , F_X, MD_X ,f_X, N=500 , numberofrepeats = 100 , c = 5000 , l = 0.005 , p = 0.99){
  
  ImplausabilityMatrix <- matrix(0 , numberofrepeats , 4)
  G_0 <- function(N){ FM_SampleGMM( X = X , N) }
  
  for(ii in 1:numberofrepeats){
    RRtimes <- FM_SampleDP(c , l , N , G_0)
    y <- FM_CalculateCDFS( RRtimes = RRtimes , xx = x )
    
    NonZeroLogical <- matrix(T , length(RRtimes) , 1)
    
    #RPeakKDEEstmate <- kde( RRtimes )
    #z <- predict(RPeakKDEEstmate , x = x )
    
    ImplausabilityMatrix[ii , 1] <- mean( abs(y[NonZeroLogical] - F_X[NonZeroLogical]) / (MD_X[NonZeroLogical] + 0.0045) , na.rm = T)
    ImplausabilityMatrix[ii , 2] <- max( abs(y[NonZeroLogical] - F_X[NonZeroLogical]) / (MD_X[NonZeroLogical]+ 0.0045) , na.rm = T)
    #ImplausabilityMatrix[ii , 3] <- (( (z[NonZeroLogical] - f_X[NonZeroLogical])/f_X[NonZeroLogical] ) )[which.max(z[NonZeroLogical])]
    ImplausabilityMatrix[ii , 3] <- 0
    #ImplausabilityMatrix[ii , 4] <- sum(log(FM_EvaluateDenistyEstimate(RRtimes ,X )))
    ImplausabilityMatrix[ii , 4] <- 0
  }
  
  return(cbind(FM_CalculateCDFS(ImplausabilityMatrix[,1] , xx ) , FM_CalculateCDFS(ImplausabilityMatrix[,2] , xx ) ) )
  
}
CTEM_CreateTrainingsetPwaves <- function(XPwave,PriorNonImplausibleSet,E_Beta,V_Beta,numbertrainingpoints = 1100,numberofrepitions = 500000,clength = 0.1 , q = 0.99 ){
  PwaveImplausibilityThresholds <- matrix(0 , numbertrainingpoints , 2)
  PwaveImplausibilityEDF <- array(0 , c(numbertrainingpoints , 2,length(seq(0,2,0.005))) )
  
  for(i in 1:numbertrainingpoints){
    #for(i in 1:dim( PriorNonImplausibleSet)[1]){
    if(i == 1){
      i <- dim(XPwave)[1]
      outputStructure <- FM_CalculateEDFPwaves(Xstar = XPwave[i,] ,
                                               X=PriorNonImplausibleSet[i,] , 
                                               E_Beta , 
                                               V_Beta, 
                                               numbersamples = numberofrepitions, 
                                               clength =clength , 
                                               q =q )
    i <- 1
}else{
    outputStructure <- FM_CalculateEDFPwaves(Xstar = XPwave[i,] ,
                                             X=PriorNonImplausibleSet[i,] , 
                                             E_Beta , 
                                             V_Beta, 
                                             numbersamples = numberofrepitions, 
                                             clength =clength , 
                                             q =q )
}
        PwaveImplausibilityThresholds[i,1] <- outputStructure[[1]]
    PwaveImplausibilityThresholds[i,2] <- outputStructure[[2]]
    PwaveImplausibilityEDF[i,1,] <- outputStructure[[3]]
    PwaveImplausibilityEDF[i,2,] <- outputStructure[[4]]
    #DP_WaitBar(i/dim( PriorNonImplausibleSet)[1])
    DP_WaitBar(i/numbertrainingpoints)
  }
  return(setNames(list(PwaveImplausibilityThresholds , PwaveImplausibilityEDF,PriorNonImplausibleSet[1:numbertrainingpoints,] ) , c('Thresholds' , 'EDF','Points') ) )
}