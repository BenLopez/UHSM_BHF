
numberofsimulations = 1100
numberofrepeats = 5000
precomputedfolderpath <- DP_SelectPrecomputedFolder()

if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath,'EmulatingCriticalThresholdDataCDF.RData') ){
  load(paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataCDF.RData'))
  Im_Crit <- EmulatingCriticalThresholdData[[1]]
  X <- EmulatingCriticalThresholdData[[2]]
}else{
  xx <- seq(0,5,0.01)
  outputstruct <- CTEm_CreateTrainingSetCDF(xx,xTotal,PriorNonImplausibleSetTotal,F_total,f_total,MD_Total  ,numberofsimulations = numberofsimulations, c = c , l = l , N = n ,numberofrepeats = numberofrepeats)
  EmulatingCriticalThresholdData <- list( Im_crit = outputstruct , X = PriorNonImplausibleSetTotal[1:numberofsimulations , ])
  save( EmulatingCriticalThresholdData , file = paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataCDF.RData') )
}




