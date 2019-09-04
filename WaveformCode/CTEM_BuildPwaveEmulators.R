numberofsimulations = 1100
numbersamples = 500000

CTEM_CreateTrainingsetPwaves <- function(E_Beta,V_Beta,numberofsimulations = 1100,numberofrepitions = 500000, clength = 0.1 , q = 0.99){
  PwaveImplausibilityThresholds <- matrix(0 , numberofsimulations , 2)
  PwaveImplausibilityEDF <- array(0 , c(numberofsimulations , 2,201) )
  
  for(i in 1:numberofsimulations){
    outputStructure <- FM_CalculateEDFPwaves(Xstar , PriorNonImplausibleSet[i,] , E_Beta , V_Beta, numbersamples = numberofrepitions, clength = clength , q = q )
    PwaveImplausibilityThresholds[i,1] <- outputStructure[[1]]
    PwaveImplausibilityThresholds[i,2] <- outputStructure[[2]]
    PwaveImplausibilityEDF[i,1,] <- outputStructure[[3]]
    PwaveImplausibilityEDF[i,2,] <- outputStructure[[4]]
    DP_WaitBar(i/1100)
  }
  return( list(PwaveImplausibilityThresholds , PwaveImplausibilityEDF , PriorNonImplausibleSet[1:numberofsimulations,]) )
}

if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath,'EmulatingCriticalThresholdData.RData') ){
  load(paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdData.RData'))
  Im_Crit <- EmulatingCriticalThresholdData[[1]]
  X <- EmulatingCriticalThresholdData[[2]]
}else{
  outputstruct <- CTEm_CreateTrainingSet( x = xTotal,
                                          PriorNonImplausibleSetTotal = PriorNonImplausibleSetTotal ,
                                          F_total = F_total ,
                                          f_total = f_total ,
                                          MD_Total = MD_Total ,
                                          numberofsimulations = numberofsimulations ,
                                          c = c , 
                                          l = l , 
                                          N = n ,
                                          numberofrepeats = numberofrepeats )
  EmulatingCriticalThresholdData <- list( Im_crit = outputstruct , X = PriorNonImplausibleSetTotal[1:numberofsimulations , ])
  Im_Crit <- EmulatingCriticalThresholdData[[1]]
  X <- EmulatingCriticalThresholdData[[2]]
  save( EmulatingCriticalThresholdData , file = paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdData.RData') )
}



