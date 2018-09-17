ECGs <- DP_LoadReducedECGs(path , Patinettotest , numberrep = numberrep , FilestoProcess = FilestoProcess)
RPeaksStruct <- DP_LoadRpeaksfile(path , Patinettotest)
if( BCOptions$DataType == 'DistributionSummaries' ){
  DistributionSummaries <- DP_LoadDistributionSummaries(path , Patinettotest)
}
if( BCOptions$DataType == 'CDFs' ){
  CDFs <-  DP_LoadFile(path , PatientsID = Patinettotest , Name = paste0(Patinettotest , '_CDFs' )) 
}

MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = Patinettotest)


if( BCOptions$DataType == 'DistributionSummaries' ){
  if( BCOptions$DensityEstimationlocal == 'GMM' ){
    AdjustedBeliefs <- BC_EstimateGlobalandLocalParametersDisSum(DistributionSummaries)
    Implausability <- BC_CalulateCulmulativeImplausability( Implausability = BC_TestingCalulateRegularisedImplausabilty(AdjustedBeliefs$W , SecondOrderStruct = LocalSecondOrderStruct , ImSecondOrderStruct) , n = BCParameters$TS_Likelihood_clique )
    GlobalImplausabilities <- BC_GlobalTestingCalulateRegularisedImplausabilty( AdjustedBeliefs$M , GlobalSecondOrderStruct , GlobalImSecondOrderStruct)
    
    if(BCOptions$GlobalUpdate == 'Yes'){
      GlobalUpdatedProbabilities <- BC_GlobalBayesianBeliefUpdateMVN(M = AdjustedBeliefs$M , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = Priorprobabilities)
    }else{
      GlobalUpdatedProbabilities <- Priorprobabilities
    }
    
    PosteriorProbabilities <- BC_BayesianBeliefUpdateGMM( W=AdjustedBeliefs$W , LocalDistributionStruct = LocalDistributionStruct , Probabilities = GlobalUpdatedProbabilities , n = BCParameters$TS_Likelihood_clique)
    #PosteriorProbabilities[  Implausability[,2]< 0.3 , 1] <- 0
    #PosteriorProbabilities[  Implausability[,2]< 0.3 , 2] <- 1
    PosteriorProbabilities[ (Implausability[,1] > 3)*(Implausability[,2] > 3) ==1 , 1] <- GlobalUpdatedProbabilities$B
    PosteriorProbabilities[ (Implausability[,1] > 3)*(Implausability[,2] > 3) ==1 , 2] <- GlobalUpdatedProbabilities$`B^c`
    PosteriorProbabilities[ Implausability[,1] > 3 , 1] <- 0
    PosteriorProbabilities[ Implausability[,1] > 3 , 2] <- 1
  }
}


if( BCOptions$DataType == 'CDFs' ){
  if( BCOptions$DensityEstimationlocal == 'GMM' ){
    AdjustedBeliefs <- BC_EstimateGlobalandLocalParametersCDFs(CDFs )
    Implausability  <- BC_CalulateCulmulativeImplausability( Implausability = BC_TestingCalulateRegularisedImplausabilty(TestData = AdjustedBeliefs$W , SecondOrderStruct = LocalSecondOrderStruct , ImSecondOrderStruct = ImSecondOrderStruct) , n = BCParameters$TS_Likelihood_clique )
    GlobalImplausabilities <- BC_GlobalTestingCalulateRegularisedImplausabilty( AdjustedBeliefs$M , GlobalSecondOrderStruct , GlobalImSecondOrderStruct)
    
    if(BCOptions$GlobalUpdate == 'Yes'){
      GlobalUpdatedProbabilities <- BC_GlobalBayesianBeliefUpdateMVN(M = AdjustedBeliefs$M , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = Priorprobabilities)
    }else{
      GlobalUpdatedProbabilities <- Priorprobabilities
    }
    
    PosteriorProbabilities <- BC_BayesianBeliefUpdateGMM( W=AdjustedBeliefs$W , LocalDistributionStruct = LocalDistributionStruct , Probabilities = GlobalUpdatedProbabilities , n = BCParameters$TS_Likelihood_clique)
    
    #mIm <- mean(Implausability[1000:5000 , 2])
    #VIm <- var(Implausability[1000:5000 , 2])
    #ImThresh <- mIm + (3*sqrt(VIm))
    
    #PosteriorProbabilities[  Implausability[,2]< ImThresh , 1] <- 0
    #PosteriorProbabilities[  Implausability[,2]< ImThresh , 2] <- 1
    PosteriorProbabilities[ (Implausability[,1] > 3)*(Implausability[,2] > 3) ==1 , 1] <- GlobalUpdatedProbabilities$B
    PosteriorProbabilities[ (Implausability[,1] > 3)*(Implausability[,2] > 3) ==1 , 2] <- GlobalUpdatedProbabilities$`B^c`
    PosteriorProbabilities[ Implausability[,1] > 3 , 1] <- 0
    PosteriorProbabilities[ Implausability[,1] > 3 , 2] <- 1
  }
}


AnnotatedAFMetaData <- BC_CreateAFAnnotationFomMetaData(t = RPeaksStruct$RRCombined$t , MetaData = MetaData)

logicaltimeseries <-  ( PosteriorProbabilities[ , 1] > BCParameters[[3]] )
logicaltimeseries[is.na(logicaltimeseries)] <- FALSE

#iqr <- rollapply(RPeaksStruct$RRCombined$RR , width = 250 , na.pad = TRUE , FUN = function(X){IQR(X, na.rm= TRUE)} )
#iqr[is.na(iqr)] <- 0

AFLocations <- ASWF_GetStartEndAF(t = RPeaksStruct$RRCombined$t , logicaltimeseries = logicaltimeseries , minutethreshold = BCParameters[[4]])
logicaltimeseries <-(Implausability[,1] > 3)*(Implausability[,2] > 3) ==1
logicaltimeseries[is.na(logicaltimeseries)] <- TRUE
BadDataLocations <- ASWF_GetStartEndAF(RPeaksStruct$RRCombined$t , logicaltimeseries = logicaltimeseries , minutethreshold = 0.1)
if(nrow(AFLocations) > 0){
AFLocations <- BC_CheckForTimeGaps(AFLocations = AFLocations , BadDataLocations =BadDataLocations , t = RPeaksStruct$RRCombined$t , ECGs = ECGs  )
}
AnnotatedAFInference <- BC_CreateAnnotationFromInference(t =RPeaksStruct$RRCombined$t , AFLocations = AFLocations)

AnnotatedBadLocations <- BC_CreateAnnotationFromInference(t = RPeaksStruct$RRCombined$t , AFLocations = BadDataLocations)
Performance <- BC_CalulatePerformance(AnnotatedAFMetaData ,  AnnotatedAFInference , AnnotatedBadLocations )
