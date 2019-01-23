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
    PosteriorProbabilities[ ((Implausability[,1] > 2)*(Implausability[,2] > 2)) ==1 , 1] <- GlobalUpdatedProbabilities$B
    PosteriorProbabilities[ ((Implausability[,1] > 2)*(Implausability[,2] > 2)) ==1 , 2] <- GlobalUpdatedProbabilities$`B^c`
    PosteriorProbabilities[ Implausability[,1] > 2 , 1] <- 0
    PosteriorProbabilities[ Implausability[,1] > 2 , 2] <- 1
    PosteriorProbabilities[ ((Implausability[,2] > 2)*((Implausability[,2]/Implausability[,1]) > 3)) == 1 , 1] <- 1
    PosteriorProbabilities[ ((Implausability[,2] > 2)*((Implausability[,2]/Implausability[,1]) > 3)) == 1 , 2] <- 0
    
    
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

if(DP_CheckIfAFPatient(MetaData = MetaData)){
TrueAFLocations <- AFD_GetStartEndAF(t = RPeaksStruct$RRCombined$t , logicaltimeseries = AnnotatedAFMetaData , minutethreshold = BCParameters[[4]])
TrueAFLocations <- BC_CheckForTimeGaps(AFLocations = TrueAFLocations , BadDataLocations = BadDataLocations , t = RPeaksStruct$RRCombined$t , ECGs = ECGs  )
AnnotatedAFMetaData <- BC_CreateAnnotationFromInference(t = RPeaksStruct$RRCombined$t , TrueAFLocations)
}

logicaltimeseries <-  ( PosteriorProbabilities[ , 1] > BCParameters[[3]] )
#logicaltimeseries <-  ( PosteriorProbabilities[ , 1] > 0.2 )
logicaltimeseries[is.na(logicaltimeseries)] <- FALSE

# Testing component
logicaltimeseries <- rollmean(logicaltimeseries , k = 500, align = c( "right") , na.pad = TRUE) > 0.5
logicaltimeseries[is.na(logicaltimeseries)] <- FALSE
AFLocations <- AFD_GetStartEndAF(t = RPeaksStruct$RRCombined$t , logicaltimeseries = logicaltimeseries , minutethreshold = BCParameters[[4]])


logicaltimeseries <- (Implausability[,1] > 3)*(Implausability[,2] > 3) ==1
logicaltimeseries[is.na(logicaltimeseries)] <- TRUE
BadDataLocations <- AFD_GetStartEndAF(RPeaksStruct$RRCombined$t , logicaltimeseries = logicaltimeseries , minutethreshold = 0.1)

if(nrow(AFLocations) > 0){
AFLocations <- BC_CheckForTimeGaps(AFLocations = AFLocations , BadDataLocations =BadDataLocations , t = RPeaksStruct$RRCombined$t , ECGs = ECGs  )
}

AnnotatedAFInference <- BC_CreateAnnotationFromInference(t =RPeaksStruct$RRCombined$t , AFLocations = AFLocations)


if(BCOptions$PAmplitudeAnalysis == 'Yes'  ){
 if( nrow(AFLocations) > 0 ){
   AFbyPwavesLogical <- matrix(0 , nrow(AFLocations) , 1)
   for(i in 1:nrow(AFLocations) ){
   PAmplitudes <- BC_CalculatePAmplitudes(AFLocations = AFLocations[i,] ,RPeaksStruct =  RPeaksStruct ,ECGs =  ECGs , AnnotatedAFInference = AnnotatedAFInference)
   if(is.na(PAmplitudes$AFPAmplitude)){
     AFbyPwavesLogical[i,] <- 1
     next
   }
   if(  abs(PAmplitudes$NAFPAmplitude - PAmplitudes$AFPAmplitude)/abs(PAmplitudes$NAFPAmplitude) > 0.45 ){
     AFbyPwavesLogical[i,] <- 1  
   if(PAmplitudes$AFPAmplitude > 0 && PAmplitudes$NAFPAmplitude > 0 && (PAmplitudes$AFPAmplitude > PAmplitudes$NAFPAmplitude) ){
     AFbyPwavesLogical[i,] <- 0  
   }
   if(PAmplitudes$AFPAmplitude < 0 && PAmplitudes$NAFPAmplitude < 0 && (abs(PAmplitudes$AFPAmplitude) > abs(PAmplitudes$NAFPAmplitude) ) ){
     AFbyPwavesLogical[i,] <- 0  
   }
   }
   if( abs(PAmplitudes$NAFPAmplitude) < 2 && abs(PAmplitudes$AFPAmplitude) < 6 ){ AFbyPwavesLogical[i,] <- 1 }
   }
   AFLocations <- AFLocations[which(AFbyPwavesLogical == 1) , ]   
   }
}
AFLocations <- BC_CleanAFTimes(AFLocations , minutes = 8)
BadDataLocations <- BC_CleanAFTimes(DP_SortMatrix(BadDataLocations) , minutes = 20)

AnnotatedAFInference <- BC_CreateAnnotationFromInference(t =RPeaksStruct$RRCombined$t , AFLocations = AFLocations)
AnnotatedBadLocations <- BC_CreateAnnotationFromInference(t = RPeaksStruct$RRCombined$t , AFLocations = BadDataLocations)
Performance <- BC_CalulatePerformance(AnnotatedAFMetaData = AnnotatedAFMetaData , AnnotatedAFInference =  AnnotatedAFInference , BadDataLocations = AnnotatedBadLocations )
