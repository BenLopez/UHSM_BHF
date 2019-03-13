{
  if(file.exists('CheckforDefaultsScript.R')){
    source('CheckforDefaultsScript.R')
  }else{
    pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
    source("LibrariesAndSettings.R" , print.eval  = TRUE )
    DP_LoadPatientIndex()
    DP_ChooseDataReps()
    FilestoProcess <- DP_ChooseECGstoProcess() 
    HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
  }
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

{
{
DataStructure <- array(0 , c(length(listAllPatients) , 8 , 400))

for(i in 1:length(listAllPatients) ){
  PatientID <- listAllPatients[[i]]
  if(DP_CheckFileExists(path , PatientID , Name = paste0(PatientID , '_PwaveSummaries' ))){
    tmp <-  DP_LoadFile(path , PatientID , Name = paste0(PatientID , '_PwaveSummaries' ))
    tmp <- as.matrix( tmp[,tmp[1,]!=0  ] )
    if(ncol(tmp) == 0 ){next}
    if(dim(tmp)[1] !=8){next}
    DataStructure[i , , (dim(DataStructure)[3] - dim(tmp)[2]+1):dim(DataStructure)[3] ] <- tmp
  }
}
}

# Remove rows with all zero rows
PatientsInStudy <- listAllPatients[-which(apply(DataStructure[, 1 , ] , 1 , function(X){sum(X!=0)} ) < 20 )]
DataStructure <- DataStructure[-which(apply(DataStructure[, 1 , ] , 1 , function(X){sum(X!=0)} ) < 20 ) , ,  ]
AFLogical <- apply(as.matrix(PatientsInStudy) , 1 , function(X){DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017  , X))})
alpha = c(sum(AFLogical)/length(AFLogical) , (1-(sum(AFLogical)/length(AFLogical))))
rownames(DataStructure) <- PatientsInStudy
#PwaveIM <- DataStructure[ , 12 , ] 
#DataStructure <- DataStructure[ , 1:11 , ]
#colnames(DataStructure) <- BLBF_CreateNamesPWaveVariables()

#DataStructure[, 1 , ] <- abs(DataStructure[, 1 , ])
#DataStructure[, 10 , ] <- log(DataStructure[, 10 , ])
#DataStructure[, 11 , ] <- log(DataStructure[, 11 , ])
#DataStructure[, 5 , ] <- log(DataStructure[, 5 , ])
#DataStructure[is.infinite(DataStructure)] <- 0
DataStructure[ , 7 , ] <- abs(DataStructure[ , 7 , ] )
DataStructure[ , 8 , ] <- sqrt(DataStructure[ , 8 , ] )
DataStructure[ , 5 , ] <- sqrt(DataStructure[ , 5 , ] )


saveddata <- DataStructure
savedAFLogical <- AFLogical
}


indextoview = 400
#AFLogical <- savedAFLogical[PwaveIM[,indextoview] <2]
#DataStructure <- saveddata[PwaveIM[,indextoview] <2 , ,]

#AFLogical <- savedAFLogical
#DataStructure <- saveddata

BC_PlotCompareTwoHists(DataStructure[AFLogical ==1 ,  , indextoview ] , DataStructure[AFLogical ==0 , , indextoview ] )
sampleofpoints <- sample(which(AFLogical ==0 ) , size = sum(AFLogical ==1) , replace = F)
BC_PlotPairsFromTwoVariables(DataStructure[AFLogical ==1 ,  , indextoview] , DataStructure[sampleofpoints , , indextoview], alpha = 0.2)

DataForLogisticRegression <- data.frame(AFLogical = as.factor(AFLogical) , DataStructure[ ,  , indextoview])

summary(glm(formula = AFLogical ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8,family=binomial(link='logit') , data=DataForLogisticRegression))

PosteriorProbability <- BLBC_FitBayesLinearBayesClassifier(Data = DP_RemoveNaRows(DataStructure[ ,  , indextoview]) , Labels = AFLogical[DP_FindNARows(DataStructure[ , , indextoview])] )


IndividualProbabilityUpdateStruct <- matrix(0 , length(AFLogical) , 30)
for(i in c(0:19)){
  IndividualProbabilityUpdateStruct[,i+1] <- BLBF_CalculatePosteriorProbabilityNBLA(DataStructure = DataStructure , AFLogical = AFLogical   , VariablesToView = c(1:8) ,  indextoview = 381 + i )
}
IndependentProbabilityUpdateStruct <- matrix(0 , length(AFLogical) , 30)
for(i in c(0:29)){
  if(i == 0){
    IndependentProbabilityUpdateStruct[PwaveIM[, 371 + i] > 2,i+1] <- sum(AFLogical)/length(AFLogical)
    IndependentProbabilityUpdateStruct[PwaveIM[, 371 + i] < 2,i+1] <- BLBF_CalculatePosteriorProbabilityNBLA(DataStructure = DataStructure[PwaveIM[, 371 + i] < 2 , , ] , AFLogical = AFLogical[PwaveIM[, 371 + i] < 2]   , VariablesToView = c(1:11) ,  indextoview = 371 + i )
  }
    if(i > 0){
    IndependentProbabilityUpdateStruct[PwaveIM[, 371 + i] > 2,i+1] <- IndependentProbabilityUpdateStruct[PwaveIM[, 371 + i] > 2,i]   
    IndependentProbabilityUpdateStruct[PwaveIM[, 371 + i] < 2,i+1] <- BLBF_CalculatePosteriorProbabilityNBLA(DataStructure = DataStructure[PwaveIM[, 371 + i] < 2 , , ] , AFLogical = AFLogical[PwaveIM[, 371 + i] < 2]  , indextoview = 371 + i , PriorProbability = IndependentProbabilityUpdateStruct[PwaveIM[, 371 + i] < 2,i])}
}

BLBF_PlotProbabilityTrajectories( IndividualProbabilityUpdateStruct  , AFLogical )
title('Individual Probability Trajectories')
BLBF_PlotSenSpecPerformanceSweep( IndividualProbabilityUpdateStruct , AFLogical   )
title('ROC Curves Over Time')

BLBF_PlotProbabilityTrajectories( t(apply(IndividualProbabilityUpdateStruct ,1 ,DP_cummean) ) , AFLogical )
title('Mean Probabilities')
BLBF_PlotSenSpecPerformanceSweep( ProbabilityUpdateStruct = t(apply(IndividualProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical = AFLogical )
title('ROC Curves Over Time')


BLBF_PlotProbabilityTrajectories( IndependentProbabilityUpdateStruct  , AFLogical )
title('Independent Probability Trajectories')
BLBF_PlotSenSpecPerformanceSweep( IndependentProbabilityUpdateStruct , AFLogical   )
title('ROC Curves Over Time')

BLBF_PlotProbabilityTrajectories( t(apply(IndependentProbabilityUpdateStruct ,1 ,DP_cummean) ) , AFLogical )
title('Mean Probabilities')
BLBF_PlotSenSpecPerformanceSweep( ProbabilityUpdateStruct = t(apply(IndependentProbabilityUpdateStruct ,1 ,DP_cummean)) , AFLogical = AFLogical )
title('ROC Curves Over Time')

