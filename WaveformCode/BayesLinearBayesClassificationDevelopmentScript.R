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

BLBC_AFibDetExtractData <- function(listAllPatients  , numberofpoints = 100){
 
ReformulateMatrix <- function(X , numberofpoints){
    output <- matrix(0 , floor(dim(X)[1]/numberofpoints) , 11*numberofpoints)
    for(i in 1:dim(output)[1]){
      output[i , ] <- as.matrix(X[(((i-1)*numberofpoints) +1) :(i*numberofpoints) , 1:11])
    }
    return(output)  
  }  
  
  StartIndex <- 1
  DataArray <- array(NA , c(21000000 / numberofpoints , (11*numberofpoints)) )
  AFLogicalArray <- array(NA , c(21000000 / numberofpoints , 1) )
  InitialArray <- array(NA , c(21000000 / numberofpoints , (11)) )
  
  #length(listAllPatients)
  for(ii in 1:length(listAllPatients) ){
    patientID <- listAllPatients[[ii]]
    MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , patientID)
    if(DP_CheckDistributionSummariesExists(path , listAllPatients[[ii]]) == FALSE ){
      next}else{
        DistributionSummaries <- data.frame(DP_LoadDistributionSummaries(path , listAllPatients[[ii]]))
        DistributionSummaries$Var.Modes[DistributionSummaries$Num.Modes == 1] =0
        DistributionSummaries$var.Mode.Density[DistributionSummaries$Num.Modes == 1] = 0
      }
    
    if( DP_CheckIfAFPatient(MetaData)  ){
      AFAnnotation <- BC_CreateAFAnnotationFomMetaData(DistributionSummaries$time , MetaData = MetaData) ==1
      AFDistributionSummaries <-  ReformulateMatrix(DistributionSummaries[AFAnnotation , ] , numberofpoints)
      NAFDistributionSummaries <- ReformulateMatrix(DistributionSummaries[AFAnnotation == F , ] , numberofpoints)
      
      DistributionSummaries <- rbind(AFDistributionSummaries , NAFDistributionSummaries)
      AFAnnotation <- rbind(as.matrix(rep(T , dim(AFDistributionSummaries)[1]) ) , as.matrix( rep(F , dim(NAFDistributionSummaries)[1])) )
      
      DataArray[StartIndex:(StartIndex + dim(DistributionSummaries)[1] - 1) ,  ] <- DistributionSummaries
      AFLogicalArray[StartIndex:(StartIndex + dim(DistributionSummaries)[1] - 1)] <- AFAnnotation
      InitialArray[StartIndex:(StartIndex + dim(DistributionSummaries)[1] - 1) ,  ] <- apply(NAFDistributionSummaries[500:1000 , ] , 2 , function(X){mean(X , na.rm = T)} )
      StartIndex <- StartIndex + dim(DistributionSummaries)[1] + 1
      
    }else{
      NAFAnnotation <- rep(FALSE , dim(DistributionSummaries)[1] )
      NAFDistributionSummaries <- ReformulateMatrix(DistributionSummaries[NAFAnnotation == F , ] , numberofpoints)
      NAFAnnotation <- rep(FALSE , dim(DistributionSummaries)[1] )
      
      DataArray[StartIndex:(StartIndex + dim(DistributionSummaries)[1] - 1) ,  ] <- NAFDistributionSummaries
      AFLogicalArray[StartIndex:(StartIndex + dim(DistributionSummaries)[1] - 1)] <- NAFAnnotation
      InitialArray[StartIndex:(StartIndex + dim(DistributionSummaries)[1] - 1) ,  ] <- apply(NAFDistributionSummaries[500:1000 , ] , 2 , function(X){mean(X , na.rm = T)} )
      StartIndex <- StartIndex + dim(DistributionSummaries)[1] + 1
    }
    DP_WaitBar(ii/length(listAllPatients) )
  }
  
  AFLogicalArray <- AFLogicalArray[DP_FindNARows(DataArray)  ]
  InitialArray <- InitialArray[DP_FindNARows(DataArray) , ]
  DataArray <- DataArray[DP_FindNARows(DataArray) , ]
  return(setNames(list(DataArray , AFLogicalArray , InitialArray) , c('Data' , 'AFLogical' , 'InitialArray')) )
  
}


ExtractedData <- BLBC_AFibDetExtractData(listAllPatients , numberofpoints = 1)

mAF <- apply(ExtractedData$Data[ExtractedData$AFLogical , ] , 2 , mean)
cAF <- cov(ExtractedData$Data[ExtractedData$AFLogical , ])
cAF <- solve( DP_AddNugget(cAF , 1e-06*diag(diag(cAF) )))

mNAF <- apply(ExtractedData$Data[ExtractedData$AFLogical ==F, ] , 2 , mean)
cNAF <- cov(ExtractedData$Data[ExtractedData$AFLogical ==F, ])
cNAF <- solve(DP_AddNugget(cNAF , 1e-06*diag(diag(cNAF) )))

ImAF <- apply(ExtractedData$Data , 1 , function(X){log((X - mAF)%*%(cAF)%*%(X - mAF))} )
ImNAF <- apply(ExtractedData$Data , 1 , function(X){log((X - mNAF)%*%(cNAF)%*%(X - mNAF))} )

ImAF <- (ImAF - mean(ImAF[ExtractedData$AFLogical]))/sqrt(var(ImAF[ExtractedData$AFLogical]))
ImNAF <- (ImNAF - mean(ImNAF[ExtractedData$AFLogical == F]))/sqrt(var(ImNAF[ExtractedData$AFLogical == F]))

x11()
par(mfrow = c(1 , 2))
BC_PlotCompareSingleHists(ImAF[ExtractedData$AFLogical == 0] , ImAF[ExtractedData$AFLogical] , main = 'ImAF')
BC_PlotCompareSingleHists(ImNAF[ExtractedData$AFLogical == 0] , ImNAF[ExtractedData$AFLogical], main = 'NAF')

ImMatrix <- cbind(ImAF , ImNAF)

SampleofPoints1 <- sample(which(ExtractedData$AFLogical) , 10000)
SampleofPoints2 <- sample(which(ExtractedData$AFLogical ==0) , 10000)
  
x11()
plot(ImMatrix[SampleofPoints1 , 1] , ImMatrix[SampleofPoints1 , 2] , col = rgb(1 , 0 , 0 , alph = 0.01) , xlim = c(-3,8) , ylim = c(-4 , 5) , pch = 16 , xlab = 'Implausibility AFib' , ylab = 'Implausibility Not-AFib')
points(ImMatrix[SampleofPoints2 , 1] , ImMatrix[SampleofPoints2 ,2] , col = rgb(0 , 0 , 1 , alph = 0.01) , pch = 16)
#abline(-3 , 0)
abline(3 , 0)
#abline( v = -3 )
abline( v = 3 )
title('Implausibility Pairs')


SampleofPoints1 <- sample(which( ((ExtractedData$AFLogical)*(ImAF <3)) == 1 ) , 100000 , replace = F)
AFModel <- kde( ImMatrix[SampleofPoints1 , ]  )

SampleofPoints1 <- sample(which( ((ExtractedData$AFLogical ==F)*(ImAF <3)) == 1 ) , 100000, replace = F)
NAFModel <- kde( ImMatrix[SampleofPoints1 , ]  )

plot(NAFModel , col = 'blue' , add = T)
plot(AFModel , col = 'red' , add = T)

F_1 <- BLBF_CalculateDensities(AFModel = AFModel , NAFModel = NAFModel , ImMatrix = ImMatrix[ImAF <3 , ] ,AFLogical = ExtractedData$AFLogical )

#alpha <- c(sum(ExtractedData$AFLogical)/length(ExtractedData$AFLogical) ,  1- sum(ExtractedData$AFLogical)/length(ExtractedData$AFLogical))

alpha <- c(sum(ExtractedData$AFLogical[ImAF <3])/length(ExtractedData$AFLogical[ImAF <3]) ,  1- sum(ExtractedData$AFLogical[ImAF <3])/length(ExtractedData$AFLogical[ImAF <3]))

PosteriorProbabilities <- (alpha[1]*F_1[,1]) / ((alpha[1]*F_1[,1]) + (alpha[2]*F_1[,2]))

PosteriorProbabilities[is.na(PosteriorProbabilities)] <- alpha[1] 

rm(ImMatrix , ImAF , PosteriorProbabilities , F_1 , ImNAF , SampleofPoints1 ,SampleofPoints2)

x11()
EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbabilities , ExtractedData$AFLogical[ImAF < 3]), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities'))

x11()
PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbabilities , ExtractedData$AFLogical[ImAF < 3]))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('Forecasting ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Forecasting PPV vs NPV Curves')
grid.arrange(ROCplot , NPVPPVPlot , nrow= 2)


AFSOS <- BLBF_CalculateSecondOrderStruct( cbind(ExtractedData$InitialArray[ExtractedData$AFLogical , ] , ExtractedData$Data[ExtractedData$AFLogical , ] ) )
SampleofPoints1 <- sample(which(ExtractedData$AFLogical==0) , 1000000 )
NAFSOS <- BLBF_CalculateSecondOrderStruct( cbind(ExtractedData$InitialArray[SampleofPoints1, ] , ExtractedData$Data[SampleofPoints1 , ] ) )


ImMatrix <- matrix(0 , a , 2)
for(i in 1:a){
j <- i
AFAdjustedBeliefs = BLBF_CalculateForecastAdjustedBeliefs(SOS = AFSOS , ExtractedData$InitialArray[j , ])
NAFAdjustedBeliefs = BLBF_CalculateForecastAdjustedBeliefs(SOS = NAFSOS , ExtractedData$InitialArray[j , ])

ImMatrix[i,1]  = BLBF_CalculateAdjustedDiscrepancy(AFAdjustedBeliefs , ExtractedData$Data[j , ])
ImMatrix[i,2] = BLBF_CalculateAdjustedDiscrepancy(NAFAdjustedBeliefs , ExtractedData$Data[j , ])

DP_WaitBar(i/a)

}


ImMatrix[,1] <- (ImMatrix[,1] - mean(ImMatrix[ExtractedData$AFLogical ==1 , 1]))/sqrt( var(ImMatrix[ExtractedData$AFLogical ==1 , 1] ) )
ImMatrix[,2] <- (ImMatrix[,2] - mean(ImMatrix[ExtractedData$AFLogical ==0 , 2]))/sqrt( var(ImMatrix[ExtractedData$AFLogical ==0 , 2] ) )

x11()
par(mfrow = c(1 , 2))
BC_PlotCompareSingleHists(ImMatrix[ExtractedData$AFLogical == 0 , 1] , ImMatrix[ExtractedData$AFLogical , 1] , main = 'ImAF')
BC_PlotCompareSingleHists(ImMatrix[ExtractedData$AFLogical == 0 , 2] , ImMatrix[ExtractedData$AFLogical , 2], main = 'NAF')


SampleofPoints1 <- sample(which(ExtractedData$AFLogical ) , 10000 , replace = F)
SampleofPoints2 <- sample(which(ExtractedData$AFLogical ==0) , 10000 , replace = F)

x11()
plot(ImMatrix[SampleofPoints1 , 1] , ImMatrix[SampleofPoints1 , 2] , col = rgb(1 , 0 , 0 , alph = 0.01) , xlim = c(-3,8) , ylim = c(-4 , 5) , pch = 16 , xlab = 'Implausibility AFib' , ylab = 'Implausibility Not-AFib')
points(ImMatrix[SampleofPoints2 , 1] , ImMatrix[SampleofPoints2 ,2] , col = rgb(0 , 0 , 1 , alph = 0.01) , pch = 16)
#abline(-3 , 0)
abline(3 , 0)
#abline( v = -3 )
abline( v = 3 )
title('Implausibility Pairs')

SampleofPoints1 <- sample(which( ((ExtractedData$AFLogical)*(ImMatrix[,1] <3)) == 1 ) , 100000 , replace = F)
AFModel <- kde( ImMatrix[SampleofPoints1 , ]  )

SampleofPoints1 <- sample(which( ((ExtractedData$AFLogical ==F)*(ImMatrix[,1] <3)) == 1 ) , 100000, replace = F)
NAFModel <- kde( ImMatrix[SampleofPoints1 , ]  )

plot(NAFModel , col = 'blue' , add = T)
plot(AFModel , col = 'red' , add = T)

F_1 <- BLBF_CalculateDensities(AFModel = AFModel , NAFModel = NAFModel , ImMatrix = ImMatrix[ImMatrix[,1] <3 , ]  )

#alpha <- c(sum(ExtractedData$AFLogical)/length(ExtractedData$AFLogical) ,  1- sum(ExtractedData$AFLogical)/length(ExtractedData$AFLogical))

alpha <- c(sum(ExtractedData$AFLogical[ImMatrix[,1] <3])/length(ExtractedData$AFLogical[ImMatrix[,1] <3]) ,  1- sum(ExtractedData$AFLogical[ImMatrix[,1] <3])/length(ExtractedData$AFLogical[ImMatrix[,1] <3]))

IntialData <- distinct(data.frame(ExtractedData$InitialArray , ExtractedData$AFLogical) ) 


PosteriorProbabilities <- (alpha[1]*F_1[,1]) / ((alpha[1]*F_1[,1]) + (alpha[2]*F_1[,2]))
PosteriorProbabilities[is.na(PosteriorProbabilities)] <- alpha[1] 

x11()
EmpericalProbabilityStructure2 <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(PosteriorProbabilities , ExtractedData$AFLogical[ImMatrix[,1] < 3]), BinWidth = 0.1))
print(BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure2) + ggtitle('Emperical Probabilities'))

x11()
PerformanceSweep <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(PosteriorProbabilities , ExtractedData$AFLogical[ImMatrix[,1] < 3]))
ROCplot <- BC_PlotsCreateROC(PerformanceSweep) + ggtitle('Forecasting ROC Curves') + geom_abline(intercept = 0 , slope = 1)
NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep) + geom_vline(xintercept = alpha[2]) + geom_hline(yintercept = alpha[1])+  ggtitle('Forecasting PPV vs NPV Curves')
grid.arrange(ROCplot , NPVPPVPlot , nrow= 2)
