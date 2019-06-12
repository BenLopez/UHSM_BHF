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
  #listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

{
  MasterPreOpData <- POM_CreateDataStructure(PatIndex2017 , DetIndex2017 , BioChemIndex2017)
  ReducedMasterPreOpData <- MasterPreOpData[MasterPreOpData$Pre_OperativeHeartRhythm != 'Atrial fibrillation/flutter' , ]
  BLBModelStruct <- POM_CreateBLBModellingStructure(MasterPreOpData =  ReducedMasterPreOpData  , listofcovariates = POM_CreateDefaultSetofCov() )
  PosteriorProbability <- POM_CalculateBLBProbability(MasterPreOpData = ReducedMasterPreOpData , BLBModelStruct, listofcovariates = POM_CreateDefaultSetofCov() , ReducedMasterPreOpData$PseudoId )
  BLBCalibrationStruct <- POM_CreateBLBCalibrationStructure(PosteriorProbability , ReducedMasterPreOpData$AFLogical )
  #POM_CalculatePreOpProbability()
}


POM_ExractPostOpFromBioChem <- function(BioChemIndex2017){

PostOpBioChem <- matrix(0 , length(BioChemIndex2017) , )

for(i in 1:length(BioChemIndex2017)){
PostOpBioChem[i , 1:8] <- as.numeric(BioChemIndex2017[[i]]$TimeSeriesData[1 , 2:9])
}


colnames(PostOpBioChem) <- c(names(BioChemIndex2017[[i]]$TimeSeriesData[1 , 2:9]) , 'AFLogical' )
rownames(PostOpBioChem) <- apply(as.matrix(names(BioChemIndex2017)) , 1 ,names(BioChemIndex2017))
return(PostOpBioChem)
}

PostOpBioChem <- PostOpBioChem[apply(PostOpBioChem , 1 , sum) !=0 , ]
PostOpBioChem <- PostOpBioChem[duplicated(rownames(PostOpBioChem)) == F , ]

PostOpBioChem <- DP_RemoveNaRows(PostOpBioChem)

PostOpBioChem <- data.frame(PostOpBioChem)
PostOpBioChem$AFLogical <- as.factor(PostOpBioChem$AFLogical)

model <- (glm(formula = AFLogical ~ Na + K + Urea + Creatinine + CRP +Albumin + Bilirubin + Mg  ,family=binomial(link='logit') , data=PostOpBioChem ))
summary( model )

model <- (glm(formula = AFLogical ~  as.factor(CRP > 100)  ,family=binomial(link='logit') , data=PostOpBioChem ))
summary( model )

AFLogical <- PostOpBioChem$AFLogical == 1
BC_PlotCompareSingleHists(model$fitted.values[AFLogical == F] , model$fitted.values[AFLogical == T])

DataMatrix <- PostOpBioChem[ , indexesofcovariates]
BC_PlotCompareSingleHists(DataMatrix[AFLogical == F , 1] , DataMatrix[AFLogical == T , 1])


listofcovariates <- c('CRP'  , 'Mg')
indexesofcovariates <-  which(names(PostOpBioChem) %in% listofcovariates )

DataMatrix <- PostOpBioChem[ , indexesofcovariates]

BC_PlotCompareTwoHists(DataMatrix[AFLogical == 1 ,] ,DataMatrix[AFLogical == 0 , ]  )
BC_PlotPairsFromTwoVariables(DataMatrix[AFLogical == 1 ,] ,DataMatrix[which(AFLogical == 0)[sample(1:sum(AFLogical == 0) , sum(AFLogical == 1))] ,]  , alpha = 0.1)

PosteriorProbability <- BLBC_FitBayesLinearBayesClassifier(Data =  DataMatrix , Labels = AFLogical )

listofcovariates <- c('Age'  , 'CPB' ,  'AdditiveEUROScore' , 'PreOpNa', 'PreopUrea' , 'PreOpAlb'  , 'PreopCreat')
indexesofcovariates <-  which(names(ReducedMasterPreOpData) %in%listofcovariates )

ReducedMasterPreOpData <- ReducedMasterPreOpData[ReducedMasterPreOpData$PseudoId %in% rownames(DataMatrix) ,  ]
AFLogical <- AFLogical[rownames(DataMatrix) %in% ReducedMasterPreOpData$PseudoId] 
DataMatrix <- DataMatrix[rownames(DataMatrix) %in% ReducedMasterPreOpData$PseudoId , ]
ReducedMasterPreOpData <- ReducedMasterPreOpData[order(ReducedMasterPreOpData$PseudoId) , ]
AFLogical <- AFLogical[order(rownames(DataMatrix))] 
DataMatrix <- DataMatrix[order(rownames(DataMatrix)) , ]


plot(DataMatrix$Urea[AFLogical] , ReducedMasterPreOpData$PreopUrea[AFLogical] , col = rgb(1 , 0, 0, alpha = 0.2)  , pch =16)
tmp <- sample(which(AFLogical ==0) , sum(AFLogical))
points(DataMatrix$Urea[tmp] , ReducedMasterPreOpData$PreopUrea[tmp] , col = rgb(0 , 0, 1, alpha = 0.2)  , pch =16)
abline(0,1)


plot(DataMatrix$Na[AFLogical] , ReducedMasterPreOpData$PreOpNa[AFLogical] , col = rgb(1 , 0, 0, alpha = 0.2)  , pch =16)
tmp <- sample(which(AFLogical ==0) , sum(AFLogical))
points(DataMatrix$Na[tmp] , ReducedMasterPreOpData$PreOpNa[tmp] , col = rgb(0 , 0, 1, alpha = 0.2)  , pch =16)
abline(0,1)


plot(DataMatrix$CRP[AFLogical] , ReducedMasterPreOpData$PreOpCRP[AFLogical] , col = rgb(1 , 0, 0, alpha = 0.2)  , pch =16)
tmp <- sample(which(AFLogical ==0) , sum(AFLogical))
points(DataMatrix$CRP[tmp] , ReducedMasterPreOpData$PreOpCRP[tmp] , col = rgb(0 , 0, 1, alpha = 0.2)  , pch =16)
abline(0,1)


plot(DataMatrix$Mg[AFLogical] , ReducedMasterPreOpData$PreopMg[AFLogical] , col = rgb(1 , 0, 0, alpha = 0.2)  , pch =16)
tmp <- sample(which(AFLogical ==0) , sum(AFLogical))
points(DataMatrix$Mg[tmp] , ReducedMasterPreOpData$PreopMg[tmp] , col = rgb(0 , 0, 1, alpha = 0.2)  , pch =16)
abline(0,1)

TotalMatrix <- cbind(DataMatrix , ReducedMasterPreOpData[ , indexesofcovariates])
TotalMatrix <- TotalMatrix[ -415 , ]
AFLogical <- AFLogical[-415  ]
tmp <- sample(which(AFLogical ==0) , sum(AFLogical))
BC_PlotPairsFromTwoVariables(TotalMatrix[AFLogical , ] , TotalMatrix[tmp , ] , alpha = 0.1)

PosteriorProbability <- BLBC_FitBayesLinearBayesClassifier(Data = as.matrix(TotalMatrix) , Labels = AFLogical )
