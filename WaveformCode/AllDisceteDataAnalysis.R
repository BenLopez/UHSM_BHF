{
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
    set.seed( 1 )
  }
}

##### Preoperation Data #####

# PatIndex2017, DetIndex2017 , BioChemIndex2017, HaemIndex2017 
#

source('ADDAProcessData.R')
#### Begin Data Analysis #####

# Pre_op Catagorical Data Structure

{
  NamesofPeropCategoricalVariables <- c( 'Gender' ,
                                         'Valve',
                                         'CABG',
                                         'Aortic',
                                         'Complex',
                                         "EjectionFractionCategory",
                                         "NYHAGrade" ,
                                         "AnginaGrade",
                                         'Urgency',
                                         "Active.Endocarditis",
                                         "HTN",
                                         "HistoryOfNeurologicalDysfunction",
                                         "HistoryOfPulmonaryDisease",
                                         "Diabetes",
                                         "Thoracic.Aorta",
                                         "PreviousCardiacSurgery",
                                         "Recent.MI",
                                         "PreOpSupport",
                                         "VentilatedPreOperation",
                                         "CardiogenicShock_pre_Operation",
                                         "IntravenousInotropesPriorToAnaesthesia",
                                         "IntravenousNitratesOrAnyHeparin",
                                         "ExtracardiacArteriopathy",
                                         "Planned.Valve.Surgery")
  
  NameofVariable <- NamesofPeropCategoricalVariables[1]
  PreOpCategoricalUnivariateResults <- POM_CategoricalUnivariateAnalysis(CategoricalData = MasterData[ , which(names(MasterData) == NameofVariable)] ,
                                                                         AFLogical = MasterData$AFLogical , MasterData , NameofVariable )
  for(i in 2:length(NamesofPeropCategoricalVariables)){
    NameofVariable <- NamesofPeropCategoricalVariables[i]
    PreOpCategoricalUnivariateResults <- rbind(PreOpCategoricalUnivariateResults,
                                               POM_CategoricalUnivariateAnalysis(CategoricalData = MasterData[ , which(names(MasterData) == NameofVariable)] ,
                                                                                AFLogical = MasterData$AFLogical , MasterData , NameofVariable ))
  }
  
  PreOpCategoricalUnivariateResults[is.na(PreOpCategoricalUnivariateResults)] <- 'NA'
  PreOpCategoricalUnivariateResults[ which(is.na(PreOpCategoricalUnivariateResults[,9])) , 9] <- 'NA'
  PreOpCategoricalUnivariateResults[ -which(PreOpCategoricalUnivariateResults[,9] < (1e-06)) , 9] <- 'NA' # no idea about -which() here
  PreOpCategoricalUnivariateResults[ which(PreOpCategoricalUnivariateResults[,10] > (1e02)) , 10] <- 'NA'
  PreOpCategoricalUnivariateResults[ which(PreOpCategoricalUnivariateResults[,8] > (1e02)) , 8] <- 'NA'
  
  PreOpCategoricalUnivariateResultsLatex <- xtable(PreOpCategoricalUnivariateResults)
  
  # Pre_op Continuous Data Structure
  
  
  NamesofPreopContinuousVariables <-c('Age' ,
                                      'Weight' , 
                                      'BMI',
                                      "LogisticEUROScore", 
                                      'PreOpNa',
                                      'PreOpK',
                                      'PreopMg',
                                      'PreopUrea' ,
                                      'PreOpCRP',
                                      'PreOpAlb' ,
                                      'PreopBili' ,
                                      'PreopCreat',
                                      "PreopHb" ,
                                      "PreopPLT" ,
                                      "PreopWBC",
                                      "PreopPT",
                                      "PreopAPTT")
  
  NameofVariable <- NamesofPreopContinuousVariables[1]
  ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
  PreOpContinuousUnivariateResults <- POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable)
  for(i in 2:length(NamesofPreopContinuousVariables)){
    NameofVariable <- NamesofPreopContinuousVariables[i]
    ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
    PreOpContinuousUnivariateResults <- rbind(PreOpContinuousUnivariateResults ,POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable) )
  }
  
  PreOpContinuousUnivariateResultsLatex <- xtable(PreOpContinuousUnivariateResults)
  
  # Post_op Catagorical Data Structure
  NamesofPostopCategoricalVariables <- c(  'Filter',
                                          'IABP2')
  
  
  
  NameofVariable <- NamesofPostopCategoricalVariables[1]
  PostOpCategoricalUnivariateResults <- POM_CategoricalUnivariateAnalysis(CategoricalData = MasterData[ , which(names(MasterData) == NameofVariable)] ,
                                                                          AFLogical = MasterData$AFLogical , MasterData , NameofVariable )
  for(i in 2:length(NamesofPostopCategoricalVariables)){
    NameofVariable <- NamesofPostopCategoricalVariables[i]
    PostOpCategoricalUnivariateResults <- rbind(PostOpCategoricalUnivariateResults,
                                                POM_CategoricalUnivariateAnalysis(CategoricalData = MasterData[ , which(names(MasterData) == NameofVariable)] ,
                                                                                  MasterData$AFLogical , MasterData , NameofVariable ))
  }
  
  PostOpCategoricalUnivariateResults[is.na(PostOpCategoricalUnivariateResults)] <- 'NA'
  PostOpCategoricalUnivariateResults[ which(is.na(PostOpCategoricalUnivariateResults[,9])) , 9] <- 'NA'
  PostOpCategoricalUnivariateResults[ -which(PostOpCategoricalUnivariateResults[,9] < (1e-06)) , 9] <- 'NA' # no idea about -which() here
  PostOpCategoricalUnivariateResults[ which(PostOpCategoricalUnivariateResults[,10] > (1e02)) , 10] <- 'NA'
  PostOpCategoricalUnivariateResults[ which(PostOpCategoricalUnivariateResults[,8] > (1e02)) , 8] <- 'NA'
  
  
  PostOpCategoricalUnivariateResultsLatex <- xtable(PostOpCategoricalUnivariateResults)
  
  
  # Post_op Continuous Data Structure
  
  NamesofPostopContinuousVariables <- c('CPB',
                                        'Na', 
                                        'K',
                                        'Urea',
                                        'Creatinine',
                                        'CRP',
                                        'Albumin',
                                        'Bilirubin',
                                        'Mg',
                                        'WBC',
                                        'Hb',
                                        'Platelets',
                                        'PT',
                                        'APTT',
                                        'INR',
                                        'Fibrinogen',
                                        'HR',
                                        "ReliableART.S",
                                        "ReliableART.M",
                                        "ReliableART.D",
                                        "CVP",
                                        "SpO2",
                                        'dNa',
                                        'dK',
                                        'dMg',
                                        'dUrea' ,
                                        'dCRP',
                                        'dAlb' ,
                                        'dBili' ,
                                        'dCreat',
                                        "dHb" ,
                                        "dPLT" ,
                                        "dWBC",
                                        "dPT",
                                        "dAPTT",
                                        "SOFA",
                                        'logCASUS',
                                        'NoradrenalineStandardised',
                                        'DopamineStandardised',
                                        'AdrenalineStandardised',
                                        'MilrinoneStandardised',
                                        'DobutamineStandardised',
                                        'ArtPO2',
                                        'Lac',
                                        'FiO26num',
                                        'GCSnum',
                                        'PaO2OverFiO2')
  
  NameofVariable <- NamesofPostopContinuousVariables[1]
  ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
  PostOpContinuousUnivariateResults <- POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable)
  for(i in 2:length(NamesofPostopContinuousVariables)){
    NameofVariable <- NamesofPostopContinuousVariables[i]
    ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
    PostOpContinuousUnivariateResults <- rbind(PostOpContinuousUnivariateResults ,POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable) )
  }
  PostOpContinuousUnivariateResults[is.na(PostOpContinuousUnivariateResults)] <- 'NA'
  PostOpContinuousUnivariateResults[ which(is.na(PostOpContinuousUnivariateResults[,9])) , 9] <- 'NA'
  PostOpContinuousUnivariateResults[ -which(PostOpContinuousUnivariateResults[,9] < (1e-06)) , 9] <- 'NA' # no idea about -which() here
  PostOpContinuousUnivariateResults[ which(PostOpContinuousUnivariateResults[,10] > (1e02)) , 10] <- 'NA'
  PostOpContinuousUnivariateResults[ which(PostOpContinuousUnivariateResults[,8] > (1e02)) , 8] <- 'NA'
  
  
  PostOpContinuousUnivariateResultsLatex <- xtable(PostOpContinuousUnivariateResults)
}


{
  print( PreOpCategoricalUnivariateResultsLatex )
  print( PreOpContinuousUnivariateResultsLatex )
  print( PostOpCategoricalUnivariateResultsLatex ) 
  print( PostOpContinuousUnivariateResultsLatex ) 
}

#### Group Models #####


#### Comparison of risk models with models built using underlying parameters ####

# LES

{NamesofLogisticEuroScoreVariables <- c( 'Age',
                                         'Gender',
                                         'PreopCreat',
                                         'Valve',
                                         'CABG',
                                         'Aortic',
                                         'Complex',
                                         "EjectionFractionCategory",
                                         "AnginaGrade",
                                         'Urgency',
                                         "Active.Endocarditis",
                                         "HTN",
                                         "HistoryOfNeurologicalDysfunction",
                                         "HistoryOfPulmonaryDisease",
                                         "Thoracic.Aorta",
                                         "PreviousCardiacSurgery",
                                         "Recent.MI",
                                         "PreOpSupport",
                                         "VentilatedPreOperation",
                                         "CardiogenicShock_pre_Operation",
                                         "IntravenousInotropesPriorToAnaesthesia",
                                         "IntravenousNitratesOrAnyHeparin",
                                         "ExtracardiacArteriopathy",
                                         "Planned.Valve.Surgery"
)

EuroComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofLogisticEuroScoreVariables ,MasterData = MasterData , ControlModel  = 'LogisticEUROScore')


DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = EuroComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


}

# SOFA

{NamesofSofaScoreVariables <- c( 'FiO26num',
                                 'ArtPO2',
                                 'VentilatedInOpLogical',
                                 'GCSnum',
                                 'Platelets',
                                 'Creatinine',
                                 'Bilirubin',
                                 'NoradrenalineStandardised',
                                 'DopamineStandardised',
                                 'AdrenalineStandardised',
                                 'MilrinoneStandardised',
                                 'DobutamineStandardised',
                                 'ReliableART.M',
                                 'PaO2OverFiO2'
)
 
SOFAComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofSofaScoreVariables ,MasterData = MasterData , ControlModel  = 'SOFA')

model <-  glm( formula = SOFAComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

}

# Log Casus

{NamesofLogCasusScoreVariables <- c( 'FiO26num',
                                     'ArtPO2',
                                     'VentilatedInOpLogical',
                                     'Platelets',
                                     'Creatinine',
                                     'Bilirubin',
                                     'Filter',
                                     'Lac',
                                     'Filter',
                                     'IABP2',
                                     'PaO2OverFiO2',
                                     'GCSnum'
)
  
LogCasusComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofLogCasusScoreVariables ,MasterData = MasterData , ControlModel  = 'logCASUS')

model <-  glm( formula = LogCasusComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


}

#### Group Analysis #####
{
{
  NamesPreOpBleedingVariables <- c('PreOpCRP',
                                 'PreopHb',
                                 'PreopPLT',
                                 'PreopWBC',
                                 'PreopPT',
                                 'PreopAPTT'
)
PreOpBleedingComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesPreOpBleedingVariables ,MasterData = MasterData )

model <-  glm( formula = PreOpBleedingComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPreOpCardioVascularVariables <- c('EjectionFractionCategory',
                                       'NYHAGrade',
                                       'AnginaGrade',
                                       'Active.Endocarditis',
                                       'HistoryOfNeurologicalDysfunction',
                                       'HistoryOfPulmonaryDisease',
                                       'Diabetes',
                                       'PreviousCardiacSurgery',
                                       'Recent.MI',
                                       'IntraAorticBallPump',
                                       'VAD',
                                       'CardiogenicShock_pre_Operation',
                                       'IntravenousInotropesPriorToAnaesthesia',
                                       'IntravenousNitratesOrAnyHeparin',
                                       'ExtracardiacArteriopathy'
)
PreOpCardioComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesPreOpCardioVascularVariables ,MasterData = MasterData )

DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = PreOpCardioComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPreOpDemographicVariables <- c('Gender',
                                    'Age',
                                    'Weight',
                                    'BMI'
)
PreOpDemographicComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpDemographicVariables ,MasterData = MasterData )

DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = PreOpDemographicComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

NamesPreOpElectrolytesVariables <- c('PreOpNa',
                                     'PreOpK',
                                     'PreopMg'
)

PreOpElectrolytesComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpElectrolytesVariables ,MasterData = MasterData )
DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = PreOpElectrolytesComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

NamesPreOpInflamatoryVariables <- c('PreOpCRP' ,
                                    'PreOpAlb' ,
                                    'PreopHb' ,
                                    'PreopPLT' ,
                                    'PreopWBC' ,
                                    'PreopPT' ,
                                    'PreopAPTT'
)
PreOpInflamatoryComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpInflamatoryVariables ,MasterData = MasterData )
DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = PreOpInflamatoryComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

NamesPreOpLiverVariables <- c('PreOpAlb',
                              'PreopBili'
)
PreOpLiverComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpLiverVariables ,MasterData = MasterData )

model <-  glm( AFLogical ~  PreOpAlb + PreopBili , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

NamesPreOpOperativeVariables <- c('Valve',
                                  'CABG',
                                  'Aortic',
                                  'Complex',
                                  'Urgency',
                                  'Thoracic.Aorta',
                                  "Planned.Valve.Surgery"
                                  
)
PreOpOperativeComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpOperativeVariables ,MasterData = MasterData )

model <-  glm( formula = PreOpOperativeComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)



NamesPreOpRenalVariables <- c('PreopUrea',
                              'PreopCreat'
)
PreOpRenalComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpRenalVariables ,MasterData = MasterData )
DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = PreOpRenalComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPreOpRespiratoryVariables <- c('VentilatedPreOperation'
                                    )
#PreOpRespiratoryComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpRespiratoryVariables ,MasterData = MasterData )
#}

#### Post Op Groups
{
{
  NamesPostOpBleedingVariables <- c('Hb',
                                 'Platelets',
                                 'PT',
                                 'APTT',
                                 'INR',
                                 'Fibrinogen',
                                 'dHb',
                                 'dPLT',
                                 'dPT',
                                 'dAPTT'
)
PostOpBleedingComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpBleedingVariables,NamesPreOpBleedingVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpBleedingComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)
}
{
NamesPostOpCardioVascularVariables <- c('IABP2',
                                        'HR',
                                        'ReliableART.S',
                                        'ReliableART.M',
                                        'ReliableART.D',
                                        'CVP',
                                        'NoradrenalineStandardised',
                                        'DopamineStandardised',
                                        'AdrenalineStandardised',
                                        'MilrinoneStandardised',
                                        'DobutamineStandardised',
                                        'Lac')
PostOpCardioVascularComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpCardioVascularVariables,NamesPreOpCardioVascularVariables) ,MasterData = MasterData )
DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = PostOpCardioVascularComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

{
NamesPostOpElectrolytesVariables <- c('Na',
                                      'K',
                                      'Mg',
                                      'dNa',
                                      'dK',
                                      'dMg'
  
)
PostOpElectrolytesComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpElectrolytesVariables,NamesPreOpElectrolytesVariables) ,MasterData = MasterData )

DataForLogistic <- POM_SampledImputation(MasterData)
model <-  glm( formula = PostOpElectrolytesComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPostOpInflamatoryVariables <- c('Hb',
                                    'Platelets',
                                    'PT',
                                    'APTT',
                                    'INR',
                                    'Fibrinogen',
                                    'CRP',
                                    'Albumin',
                                    'WBC',
                                    'dHb',
                                    'dPLT',
                                    'dPT',
                                    'dAPTT',
                                    'dCRP',
                                    'dWBC',
                                    'dAlb'
)

PostOpInflamatoryComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpInflamatoryVariables,NamesPreOpInflamatoryVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpInflamatoryComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPostOpLiverVariables <- c('Albumin',
                               'Bilirubin',
                               'dAlb',
                               'dBili')
PostOpLiverComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpLiverVariables,NamesPreOpLiverVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpLiverComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

#NamesPostOpNeuroVariables <- c('GCSnum'
#                               )
#PostOpNeuroComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpNeuroVariables) ,MasterData = MasterData )


NamesPostOpOperativeVariables <- c('CPB',
                                  'IABP2')

PostOpOperativeComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpOperativeVariables,NamesPreOpOperativeVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpOperativeComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPostOpRenalVariables <- c('CVVHDyalysis',
                               'Urea',
                               'Creatinine',
                               'dUrea',
                               'dCreat')
PostOpRenalComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpRenalVariables,NamesPreOpRenalVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpRenalComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

NamesPostOpRespiratoryVariables <- c('SpO2',
                                     'ArtPO2',
                        model <-  glm( formula = PostOpRespiratoryComparisonOutputs$ModelAUC, family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)             'FiO26num',
                                     'PaO2OverFiO2')

PostOpRespiratoryComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpRespiratoryVariables,NamesPreOpRespiratoryVariables) ,MasterData = MasterData )


}
}
}

#### Pre_op Predictive Model #####

{NamesofPeropVariables <- c(NamesofPeropCategoricalVariables , NamesofPreopContinuousVariables[c(-2,-4)] )
IndiciesofPeropVariables <- which(names(MasterData) %in% NamesofPeropVariables)
NamesofPeropVariables <- names(MasterData)[IndiciesofPeropVariables]

{matrixforstep <- diag(length(IndiciesofPeropVariables))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAUC))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAUC))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAUC))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAUC))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAUC))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAUC))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAUC))
matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = EuroComparisonOutputs$ModelAUC))

#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = EuroComparisonOutputs$ModelAIC))
}

PreOpStepOutput <- FC_StepWiseForwardAUC(PreoperativeIndices = IndiciesofPeropVariables , MasterData , matrixforstep)

AUC2 <- FC_CalculateCrossValidatedROC(formulaformodel = PreOpStepOutput[[2]] , PreoperativeIndices = IndiciesofPeropVariables , MasterData)

#### Post_op Predictive Model #####

NamesofPostopVariables <- c(NamesofPostopCategoricalVariables ,NamesofPostopContinuousVariables[c(-37,-36)] )
NamesofPostopVariables <- c(NamesofPostopVariables , NamesofPeropVariables)
IndiciesofPostopVariables <- which(names(MasterData) %in% NamesofPostopVariables)
NamesofPostopVariables <- names(MasterData)[IndiciesofPostopVariables]

{matrixforstep <- diag(length(IndiciesofPostopVariables))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = EuroComparisonOutputs$ModelAUC))

#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAIC))
#matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = EuroComparisonOutputs$ModelAIC))

matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = SOFAComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = LogCasusComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpCardioVascularComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpElectrolytesComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpBleedingComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpInflamatoryComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpLiverComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpOperativeComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRenalComparisonOutputs$ModelAUC))
matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRespiratoryComparisonOutputs$ModelAUC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpStepOutput[[2]]))

#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = SOFAComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = LogCasusComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpCardioVascularComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpElectrolytesComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpBleedingComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpInflamatoryComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpLiverComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpOperativeComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRenalComparisonOutputs$ModelAIC))
#matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRespiratoryComparisonOutputs$ModelAIC))

}
PostOpStepOutput <- FC_StepWiseForwardAUC(PreoperativeIndices = IndiciesofPostopVariables , MasterData , matrixforstep)
}

##### Model validation
{
{
model <- glm(formula = PreOpStepOutput[[2]]
             , family = binomial(link = "logit"),  data = DataForLogistic)
summary(model)
xtable(model)

  LogisticProbility <- predict(model , DataForLogistic , type = c('response'))
  
  PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(as.matrix(LogisticProbility) ,MasterData$AFLogical == 1))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1) + xlab('1 - Specificity')
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) #+ geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')
  
  print(BC_CalculateAreaUnderCurve(PerformanceSweep3))
  
  phist <- ggplot(data.frame( p =LogisticProbility[MasterData$AFLogical ==0]) , aes(p)) + geom_density(col = rgb(0,0,1,alpha = 0.5)) +
    geom_density(data = data.frame( p =LogisticProbility[MasterData$AFLogical ==1]) , aes(p), col = rgb(1,0,0,alpha = 0.5)) + 
    ggtitle('Histogram of Model Output')
  
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$AFLogical == 1), BinWidth = 0.1))
  EmpericalProbabilityPlot <- BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities')
  
  
  
  x11(20,20)
  grid.arrange(phist , EmpericalProbabilityPlot , ROCplot , NPVPPVPlot + ylim(0.2,1)  , ncol =2 , nrow = 2)
}

{
  
  model <- glm(formula = PostOpStepOutput[[2]]
               , family = binomial(link = "logit"),  data = DataForLogistic)
  summary(model)
  xtable(model)
  
  LogisticProbility <- predict(model , DataForLogistic , type = c('response'))
  
  PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(as.matrix(LogisticProbility) ,MasterData$AFLogical == 1))
  ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1) + xlab('1 - Specificity')
  NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) #+ geom_vline(xintercept =( 1-0.19)) + geom_hline(yintercept = 0.19)+  ggtitle('PPV vs NPV Curves')
  
  print(BC_CalculateAreaUnderCurve(PerformanceSweep3))
  
  phist <- ggplot(data.frame( p =LogisticProbility[MasterData$AFLogical ==0]) , aes(p)) + geom_density(col = rgb(0,0,1,alpha = 0.5)) +
    geom_density(data = data.frame( p =LogisticProbility[MasterData$AFLogical ==1]) , aes(p), col = rgb(1,0,0,alpha = 0.5)) + 
    ggtitle('Histogram of Model Output')
  
  EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$AFLogical == 1), BinWidth = 0.1))
  EmpericalProbabilityPlot <- BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities')
  
  
  
  x11(20,20)
  grid.arrange(phist , EmpericalProbabilityPlot , ROCplot , NPVPPVPlot + ylim(0.2,1) , ncol =2 , nrow = 2)
}
}

#### Comparison of Models ####

PreOpModelOutputs <- POM_CalulateProbabilitiesFromModel(PreOpBleedingComparisonOutputs$ModelAUC, DataForLogistic)
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PreOpDemographicComparisonOutputs$ModelAUC, DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PreOpCardioComparisonOutputs$ModelAUC, DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PreOpElectrolytesComparisonOutputs$ModelAUC, DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PreOpInflamatoryComparisonOutputs$ModelAUC, DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PreOpLiverComparisonOutputs$ModelAUC, DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PreOpOperativeComparisonOutputs$ModelAUC, DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(as.formula('AFLogical ~ PreopUrea'), DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(EuroComparisonOutputs$ModelAUC, DataForLogistic))
PreOpModelOutputs <- cbind(PreOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PreOpStepOutput[[2]], DataForLogistic))


PreOpC <- cov2cor(cov(PreOpModelOutputs))
rownames(PreOpC) <- c('Bleeding' , 'Demograhic' , 'Cardiovascular' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal', 'EuroScore', 'PreopModel')
colnames(PreOpC) <- c('Bleeding' , 'Demograhic' , 'Cardiovascular' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal', 'EuroScore', 'PreopModel')


BC_PlotPairsFromTwoVariables(PreOpModelOutputs[AFLogical ==1,1:8] ,PreOpModelOutputs[sample(which(AFLogical ==0),260),1:8],
             alpha = 0.1 ,
             labels = c('Bleeding' , 'Demograhic' , 'Cardiovascular' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal'))

PostOpModelOutputs <- POM_CalulateProbabilitiesFromModel(PostOpBleedingComparisonOutputs$ModelAUC, DataForLogistic)
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpCardioVascularComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpElectrolytesComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpInflamatoryComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpLiverComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpOperativeComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpRenalComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpRespiratoryComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(SOFAComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind(PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpStepOutput[[2]], DataForLogistic))


BC_PlotPairsFromTwoVariables(PostOpModelOutputs[AFLogical ==1,1:9] ,PostOpModelOutputs[sample(which(AFLogical ==0),260),1:9],
                             alpha = 0.1 ,
                             labels = c('Bleeding' , 'Demograhic' , 'Cardiovascular' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory', 'SOFA', 'PostOpModel'))

PostOpC <- cov2cor(cov(PostOpModelOutputs))
rownames(PostOpC) <- c('Bleeding' , 'Demograhic' , 'Cardio' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory')
colnames(PostOpC) <- c('Bleeding' , 'Demograhic' , 'Cardio' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory')

BC_PlotPairsFromTwoVariables(PostOpModelOutputs[AFLogical ==1,1:9] ,PostOpModelOutputs[sample(which(AFLogical ==0),260),1:9],
                             alpha = 0.1 ,
                             labels = c('Bleeding' , 'Demograhic' , 'Cardio' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory', 'SOFA', 'PostOpModel'))
Total <- cbind(PreOpModelOutputs[, c(1:8 , 10)] , PostOpModelOutputs[, c(1:8 , 10)])
x11()
SampleofPoints<- sample(which(AFLogical ==0),260)
  p5 <- ggplot(data.frame(x = Total[AFLogical ==1,c(9)] , y= Total[AFLogical ==1,c(18)] ) , aes(x , y)) + 
    geom_point(col = rgb(1,0,0,alpha = 0.25)) +
    geom_point(data = data.frame(x =Total[SampleofPoints,c(9)] , y = Total[SampleofPoints,c(18)] ) , aes(x , y), col = rgb(0,0,1,alpha = 0.25)) +
    xlab('Pre-operative Model Output') +
    ylab('Post-operative Model Output') +
    xlim(0,1)+
    ylim(0,1)+
    geom_abline(intercept = 0,slope = 1)
print(p5)

BC_PlotPairsFromTwoVariables(Total[AFLogical ==1,] ,Total[sample(which(AFLogical ==0),260),],
                             alpha = 0.1 , labels = c(c('Bleeding' , 'Demographic' , 'Cardiovascular' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal') , c('Bleeding'  , 'Cardio' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory' , 'Post-op Model')))


BC_PlotPairsFromTwoVariables(Total[AFLogical ==1,] ,Total[sample(which(AFLogical ==0),260),],
                               alpha = 0.1 , labels = c(c('Bleeding' , 'Demographic' , 'Cardiovascular' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Pre-op Model') , c('Bleeding' , 'Cardio' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory','Post-op Model')),
                             horInd = c(10:18),
                             verInd = c(1:9))
}
  