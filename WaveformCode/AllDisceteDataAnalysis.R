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
                                         'ProcDetails' ,
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
                                                                                 MasterData$AFLogical , MasterData , NameofVariable ))
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
  NamesofPostopCategoricalVariables <- c( 'IntubatedInSurgury' ,
                                          'Filter',
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
                                        'VasopressinStandardised',
                                        'DobutamineStandardised',
                                        'ArtPO2',
                                        'Lac',
                                        'FiO26num',
                                        'PAR',
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
{
# LES

{NamesofLogisticEuroScoreVariables <- c( 'Age',
                                         'Gender',
                                         'PreopCreat',
                                         'ProcDetails' ,
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
                                 'VasopressinStandardised',
                                 'DobutamineStandardised',
                                 'ReliableART.M',
                                 'PaO2OverFiO2'
)
 
SOFAComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofSofaScoreVariables ,MasterData = MasterData , ControlModel  = 'SOFA')
}

# Log Casus

{NamesofLogCasusScoreVariables <- c( 'FiO26num',
                                     'ArtPO2',
                                     'VentilatedInOpLogical',
                                     'Platelets',
                                     'Creatinine',
                                     'Bilirubin',
                                     'Filter',
                                     'PAR',
                                     'Lac',
                                     'Filter',
                                     'IABP2',
                                     'PaO2OverFiO2',
                                     'GCSnum'
)
  
LogCasusComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofLogCasusScoreVariables ,MasterData = MasterData , ControlModel  = 'logCASUS')
}

#### Group Analysis #####
{
  NamesPreOpBleedingVariables <- c('PreOpCRP',
                                 'PreopHb',
                                 'PreopPLT',
                                 'PreopWBC',
                                 'PreopPT',
                                 'PreopAPTT'
)
PreOpBleedingComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesPreOpBleedingVariables ,MasterData = MasterData )


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


NamesPreOpDemographicVariables <- c('Gender',
                                    'Age',
                                    'Weight',
                                    'BMI'
)
PreOpDemographicComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpDemographicVariables ,MasterData = MasterData )

NamesPreOpElectrolytesVariables <- c('PreOpNa',
                                     'PreOpK',
                                     'PreopMg'
)

PreOpElectrolytesComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpElectrolytesVariables ,MasterData = MasterData )
NamesPreOpInflamatoryVariables <- c('PreOpCRP',
                                    'PreOpAlb',
                                    'PreopHb',
                                    'PreopPLT',
                                    'PreopWBC',
                                    'PreopPT',
                                    'PreopAPTT'
)
PreOpInflamatoryComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpInflamatoryVariables ,MasterData = MasterData )

NamesPreOpLiverVariables <- c('PreOpAlb',
                              'PreopBili'
)
PreOpLiverComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpLiverVariables ,MasterData = MasterData )

NamesPreOpOperativeVariables <- c('ProcDetails',
                                  'Urgency',
                                  'Thoracic.Aorta'
)
PreOpOperativeComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpOperativeVariables ,MasterData = MasterData )

NamesPreOpRenalVariables <- c('PreopUrea',
                              'PreopCreat'
)
PreOpRenalComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpRenalVariables ,MasterData = MasterData )

NamesPreOpRespiratoryVariables <- c('VentilatedPreOperation'
                                    )
PreOpRespiratoryComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpRespiratoryVariables ,MasterData = MasterData )
}

#### Post Op Groups
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

model <-  glm( formula = PostOpBleedingComparisonOutputs$ModelAIC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


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
                                        'VasopressinStandardised',
                                        'DobutamineStandardised',
                                        'Lac',
                                        'PAR')
PostOpCardioVascularComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpCardioVascularVariables,NamesPreOpCardioVascularVariables) ,MasterData = MasterData )


NamesPostOpElectrolytesVariables <- c('Na',
                                      'K',
                                      'Mg',
                                      'dNa',
                                      'dK',
                                      'dMg'
  
)
PostOpElectrolytesComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpElectrolytesVariables,NamesPreOpElectrolytesVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpElectrolytesComparisonOutputs$ModelAIC , family=binomial(link='logit') , data=DataForLogistic)
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

model <-  glm( formula = PostOpInflamatoryComparisonOutputs$ModelAIC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPostOpLiverVariables <- c('Albumin',
                               'Bilirubin',
                               'dAlb',
                               'dBili')
PostOpLiverComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpLiverVariables,NamesPreOpLiverVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpLiverComparisonOutputs$ModelAIC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)

NamesPostOpNeuroVariables <- c('GCSnum'
                               )
PostOpNeuroComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpNeuroVariables) ,MasterData = MasterData )


NamesPostOpOperativeVariables <- c('CPB',
                                  'IABP2')

PostOpOperativeComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpOperativeVariables,NamesPreOpOperativeVariables) ,MasterData = MasterData )

model <-  glm( formula = PostOpOperativeComparisonOutputs$ModelAIC , family=binomial(link='logit') , data=DataForLogistic)
summary(model)  
xtable(model)


NamesPostOpRenalVariables <- c('CVVHDyalysis',
                              'Urea',
                              'Creatinine',
                              'dUrea',
                              'dCreat')
PostOpRenalComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpRenalVariables,NamesPreOpRenalVariables) ,MasterData = MasterData )

NamesPostOpRespiratoryVariables <- c('IntupatedInSurgery',
                                     'SpO2',
                                     'ArtPO2',
                                     'FiO26num',
                                     'PaO2OverFiO2')

PostOpRespiratoryComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpRespiratoryVariables,NamesPreOpRespiratoryVariables) ,MasterData = MasterData )
}
}


#### Pre_op Predictive Model #####


#NamesofPeropVariables <- c(NamesofPeropCategoricalVariables ,NamesofPreopContinuousVariables )
#IndiciesofPeropVariables <- which(names(MasterData) %in% NamesofPeropVariables)

#PreOpStepOutput <- FC_StepWiseForwardAUC(PreoperativeIndices = IndiciesofPeropVariables , MasterData)


#### Post_op Predictive Model #####


