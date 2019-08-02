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

{NamesofPeropCategoricalVariables <- c( 'Gender' ,
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
PreOpCategoricalUnivariateResults[ is.na(PreOpCategoricalUnivariateResults[,9]) , 9] <- 'NA'
PreOpCategoricalUnivariateResults[ PreOpCategoricalUnivariateResults[,10] > 1e06 , 10] <- 'NA'

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
PostOpCategoricalUnivariateResults[ is.na(PostOpCategoricalUnivariateResults[,9]) , 9] <- 'NA'
PostOpCategoricalUnivariateResults[ PostOpCategoricalUnivariateResults[,10] > 1e06 , 10] <- 'NA'

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

NamesofLogisticEuroScoreVariables <- c( 'Age',
                                        'Gender',
                                        'PreopCreat',
                                        'ProcDetails' ,
                                        "EjectionFractionCategory",
                                        "AnginaGrade",
                                        'Urgency'
)

# SOFA

NamesofSofaScoreVariables <- c( 'FiO26num',
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

IndiciesofSofaScoreVariables <- which(names(MasterData) %in% NamesofSofaScoreVariables)

formulaformodel <- FB_CreateFormula('AFLogical' , IndiciesofSofaScoreVariables , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  ,  ]))
model <- (glm(formula = formulaformodel,family=binomial(link='logit') , data=DataForLogistic))
stepaicoutput <- stepAIC(model)
model <- (glm(formula = stepaicoutput$formula,family=binomial(link='logit') , data=DataForLogistic))
summary( model )
AUCSOFA1 <- FC_CalculateCrossValidatedROC(formulaformodel = stepaicoutput$formula , PreoperativeIndices = IndiciesofSofaScoreVariables , MasterData)

StepOutput <- FC_StepWiseForwardAUC(IndiciesofSofaScoreVariables , MasterData)
AUCSOFA2 <- FC_CalculateCrossValidatedROC(formulaformodel = StepOutput[[2]] , PreoperativeIndices = IndiciesofSofaScoreVariables , MasterData)

formulaformodel <- FB_CreateFormula('AFLogical' , which(names(MasterData) == 'SOFA') , MasterData)
AUCSOFA3 <- FC_CalculateCrossValidatedROC(formulaformodel ,which(names(MasterData) == 'SOFA'), MasterData)

# Log Casus

NamesofLogCasusScoreVariables <- c( 'FiO26num',
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
                                'PaO2OverFiO2'
)

IndiciesofLogCASUSScoreVariables <- which(names(MasterData) %in% NamesofLogCasusScoreVariables)

formulaformodel <- FB_CreateFormula('AFLogical' , IndiciesofLogCASUSScoreVariables , MasterData)

DataForLogistic <- POM_SampledImputation(MasterPreOpData = data.frame(MasterData[  ,  ]))
model <- (glm(formula = formulaformodel,family=binomial(link='logit') , data=DataForLogistic))
stepaicoutput <- stepAIC(model)
model <- (glm(formula = stepaicoutput$formula,family=binomial(link='logit') , data=DataForLogistic))
summary( model )
AUClogCASUS1 <- FC_CalculateCrossValidatedROC(formulaformodel = stepaicoutput$formula , PreoperativeIndices = IndiciesofLogCASUSScoreVariables , MasterData)

StepOutput <- FC_StepWiseForwardAUC(IndiciesofLogCASUSScoreVariables , MasterData)
AUClogCASUS2 <- FC_CalculateCrossValidatedROC(formulaformodel = StepOutput[[2]] , PreoperativeIndices = IndiciesofLogCASUSScoreVariables , MasterData)

formulaformodel <- FB_CreateFormula('AFLogical' , which(names(MasterData) == 'logCASUS') , MasterData)
AUClogCASUS3 <- FC_CalculateCrossValidatedROC(formulaformodel ,which(names(MasterData) == 'logCASUS'), MasterData)


#### Pre_op Predictive Model #####


NamesofPeropVariables <- c(NamesofPeropCategoricalVariables ,NamesofPreopContinuousVariables )
IndiciesofPeropVariables <- which(names(MasterData) %in% NamesofPeropVariables)

PreOpStepOutput <- FC_StepWiseForwardAUC(PreoperativeIndices = IndiciesofPeropVariables , MasterData)


#### Post_op Predictive Model #####


