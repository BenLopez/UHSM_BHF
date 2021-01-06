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

source('ADDAProcessDataReadmission2.R')
#### Begin Data Analysis #####

# Adjust functions to change control variables. 
{  POM_CategoricalUnivariateAnalysis <- function(CategoricalData , AFLogical , MasterData , NameofVariable ){
  
  CategoricalData <- as.factor(CategoricalData)
  
  output <- data.frame(matrix(NA , length(levels(CategoricalData)) , length(POM_CategoricalCreateFieldNames()) ))
  output[,1] <- levels(CategoricalData)
  names(output) <- POM_CategoricalCreateFieldNames()
  
  for(i in 1:length(levels(CategoricalData))){
    # basic analysis  
    tmp <- as.factor(CategoricalData == levels(CategoricalData)[i] )
    output[i , 2] <- as.numeric(summary(tmp)[2])
    output[i , 3] <- as.numeric(summary(tmp[AFLogical == 0])[2])
    output[i , 4] <- as.numeric(summary(tmp[AFLogical == 1])[2])
    output[i , 5] <- summary(tmp[AFLogical == 0])[2]/sum(AFLogical == 0)*100
    output[i , 6] <- summary(tmp[AFLogical == 1])[2]/sum(AFLogical == 1)*100
    testtmp <- prop.test( c(summary(tmp[AFLogical == 1])[1] ,summary(tmp[AFLogical == 0])[1]) , c(  sum(AFLogical == 1), sum(AFLogical == 0) ))
    output[i , 7] <- signif(testtmp$p.value , 3)
    # logistic regression analysis
  }
  
  Formula1 <- FB_CreateFormula('AFLogical' , unique(c(which(names(MasterData) == NameofVariable), which(names(MasterData) == 'Valve'), which(names(MasterData) == 'CABG'), which(names(MasterData) == 'Aortic') ,  which(names(MasterData) == 'Complex') ,which(names(MasterData) == 'logCASUS'),which(names(MasterData) == 'CPB')  ) ) , MasterData )
  basevariables <-  unique(c(which(names(MasterData) == NameofVariable), which(names(MasterData) == 'Valve'), which(names(MasterData) == 'CABG'), which(names(MasterData) == 'Aortic') ,  which(names(MasterData) == 'Complex') ,which(names(MasterData) == 'LogisticEUROScore'),which(names(MasterData) == 'CPB')  ) )
  if(sum(basevariables == which(names(MasterData) == NameofVariable))>0){basevariables <- basevariables[-which(basevariables == which(names(MasterData) == NameofVariable))]}
  Formula2 <- FB_CreateFormula('AFLogical' , basevariables , MasterData )
  
  model <- glm(Formula1 ,
               family=binomial(link='logit'),
               data = MasterData)
  
  
  output[ , 8] <- exp(coef(model)[1:length(levels(CategoricalData))])
  output[ , 9:10] <- exp(confint(model, trace = FALSE) )[1:length(levels(CategoricalData)),]
  output[ , 11] <- summary(model)$coefficients[1:length(levels(CategoricalData)),4]  
  output[ , 2:11] <- signif(output[,2:11] , 4)
  
  output[1 , 8:11] <- NA
  output[1 , 12] <- -(FC_CalculateCrossValidatedROC(Formula1 , which(FC_ExtractIndiciesFromFormula(names(MasterData),Formula1) ==1) , MasterData = MasterData)-FC_CalculateCrossValidatedROC(Formula2 , which(FC_ExtractIndiciesFromFormula(names(MasterData),Formula1) ==1) , MasterData = MasterData))
  return(output)
}
POM_ContinuousUnivariateAnalysis <- function(ContinuousData , AFLogical , MasterData , NameofVariable){
  ContinuousData <- as.numeric(ContinuousData)
  output <- data.frame(matrix(NA , 1 , length(POM_ContinuousCreateFieldNames()) ))
  output[,1] <- NameofVariable
  names(output) <- POM_ContinuousCreateFieldNames()
  
  output[,2] <- sum(!is.na(ContinuousData))
  output[,3] <- mean(ContinuousData[AFLogical ==0] , na.rm = T) 
  output[,4] <- sqrt(var(ContinuousData[AFLogical ==0] , na.rm = T))
  output[,5] <- mean(ContinuousData[AFLogical ==1] , na.rm = T) 
  output[,6] <- sqrt(var(ContinuousData[AFLogical ==1] , na.rm = T))
  
  tmptest <- t.test(ContinuousData[AFLogical ==1] , ContinuousData[AFLogical ==0])
  output[,7] <- tmptest$p.value
  
  Formula1 <- FB_CreateFormula('AFLogical' , unique(c(which(names(MasterData) == NameofVariable), which(names(MasterData) == 'Valve'), which(names(MasterData) == 'CABG'), which(names(MasterData) == 'Aortic') ,  which(names(MasterData) == 'Complex') ,which(names(MasterData) == 'logCASUS'),which(names(MasterData) == 'CPB')  ) ) , MasterData )
  basevariables <-  unique(c(which(names(MasterData) == NameofVariable), which(names(MasterData) == 'Valve'), which(names(MasterData) == 'CABG'), which(names(MasterData) == 'Aortic') ,  which(names(MasterData) == 'Complex') ,which(names(MasterData) == 'LogisticEUROScore'),which(names(MasterData) == 'CPB')  ) )
  if(sum(basevariables == which(names(MasterData) == NameofVariable))>0){basevariables <- basevariables[-which(basevariables == which(names(MasterData) == NameofVariable))]}
  Formula2 <- FB_CreateFormula('AFLogical' , basevariables , MasterData )
  model <- glm(FB_CreateFormula('AFLogical' , unique(c(which(names(MasterData) == NameofVariable), which(names(MasterData) == 'ProcDetails') ,which(names(MasterData) == 'LogisticEUROScore'),which(names(MasterData) == 'CPB')  ) ) , MasterData ) ,
               family=binomial(link='logit'),
               data = MasterData)
  output[ , 8] <- exp(coef(model))[2]
  output[ , 9:10] <- exp(confint(model) )[2,]
  output[ , 11] <- summary(model)$coefficients[2,4]  
  output[ , 2:11] <- signif(output[,2:11] , 4)
  output[ , 12] <- -(FC_CalculateCrossValidatedROC(Formula1 , which(FC_ExtractIndiciesFromFormula(names(MasterData),Formula1) ==1) , MasterData = MasterData)-FC_CalculateCrossValidatedROC(Formula2 , which(FC_ExtractIndiciesFromFormula(names(MasterData),Formula1) ==1) , MasterData = MasterData))
  
  return(output)
}}


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
                                         "CardiogenicShock_pre_Operation",
                                         "IntravenousInotropesPriorToAnaesthesia",
                                         "IntravenousNitratesOrAnyHeparin",
                                         "ExtracardiacArteriopathy",
                                         "Planned.Valve.Surgery",
                                         "Pre_OperativeHeartRhythm"
  )
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
  NamesofPostopCategoricalVariables <- c('IABP2',
                                         "AKIUODayOne",
                                         "Filter")
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
                                        'PaO2OverFiO2',
                                        'Urine')
  NamesofLogisticEuroScoreVariables <- c( 'Age',
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
                                          "CardiogenicShock_pre_Operation",
                                          "IntravenousInotropesPriorToAnaesthesia",
                                          "IntravenousNitratesOrAnyHeparin",
                                          "ExtracardiacArteriopathy",
                                          "Planned.Valve.Surgery"
  )
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
                                  'DobutamineStandardised',
                                  'ReliableART.M',
                                  'PaO2OverFiO2',
                                  'Filter'
  )
  NamesofLogCasusScoreVariables <- c( 'FiO26num',
                                      'ArtPO2',
                                      'VentilatedInOpLogical',
                                      'Platelets',
                                      'Creatinine',
                                      'Bilirubin',
                                      'Lac',
                                      'IABP2',
                                      'PaO2OverFiO2',
                                      'GCSnum'
  )
  NamesPreOpBleedingVariables <- c('PreOpCRP',
                                   'PreopHb',
                                   'PreopPLT',
                                   'PreopWBC',
                                   'PreopPT',
                                   'PreopAPTT'
  )
  NamesPreOpCardioVascularVariables <- c('EjectionFractionCategory',
                                         'NYHAGrade',
                                         'AnginaGrade',
                                         'Active.Endocarditis',
                                         'HistoryOfNeurologicalDysfunction',
                                         'HistoryOfPulmonaryDisease',
                                         'Diabetes',
                                         'PreviousCardiacSurgery',
                                         'Recent.MI',
                                         'VAD',
                                         'CardiogenicShock_pre_Operation',
                                         'IntravenousInotropesPriorToAnaesthesia',
                                         'IntravenousNitratesOrAnyHeparin',
                                         'ExtracardiacArteriopathy'
  )
  NamesPreOpDemographicVariables <- c('Gender',
                                      'Age',
                                      'Weight',
                                      'BMI'
  )
  NamesPreOpElectrolytesVariables <- c('PreOpNa',
                                       'PreOpK',
                                       'PreopMg'
  )
  NamesPreOpInflamatoryVariables <- c('PreOpCRP' ,
                                      'PreOpAlb' ,
                                      'PreopHb' ,
                                      'PreopPLT' ,
                                      'PreopWBC' ,
                                      'PreopPT' ,
                                      'PreopAPTT'
  )
  NamesPreOpLiverVariables <- c('PreOpAlb',
                                'PreopBili'
  )
  NamesPreOpOperativeVariables <- c('Valve',
                                    'CABG',
                                    'Aortic',
                                    'Complex',
                                    'Urgency',
                                    'Thoracic.Aorta',
                                    "Planned.Valve.Surgery"
                                    
  )
  NamesPreOpRenalVariables <- c('PreopUrea',
                                'PreopCreat'
  )
  NamesPostOpBleedingVariables <- c('Hb',
                                    'Platelets',
                                    'PT',
                                    'APTT',
                                    'INR',
                                    'dHb',
                                    'dPLT',
                                    'dPT',
                                    'dAPTT'
  )
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
  NamesPostOpElectrolytesVariables <- c('Na',
                                        'K',
                                        'Mg',
                                        'dNa',
                                        'dK',
                                        'dMg'
                                        
  )
  NamesPostOpInflamatoryVariables <- c('Hb',
                                       'Platelets',
                                       'PT',
                                       'APTT',
                                       'INR',
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
  NamesPostOpLiverVariables <- c('Albumin',
                                 'Bilirubin',
                                 'dAlb',
                                 'dBili')
  NamesPostOpOperativeVariables <- c('CPB',
                                     'IABP2')
  
  NamesPostOpRenalVariables <- c('CVVHDyalysis',
                                 'Urea',
                                 'Creatinine',
                                 'dUrea',
                                 'dCreat')
  NamesPostOpRespiratoryVariables <- c('SpO2',
                                       'ArtPO2',
                                       'FiO26num',
                                       'PaO2OverFiO2')
}


#### Begin Data Analysis #####
source('ADDAModelSelection.R')


PreOpCategoricalUnivariateResultsLatex
PostOpCategoricalUnivariateResultsLatex
PreOpContinuousUnivariateResultsLatex
PostOpContinuousUnivariateResultsLatex

PostopBleedingModelTable
PostopCardioVascularModelTable
PostopElectrolytesModelTable 
PostopInflamatoryModelTable 
PostopLiverModelTable
PostopOperativeModelTable
PostopRenalModelTable
PostopRespiratoryModelTable
Postopfullmodeltable

