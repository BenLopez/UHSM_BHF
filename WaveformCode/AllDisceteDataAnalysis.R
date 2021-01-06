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

source('ADDAProcessData.R')
#### Begin Data Analysis #####

# Pre_op Catagorical Data Structure
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
                                       "Planned.Valve.Surgery"
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
NamesofPostopCategoricalVariables <- c('Filter',
                                       "AKIUODayOne")
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
                                'PaO2OverFiO2'
)
NamesofLogCasusScoreVariables <- c( 'FiO26num',
                                    'ArtPO2',
                                    'VentilatedInOpLogical',
                                    'Platelets',
                                    'Creatinine',
                                    'Bilirubin',
                                    'Filter',
                                    'Lac',
                                    'Filter',
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
                                  'Fibrinogen',
                                  'dHb',
                                  'dPLT',
                                  'dPT',
                                  'dAPTT'
)
NamesPostOpCardioVascularVariables <- c('HR',
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
NamesPostOpLiverVariables <- c('Albumin',
                               'Bilirubin',
                               'dAlb',
                               'dBili')
NamesPostOpOperativeVariables <- c('CPB')

NamesPostOpRenalVariables <- c('CVVHDyalysis',
                               'Urea',
                               'Creatinine',
                               'dUrea',
                               'dCreat')
NamesPostOpRespiratoryVariables <- c('SpO2',
                                     'ArtPO2',
                                     'FiO26num',
                                     'PaO2OverFiO2')


source('ADDAModelSelection.R')


PreOpCategoricalUnivariateResultsLatex
PostOpCategoricalUnivariateResultsLatex
PreOpContinuousUnivariateResultsLatex
PostOpContinuousUnivariateResultsLatex

LogisticEuroscoreModelTableLatex
SOFAEuroscoreModelTableLatex
CASUSEuroscoreModelTableLatex

PreopBleedingModelTable
PreopDemographicModelTable
PreopElectrolytesModelTable 
PreopInflamatoryModelTable 
PreopLiverModelTable
PreopOperativeModelTable
PreopRenalModelTable
PreopCardioVascularModelTable

PostopBleedingModelTable
PostopCardioVascularModelTable
PostopElectrolytesModelTable 
PostopInflamatoryModelTable 
PostopLiverModelTable
PostopOperativeModelTable
PostopRenalModelTable
PostopRespiratoryModelTable

Preopfullmodeltable
Postopfullmodeltable


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
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpCardioVascularComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpElectrolytesComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpInflamatoryComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpLiverComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpOperativeComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpRenalComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpRespiratoryComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(SOFAComparisonOutputs$ModelAUC, DataForLogistic))
PostOpModelOutputs <- cbind( PostOpModelOutputs ,  POM_CalulateProbabilitiesFromModel(PostOpStepOutput[[2]], DataForLogistic))


BC_PlotPairsFromTwoVariables(PostOpModelOutputs[AFLogical ==1,1:9] ,PostOpModelOutputs[sample(which(AFLogical ==0),260),1:9],
                             alpha = 0.1 ,
                             labels = c('Bleeding' , 'Cardiovascular' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory', 'SOFA', 'PostOpModel'))

PostOpC <- cov2cor(cov(PostOpModelOutputs))
rownames(PostOpC) <- c('Bleeding' , 'Cardio' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory', 'SOFA', 'PostOpModel')
colnames(PostOpC) <- c('Bleeding' , 'Cardio' , 'Electrolytes' , 'Inflamatory' , 'Liver', 'Operative' , 'Renal','Respiratory', 'SOFA', 'PostOpModel')

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

