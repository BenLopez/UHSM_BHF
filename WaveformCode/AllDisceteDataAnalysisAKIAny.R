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

source('ADDAProcessData2.R')
#MasterData$AFLogical <- MasterData$AFLogical[sample(1:dim(MasterData)[1],dim(MasterData)[1])]
#AFLogical <- MasterData$AFLogical
{
NamesofPeropCategoricalVariables <- c( 'Gender' ,
                                       'ProcDetails',
                                       'ProcDetailsReduced',
                                       "EjectionFractionCategory",
                                       "NYHAGrade" ,
                                       "AnginaGrade",
                                       'Urgency',
                                       'UrgencyReduced',
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
                                      #'dNa',
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
                                        'ProcDetails',
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
NamesPreOpOperativeVariables <- c('ProcDetails',
                                  'ProcDetailsReduced',
                                  'UrgencyReduced',
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
source( 'ADDAModelSelection.R' )
  
#### Tables ####

tmp <- gsub('AFib' ,'AKI' , colnames(PreOpCategoricalUnivariateResultsLatex))
tmp[12] <- 'Diff CV AUC'
colnames(PreOpCategoricalUnivariateResultsLatex) <- tmp
tmp <- gsub('AFib' ,'AKI' , colnames(PreOpContinuousUnivariateResultsLatex))
tmp[12] <- 'Diff AUC'
colnames(PreOpContinuousUnivariateResultsLatex) <- tmp
tmp <- gsub('AFib' ,'AKI' , colnames(PostOpCategoricalUnivariateResultsLatex))
tmp[12] <- 'Diff AUC'
colnames(PostOpCategoricalUnivariateResultsLatex) <- tmp
tmp <- gsub('AFib' ,'AKI' , colnames(PostOpContinuousUnivariateResultsLatex ))
tmp[12] <- 'Diff AUC'
colnames(PostOpContinuousUnivariateResultsLatex ) <- tmp

PostOpCategoricalUnivariateResultsLatex
PreOpCategoricalUnivariateResultsLatex

PreOpContinuousUnivariateResultsLatex
PostOpContinuousUnivariateResultsLatex

PostopBleedingModelTable
PostopCardioVascularModelTable
PostopElectrolytesModelTable 
PostopInflamatoryModelTable 
PostopLiverModelTable
PreopOperativeModelTable
PostopRenalModelTable
PostopRespiratoryModelTable
PreopDemographicModelTable

Postopfullmodeltable