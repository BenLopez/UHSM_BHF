if(!file.exists('MasterData.RData')){
  MasterData <- POM_ExtractFirstRecords( AllDataStructure )
save(MasterData , file = 'MasterData.RData')
  }else{
  load('MasterData.RData')
}

#### Process Data ####
# Procedure details 

MasterData <- MasterData[apply(as.matrix(MasterData$NewPseudoId) , 1 , DP_ExtractOpCodeFromNewPatinetID) == '-1-1-1' , ]

#MasterData$ProcDetails <- DP_AssignSurgeryLabels(MasterData$ProcDetails)
#MasterData$ProcDetails[MasterData$ProcDetails == 'MCS' | MasterData$ProcDetails == 'PP' | MasterData$ProcDetails == 'SP' | MasterData$ProcDetails == 'Transplant'] = 'Other'

# Create AFLogical

MasterData$ConfirmedFirstNewAF[is.na(MasterData$ConfirmedFirstNewAF)] <- 'NA'
MasterData$AFLogical = ((!is.na(MasterData$FirstNewAF)) & (MasterData$ConfirmedFirstNewAF != 'CNAF'))

# Remove treated with AFib medication

AFMedication <- POM_ExtractAFMedication(AllDataStructure = AllDataStructure)
MasterData <- MasterData[-which(  ((MasterData$AFLogical ==0)&(AFMedication$DidoxinLogical == 1)) | ((MasterData$AFLogical ==0)&(AFMedication$AmiodaroneLogical == 1)) ),]

# Remove Pre_op AFib

MasterData <- MasterData[MasterData$Pre_OperativeHeartRhythm != "Atrial fibrillation/flutter" , ]

# Process vasoactive drugs
{MasterData$NoradrenalineStandardised <- (MasterData$Noradrenaline..80mcg.ml.Volume/((1/80)*60*MasterData$Weight) )
  MasterData$NoradrenalineStandardised[is.na(MasterData$NoradrenalineStandardised )] <- (MasterData$Noradrenaline..160.mcg.ml.Volume[is.na(MasterData$NoradrenalineStandardised )]/((1/160)*60*MasterData$Weight[is.na(MasterData$NoradrenalineStandardised )]) )
  MasterData$NoradrenalineStandardised[is.na(MasterData$NoradrenalineStandardised )] <- (MasterData$Noradrenaline..320.mcg.ml.Volume[is.na(MasterData$NoradrenalineStandardised )]/((1/320)*60*MasterData$Weight[is.na(MasterData$NoradrenalineStandardised )]) )
  MasterData$NoradrenalineStandardised[is.na(MasterData$NoradrenalineStandardised )] <- (MasterData$Noradrenaline..unknown.dose.Volume[is.na(MasterData$NoradrenalineStandardised )]/((1/80)*60*MasterData$Weight[is.na(MasterData$NoradrenalineStandardised )]) )
  MasterData$NoradrenalineStandardised[is.na(MasterData$NoradrenalineStandardised)] <- 0
  
  MasterData$AdrenalineStandardised <- (MasterData$Adrenaline..20mcg.ml.Volume../((1/20)*60*MasterData$Weight) )
  MasterData$AdrenalineStandardised[is.na(MasterData$AdrenalineStandardised )] <- (MasterData$Adrenaline..40mcg.ml.Volume..[is.na(MasterData$AdrenalineStandardised )]/((1/40)*60*MasterData$Weight[is.na(MasterData$AdrenalineStandardised )]) )
  MasterData$AdrenalineStandardised[is.na(MasterData$AdrenalineStandardised )] <- (MasterData$Adrenaline..80mcg.ml.Volume..[is.na(MasterData$AdrenalineStandardised )]/((1/80)*60*MasterData$Weight[is.na(MasterData$AdrenalineStandardised )]) )
  MasterData$AdrenalineStandardised[is.na(MasterData$AdrenalineStandardised )] <- (MasterData$Adrenaline..unknown.dose.Volume[is.na(MasterData$AdrenalineStandardised )]/((1/20)*60*MasterData$Weight[is.na(MasterData$AdrenalineStandardised )]) )
  MasterData$AdrenalineStandardised[is.na(MasterData$AdrenalineStandardised)] <- 0
  
  MasterData$DopamineStandardised <- (MasterData$Dopamine.200mg.50mls..Volume/((1/4000)*60*MasterData$Weight) )
  MasterData$DopamineStandardised[is.na(MasterData$DopamineStandardised )] <- (MasterData$Dopamine.unknown.dose.Volume[is.na(MasterData$DopamineStandardised )]/((1/4000)*60*MasterData$Weight[is.na(MasterData$DopamineStandardised )]) )
  MasterData$DopamineStandardised[is.na(MasterData$DopamineStandardised)] <- 0
  
  MasterData$MilrinoneStandardised <- (MasterData$Milrinone.10mg.in.50mls.Volume*200)/(60*MasterData$Weight)
  MasterData$MilrinoneStandardised[is.na(MasterData$MilrinoneStandardised )] <- (MasterData$Milrinone.20mg.in.50mls.Volume[is.na(MasterData$MilrinoneStandardised )]*400)/(60*MasterData$Weight[is.na(MasterData$MilrinoneStandardised )])
  MasterData$MilrinoneStandardised[is.na(MasterData$MilrinoneStandardised )] <- (MasterData$Milrinone.unknown.dose.Rate[is.na(MasterData$MilrinoneStandardised )]*200)/(60*MasterData$Weight[is.na(MasterData$MilrinoneStandardised )])
  MasterData$MilrinoneStandardised[is.na(MasterData$MilrinoneStandardised)] <- 0

 # MasterData$VasopressinStandardised <- (MasterData$Vasopressin..mL...Volume*0.4)/(60*MasterData$Weight)
#  MasterData$VasopressinStandardised[is.na(MasterData$VasopressinStandardised)] <- 0

  MasterData$DobutamineStandardised <- (MasterData$Dobutamine.250mg.50mls.Volume*5000)/(60*MasterData$Weight)
  MasterData$DobutamineStandardised[is.na(MasterData$DobutamineStandardised)] <-  (MasterData$Dobutamine.unknown.dose.Volume*5000)/(60*MasterData$Weight)
  MasterData$DobutamineStandardised[is.na(MasterData$DobutamineStandardised)] <- 0
  
}

# Intubation in surgery 

MasterData$IntubatedInSurgury <- abs(DP_StripTime(MasterData$DateTubed) - DP_StripTime(MasterData$FirstITUEntry)) <=0

MasterData$GCSnum[is.na(MasterData$GCSnum) & (MasterData$IntubatedInSurgury == 0)] <- 15

# Fix strings which should be numeric

ListNumericVariables <- c('Age' ,
                          'Weight' , 
                          "LogisticEUROScore" ,
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
                          "PreopAPTT",
                          'CPB',
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
                          'PostOpFluids',
                          'HR',
                          "ReliableART.S",
                          "ReliableART.M",
                          "ReliableART.D",
                          "CVP",
                          "SpO2",
                          'NoradrenalineStandardised',
                          'DopamineStandardised',
                          'AdrenalineStandardised',
                          'MilrinoneStandardised',
                          'DobutamineStandardised',
                          'ArtPO2',
                          'Lac',
                          'FiO26num',
                          'GCSnum')

MasterData$CPB[MasterData$CPB == 'not recorded'] = NA
MasterData$CPB[is.na(MasterData$CPB)] <- "00:00"
MasterData$CPB <- as.numeric(difftime(strptime(x = MasterData$CPB ,  format = '%R') ,  sort(strptime(x = MasterData$CPB ,  format = '%R'))[1],units = c('mins') ))

MasterData[ ,names(MasterData) %in% ListNumericVariables ] <- apply(MasterData[ ,names(MasterData) %in% ListNumericVariables ] , 2 , as.numeric )

# Make Factors

ListofCategoricalVariables  <- c('AFlogical',
                                 'Gender' ,
                                 'ProcDetails' ,
                                 "EjectionFractionCategory",
                                 "NYHAGrade",
                                 "AnginaGrade",
                                 "Filter",
                                 "IABP2",
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

MasterData[ ,names(MasterData) %in% ListofCategoricalVariables ] <- apply(MasterData[ ,names(MasterData) %in% ListofCategoricalVariables ] , 2 , as.factor )



# Outliers

MasterData <- POM_RemoveOutliers(MasterPreOpData = MasterData)

# Calulated variables

{ # Differences
  MasterData$dMg <- MasterData$Mg - MasterData$PreopMg
  MasterData$dAlb <- MasterData$Albumin - MasterData$PreOpAlb
  MasterData$dUrea <- MasterData$Urea - MasterData$PreopUrea
  MasterData$dNa <- MasterData$Na - MasterData$PreOpNa
  MasterData$dBili <- MasterData$Bilirubin - MasterData$PreopBili
  MasterData$dCreat <- MasterData$Creatinine - MasterData$PreopCreat
  MasterData$dK <-   MasterData$K - MasterData$PreOpK
  MasterData$dPT <-  MasterData$PT - MasterData$PreopPT
  MasterData$dCRP <- MasterData$CRP - MasterData$PreOpCRP 
  MasterData$dPLT <- MasterData$Platelets - MasterData$PreopPLT
  MasterData$dHb <- MasterData$Hb - MasterData$PreopHb 
  MasterData$dAPTT <- MasterData$APTT - MasterData$PreopAPTT
  MasterData$dWBC <- MasterData$WBC - MasterData$PreopWBC
  
  MasterData$BMI = MasterData$Weight/(MasterData$Height/100)
  MasterData$PaO2OverFiO2 <- MasterData$ArtPO2/MasterData$FiO26num
  MasterData$PAR <- (MasterData$HR*MasterData$CVP)/MasterData$ReliableART.M
  
}
# Post-op risk scores

PostOpRiskScore <- POM_ExtraPostOpScores(PostOpScoresIndex2017 = PostOpScoresIndex2017)

MasterData <- data.frame(MasterData , logCASUS = matrix(NA , dim(MasterData)[1]) ,SOFA = matrix(NA , dim(MasterData)[1]) ,RACE = matrix(NA , dim(MasterData)[1])  )
MasterData <- MasterData[!is.na(MasterData$PseudoId),] 

for(i in 1:dim(MasterData)[1]){
  if(!is.na(MasterData$PseudoId[i])){
    if(sum(MasterData$PseudoId[i] == rownames(PostOpRiskScore))>0){
      MasterData$logCASUS[i] = PostOpRiskScore[which(MasterData$PseudoId[i] == rownames(PostOpRiskScore)) , 1]
      MasterData$SOFA[i] = PostOpRiskScore[which(MasterData$PseudoId[i] == rownames(PostOpRiskScore)) , 3]
      MasterData$RACE[i] = PostOpRiskScore[which(MasterData$PseudoId[i] == rownames(PostOpRiskScore)) , 2]
    }
  }
}

# Transformations

MasterData$LogisticEUROScore <- log(MasterData$LogisticEUROScore)
MasterData$SCTSLogisticEuroSCORE <- MasterData$SCTSLogisticEuroSCORE

MasterData$CRP <- log(MasterData$CRP)
MasterData$PreOpCRP <- log(MasterData$PreOpCRP)

MasterData$RACE <- log(MasterData$RACE)
MasterData$logCASUS <- log(MasterData$logCASUS)

MasterData$PreopBili <- log(MasterData$PreopBili)
MasterData$Bilirubin <- log(MasterData$Bilirubin)


# Categorical data fixes for odd fields
{MasterData$EjectionFractionCategory[MasterData$EjectionFractionCategory == ''] <- NA
MasterData$EjectionFractionCategory[MasterData$EjectionFractionCategory == "Not measured"] <- NA
MasterData$EjectionFractionCategory[MasterData$EjectionFractionCategory == "Good (LVEF > 50%)"] <- "LVEF > 50%"
MasterData$EjectionFractionCategory[MasterData$EjectionFractionCategory == "Fair (LVEF 30-50%)"] <- "LVEF  30-50%"
MasterData$EjectionFractionCategory[MasterData$EjectionFractionCategory == "Poor (LVEF < 30%)"] <- "LVEF < 30%"


MasterData$NYHAGrade[MasterData$NYHAGrade == ''] <- NA
MasterData$NYHAGrade[MasterData$NYHAGrade == 'Not Known'] <- NA
MasterData$NYHAGrade[MasterData$NYHAGrade == "Slight Limitation (II)"] <- 'NYHA(O-II)'
MasterData$NYHAGrade[MasterData$NYHAGrade == 'Not Limiting (I)'] <- 'NYHA(O-II)'
MasterData$NYHAGrade[MasterData$NYHAGrade == "Marked Limitation (III)"] <- 'NYHA(III-IV)'
MasterData$NYHAGrade[MasterData$NYHAGrade == "Any Activity or At Rest (IV)"] <- 'NYHA(III-IV)'
MasterData$NYHAGrade[MasterData$NYHAGrade == "None"] <- 'NYHA(O-II)'

MasterData$AnginaGrade[MasterData$AnginaGrade == ''] <- NA
MasterData$AnginaGrade[MasterData$AnginaGrade == 'Unknown'] <- NA
MasterData$AnginaGrade[MasterData$AnginaGrade == "Strenuous Exertion (I)"   ] <- 'AnginaGrade(O-II)'
MasterData$AnginaGrade[MasterData$AnginaGrade == "Slight Limitation (II)"     ] <- 'AnginaGrade(O-II)'
MasterData$AnginaGrade[MasterData$AnginaGrade == "Marked Limitation (III)"     ] <- 'AnginaGrade(III-IV)'
MasterData$AnginaGrade[MasterData$AnginaGrade == "Any Activity or At Rest (IV)"     ] <- 'AnginaGrade(IV-IV)'
MasterData$AnginaGrade[MasterData$AnginaGrade == "No"     ] <- 'AnginaGrade(O-II)'


#MasterData$IntubatedInSurgury[MasterData$IntubatedInSurgury == ' TRUE'] <- 'IntupatedInSurgery'
#MasterData$IntubatedInSurgury[MasterData$IntubatedInSurgury == 'FALSE'] <- 'NotIntupatedInSurgery'

MasterData$Filter[MasterData$Filter == ' TRUE'] <- 'CVVHDyalysis(Yes)'
MasterData$Filter[MasterData$Filter == 'FALSE'] <- 'CVVHDyalysis(No)'

MasterData$IABP2[MasterData$IABP2 == 'No'] <- 'IABP2(No)'
MasterData$IABP2[MasterData$IABP2 == 'Yes'] <- 'IABP2(Yes)'


MasterData$Active.Endocarditis[MasterData$Active.Endocarditis == ''] <- NA
MasterData$Active.Endocarditis[MasterData$Active.Endocarditis == 'No'] <- 'Active.Endocarditis(No)'
MasterData$Active.Endocarditis[MasterData$Active.Endocarditis == 'Yes'] <- 'Active.Endocarditis(Yes)'

MasterData$HTN[MasterData$HTN == ''] <- NA
MasterData$HTN[MasterData$HTN == 'Unknown'] <- NA
MasterData$HTN[MasterData$HTN== 'No'] <- 'HTN(No)'
MasterData$HTN[MasterData$HTN == 'Yes'] <- 'HTN(Yes)'


MasterData$HistoryOfNeurologicalDysfunction[MasterData$HistoryOfNeurologicalDysfunction == ''] <- NA
MasterData$HistoryOfNeurologicalDysfunction[MasterData$HistoryOfNeurologicalDysfunction== 'No'] <- 'HistoryOfNeurologicalDysfunction(No)'
MasterData$HistoryOfNeurologicalDysfunction[MasterData$HistoryOfNeurologicalDysfunction == 'Yes'] <- 'HistoryOfNeurologicalDysfunction(Yes)'

MasterData$HistoryOfPulmonaryDisease[MasterData$HistoryOfPulmonaryDisease == ''] <- NA
MasterData$HistoryOfPulmonaryDisease[MasterData$HistoryOfPulmonaryDisease== 'No'] <- 'HistoryOfPulmonaryDisease(No)'
MasterData$HistoryOfPulmonaryDisease[MasterData$HistoryOfPulmonaryDisease== 'NO '] <- 'HistoryOfPulmonaryDisease(No)'
MasterData$HistoryOfPulmonaryDisease[MasterData$HistoryOfPulmonaryDisease== 'No '] <- 'HistoryOfPulmonaryDisease(No)'
MasterData$HistoryOfPulmonaryDisease[MasterData$HistoryOfPulmonaryDisease == 'Yes'] <- 'HistoryOfPulmonaryDisease(Yes)'

MasterData$Diabetes[MasterData$Diabetes == ''] <- NA
MasterData$Diabetes[MasterData$Diabetes == 'Unknown'] <- NA
MasterData$Diabetes[MasterData$Diabetes== 'No'] <- 'Diabetes(No)'
MasterData$Diabetes[MasterData$Diabetes == 'Yes'] <- 'Diabetes(Yes)'

MasterData$Thoracic.Aorta[MasterData$Thoracic.Aorta == ''] <- NA
MasterData$Thoracic.Aorta[MasterData$Thoracic.Aorta== 'No'] <- 'Thoracic.Aorta(No)'
MasterData$Thoracic.Aorta[MasterData$Thoracic.Aorta == 'Yes'] <- 'Thoracic.Aorta(Yes)'

MasterData$PreviousCardiacSurgery[MasterData$PreviousCardiacSurgery != ''] <- 'PreviousCardiacSurgery(Yes)'
MasterData$PreviousCardiacSurgery[MasterData$PreviousCardiacSurgery == ''] <- 'PreviousCardiacSurgery(No)'

MasterData$Recent.MI[MasterData$Recent.MI == ''] <- NA
MasterData$Recent.MI[MasterData$Recent.MI == 'No'] <- 'Recent.MI(No)'
MasterData$Recent.MI[MasterData$Recent.MI == 'No '] <- 'Recent.MI(No)'

MasterData$Recent.MI[MasterData$Recent.MI != 'Recent.MI(No)' & !is.na(MasterData$Recent.MI)] <- 'Recent.MI(Yes)'

MasterData$PreOpSupport[MasterData$PreOpSupport == ''] <- NA
MasterData$PreOpSupport[which(MasterData$PreOpSupport != 'No ')] <- 'PreOpSupport(Yes)'
MasterData$PreOpSupport[MasterData$PreOpSupport== 'No '] <- 'PreOpSupport(No)'

MasterData$VentilatedPreOperation[MasterData$VentilatedPreOperation == ''] <- NA
MasterData$VentilatedPreOperation[MasterData$VentilatedPreOperation== 'No'] <- 'VentilatedPreOperation(No)'
MasterData$VentilatedPreOperation[MasterData$VentilatedPreOperation == 'Yes'] <- 'VentilatedPreOperation(Yes)'


MasterData$CardiogenicShock_pre_Operation[MasterData$CardiogenicShock_pre_Operation == ''] <- NA
MasterData$CardiogenicShock_pre_Operation[MasterData$CardiogenicShock_pre_Operation== 'No'] <- 'CardiogenicShock_pre_Operation(No)'
MasterData$CardiogenicShock_pre_Operation[MasterData$CardiogenicShock_pre_Operation == 'Yes'] <- 'CardiogenicShock_pre_Operation(Yes)'


MasterData$IntravenousInotropesPriorToAnaesthesia[MasterData$IntravenousInotropesPriorToAnaesthesia == ''] <- NA
MasterData$IntravenousInotropesPriorToAnaesthesia[MasterData$IntravenousInotropesPriorToAnaesthesia== 'No'] <- 'IntravenousInotropesPriorToAnaesthesia(No)'
MasterData$IntravenousInotropesPriorToAnaesthesia[MasterData$IntravenousInotropesPriorToAnaesthesia == 'Yes'] <- 'IntravenousInotropesPriorToAnaesthesia(Yes)'


MasterData$IntravenousNitratesOrAnyHeparin[MasterData$IntravenousNitratesOrAnyHeparin == ''] <- NA
MasterData$IntravenousNitratesOrAnyHeparin[MasterData$IntravenousNitratesOrAnyHeparin== 'No'] <- 'IntravenousNitratesOrAnyHeparin(No)'
MasterData$IntravenousNitratesOrAnyHeparin[MasterData$IntravenousNitratesOrAnyHeparin == 'Yes'] <- 'IntravenousNitratesOrAnyHeparin(Yes)'


MasterData$ExtracardiacArteriopathy[MasterData$ExtracardiacArteriopathy == ''] <- NA
MasterData$ExtracardiacArteriopathy[MasterData$ExtracardiacArteriopathy== 'No'] <- 'ExtracardiacArteriopathy(No)'
MasterData$ExtracardiacArteriopathy[MasterData$ExtracardiacArteriopathy == 'Yes'] <- 'ExtracardiacArteriopathy(Yes)'


MasterData$Planned.Valve.Surgery[(MasterData$Planned.Valve.Surgery == '') ] <- 'None'
MasterData$Planned.Valve.Surgery[!(MasterData$Planned.Valve.Surgery == 'Mitral Valve') & !(MasterData$Planned.Valve.Surgery == "Aortic Valve") & !(MasterData$Planned.Valve.Surgery == "Tricuspid Valve") & !(MasterData$Planned.Valve.Surgery =='None') ] <- 'Multiple'

MasterData$ProcDetails[(MasterData$ProcDetails == 'CABG') & (MasterData$Planned.Valve.Surgery != 'None')] = "CABG and Valve"

MasterData <- cbind(MasterData , DP_AssignSurgeryLabels2(MasterData$ProcDetails , ProcedureDetails ))

MasterData$Complex <- as.vector(MasterData$Complex)
MasterData$Complex[MasterData$Complex == 'Y'] = 'Procedure Details(Complex)'
MasterData$Valve <- as.vector(MasterData$Valve)
MasterData$Valve[MasterData$Valve == 'Y'] = 'Procedure Details(Valve)'
MasterData$CABG <- as.vector(MasterData$CABG)
MasterData$CABG[MasterData$CABG == 'Y'] = 'Procedure Details(CABG)'
MasterData$Aortic <- as.vector(MasterData$Aortic)
MasterData$Aortic[MasterData$Aortic == 'Y'] = 'Procedure Details(Aortic)'
MasterData$Other <- as.vector(MasterData$Other)
MasterData$Other[MasterData$Other == 'Y'] = 'Procedure Details(Other)'

MasterData[MasterData$Other == 'N',]

}

ListofCategoricalVariables  <- c('AFlogical',
                                 'Gender' ,
                                 'ProcDetails' ,
                                 'Valve',
                                 'CABG',
                                 'Aortic',
                                 'Complex',
                                 'Other',
                                 "EjectionFractionCategory",
                                 "NYHAGrade",
                                 "AnginaGrade",
                                 "Filter",
                                 "IABP2",
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

MasterData[ ,names(MasterData) %in% ListofCategoricalVariables ] <- apply(MasterData[ ,names(MasterData) %in% ListofCategoricalVariables ] , 2 , as.factor )

MasterData$Planned.Valve.Surgery[MasterData$Valve == "Procedure Details(Valve)" & MasterData$Planned.Valve.Surgery == 'None'] <-  "Aortic Valve" 

# difference transforms
{
  MasterData$dUrea <- abs(DP_NormaliseData(MasterData$dUrea) )
  MasterData$dK <-abs(DP_NormaliseData(MasterData$dK))
  MasterData$dCreat <- abs(DP_NormaliseData(MasterData$dCreat) )
  MasterData$dNa <- abs( DP_NormaliseData(MasterData$dNa) )
  MasterData$dCRP <- abs( DP_NormaliseData(abs(DP_NormaliseData(MasterData$dCRP) )))
  MasterData$dWBC <- abs( DP_NormaliseData(MasterData$dWBC) )
  MasterData$dAlb <- abs( DP_NormaliseData(MasterData$dAlb) )
}

AFLogical <- MasterData$AFLogical
MasterData$CICUStay <- abs(DP_StripTime(MasterData$AdmitDateTime) - DP_StripTime( MasterData$DisDateTime) )

MasterData$Valve[MasterData$Planned.Valve.Surgery != 'None'] <- "Procedure Details(Valve)"