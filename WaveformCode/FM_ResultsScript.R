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

# Corrections after consultantion with experts and outputs of History Matching
{
PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z209'] = "07/10/2016 10:46"
PatIndex2017$EndFirstNewAF[PatIndex2017$PseudoId == 'z209'] = "07/10/2016 11:29"

PatIndex2017$FirstNewAF[PatIndex2017$PseudoId == 'z160'] = NA

PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z925'] = "10/06/2017 08:31"
PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z401'] = "04/12/2016 19:27"
PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z1281'] = "26/10/2017 01:04"
PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z580'] = "04/02/2017 06:10"
}

{
VectorofDifferences <- matrix(NA , length(listAllPatients) , 1)
AFLogical <- matrix(NA , length(listAllPatients))
duds <- matrix(FALSE , length(listAllPatients))

counter <- 1
for(PatientID in listAllPatients[1:length(listAllPatients)] ){
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID )
  
AFLogical[counter] <- DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID ))

if( file.exists(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData')) ){ 

load(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData'))
if(length(outputstruct)==2){
  duds[counter] <- TRUE
  counter <- counter +1
  next}
t <- FM_ExtractTimeFromHMOutput(outputstruct ) 
AFAnnotationHM <- FM_CreateAFAnnoation( outputstruct[[length(outputstruct)]] )[1:(length(outputstruct)-2)]
AFAnnotationExpert <- BC_CreateAFAnnotationFomMetaData(t , MetaData = DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientID))

if(sum(AFAnnotationExpert)>0){
VectorofDifferences[counter] <- difftime(t[which(AFAnnotationHM > 0)[1]] ,DP_StripTime(MetaData$ConfirmedFirstNewAF), units = c('hours'))
}else{
VectorofDifferences[counter] <- 0   
}
if(sum(AFAnnotationExpert)==0 & sum(AFAnnotationHM)>0){
  VectorofDifferences[counter] <- -Inf   
}
if(sum(AFAnnotationExpert)>0 & sum(AFAnnotationHM)==0){
  VectorofDifferences[counter] <-  Inf   
}
}
DP_WaitBar(counter/length(listAllPatients))
counter <- counter + 1
}

p_results <- ggplot(data.frame(x = c(1:sum(AFLogical& !is.infinite(VectorofDifferences),na.rm = T)) ,y=c(VectorofDifferences[AFLogical == T& !is.infinite(VectorofDifferences)])),aes(x,y)) + 
  geom_point(col = 'blue') +
  ggtitle('Comparison of History Matching and Expert Annotation First Warning') +xlab('Patinet Index') + ylab('Difference in Annotation Time (Hours)')
}

# Patinet wise sensitivity and specificity (when the algroihtm say they ever go into AF, so does Sam)
{ 
P <- sum( AFLogical )
N <- sum( !AFLogical )
TP <- sum( !is.infinite( VectorofDifferences[AFLogical] ) )
TN <- sum( !is.infinite( VectorofDifferences[!AFLogical] ) )
FP <- sum( is.infinite( VectorofDifferences[!AFLogical] ) )
FN <- sum( is.infinite( VectorofDifferences[AFLogical] ) )  
  
Sensitivity <- TP/P
Specificity <- TN/N
PPV <- TP/ (TP + FP)
NPV <- TN / (TN + FN)

ResultsStructure <- matrix(0,4 , 4)
colnames(ResultsStructure) <- c('Sensitivty','Specificity','PPV' , 'NPV')
rownames(ResultsStructure) <- c('FM Patient wise','a','b' , 'c')

ResultsStructure[1,] <- round(c(Sensitivity , Specificity , PPV , NPV),3)

listAllPatients[!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences > 0.25]
listAllPatients[!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences< -1]
listAllPatients[is.infinite(VectorofDifferences) & AFLogical]
listAllPatients[is.infinite(VectorofDifferences) & !AFLogical]
}

# If 10 minutes early or late class as a fail.

TP <- sum( !is.infinite(VectorofDifferences[AFLogical] ) ) -
  (sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences > 1) + sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences< -1) )
TN <- sum( !is.infinite(VectorofDifferences[!AFLogical] ) )
FP <- sum( is.infinite(VectorofDifferences[!AFLogical] ) )
FN <- sum( is.infinite(VectorofDifferences[AFLogical] ) )  +
  (sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences > 1) + sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences< -1) )

Sensitivity <- TP/P
Specificity <- TN/N
PPV <- TP/ (TP + FP)
NPV <- TN / (TN + FN)

colnames(ResultsStructure) <- c('Sensitivty','Specificity','PPV' , 'NPV')
rownames(ResultsStructure) <- c('FM Patient wise','FM Within 15','b' , 'c')

ResultsStructure[2,] <- round(c(Sensitivity , Specificity , PPV , NPV),3)


##### Nurse Analysis #####

{ VectorofDifferences <- matrix(NA , length(listAllPatients) , 1)
  AFLogical <- matrix(NA , length(listAllPatients))
  duds <- matrix(FALSE , length(listAllPatients))

  counter <- 1
for(PatientID in listAllPatients[1:length(listAllPatients)] ){
  
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID )
  AFLogical[counter] <- DP_CheckIfAFPatient(MetaData)
  
  if(!is.na(MetaData$FirstNewAF) & !is.na(MetaData$ConfirmedFirstNewAF)){
  VectorofDifferences[counter] <- difftime(DP_StripTime(MetaData$FirstNewAF) , DP_StripTime(MetaData$ConfirmedFirstNewAF) ,units = 'hours')
  }else{
  VectorofDifferences[counter] <- 0
  }
  if(is.na(MetaData$FirstNewAF) & !is.na(MetaData$ConfirmedFirstNewAF)){
  VectorofDifferences[counter] <- -Inf
  }
  if(!is.na(MetaData$FirstNewAF) & is.na(MetaData$ConfirmedFirstNewAF)){
    VectorofDifferences[counter] <- Inf
  }
  DP_WaitBar(counter/length(listAllPatients))
  counter <- counter + 1
}

p_resultsNurse <- ggplot(data.frame(x = c(1:sum(AFLogical& !is.infinite(VectorofDifferences),na.rm = T)) ,y=c(VectorofDifferences[AFLogical == T& !is.infinite(VectorofDifferences)])),aes(x,y)) + 
  geom_point(col = 'blue') +
  ggtitle('Comparison of TCs and Expert Annotation First Warning') +xlab('Patinet Index') + ylab('Difference in Annotation Time (Hours)')
}

TP <- sum( !is.infinite(VectorofDifferences[AFLogical] ) )
TN <- sum( !is.infinite(VectorofDifferences[!AFLogical] ) )
FP <- sum( is.infinite(VectorofDifferences[!AFLogical] ) )
FN <- sum( is.infinite(VectorofDifferences[AFLogical] ) )  

Sensitivity <- TP/P
Specificity <- TN/N
PPV <- TP/ (TP + FP)
NPV <- TN / (TN + FN)

colnames(ResultsStructure) <- c('Sensitivty','Specificity','PPV' , 'NPV')
rownames(ResultsStructure) <- c('FM Patient wise','FM Within 15','CT Patient wise' , 'c')

ResultsStructure[3,] <- round(c(Sensitivity , Specificity , PPV , NPV),3)

TP <- sum( !is.infinite(VectorofDifferences[AFLogical] ) ) -
  (sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences > 1) + sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences< -1) )
TN <- sum( !is.infinite(VectorofDifferences[!AFLogical] ) )
FP <- sum( is.infinite(VectorofDifferences[!AFLogical] ) )
FN <- sum( is.infinite(VectorofDifferences[AFLogical] ) )  +
  (sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences > 1) + sum(!is.infinite(VectorofDifferences) & !is.na(VectorofDifferences) & VectorofDifferences< -1) )

Sensitivity <- TP/P
Specificity <- TN/N
PPV <- TP/ (TP + FP)
NPV <- TN / (TN + FN)

colnames(ResultsStructure) <- c('Sensitivty','Specificity','PPV' , 'NPV')
rownames(ResultsStructure) <- c('FM Patient wise','FM Within 15','CT Patient wise' , 'CT Patient Within 15')

ResultsStructure[4,] <- round(c(Sensitivity , Specificity , PPV , NPV),3)


##### Closest AF eopside to Sam's #####

{VectorofDifferences <- matrix(NA , length(listAllPatients) , 1)
AFLogical <- matrix(NA , length(listAllPatients))
duds <- matrix(FALSE , length(listAllPatients))

counter <- 1
for(PatientID in listAllPatients[1:length(listAllPatients)] ){
  
  AFLogical[counter] <- DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID ))
  MetaData = DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientID)
  if( file.exists(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData')) ){ 
    
    load(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData'))
    if(length(outputstruct)==2){
      duds[counter] <- TRUE
      counter <- counter +1
      next}
    t <- FM_ExtractTimeFromHMOutput(outputstruct ) 
    AFAnnotationHM <- FM_CreateAFAnnoation( outputstruct[[length(outputstruct)]] )[1:(length(outputstruct)-2)]
    AFAnnotationExpert <- BC_CreateAFAnnotationFomMetaData(t , MetaData = DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientID))
    
    if(sum(AFAnnotationExpert)>0){
      if(sum(AFAnnotationHM > 0 )){
      tmp <- difftime(t[which(AFAnnotationHM > 0)] ,DP_StripTime(MetaData$ConfirmedFirstNewAF[1]), units = c('mins'))
      VectorofDifferences[counter] <- tmp[which.min(abs(tmp))]}
    }else{
      VectorofDifferences[counter] <- 0   
    }
    if(sum(AFAnnotationExpert)==0 & sum(AFAnnotationHM)>0){
      VectorofDifferences[counter] <- -Inf   
    }
    if(sum(AFAnnotationExpert)>0 & sum(AFAnnotationHM)==0){
      VectorofDifferences[counter] <-  Inf   
    }
  }
  DP_WaitBar(counter/length(listAllPatients))
  counter <- counter + 1
}
p_results2 <- ggplot(data.frame(x = c(1:sum(AFLogical& !is.infinite(VectorofDifferences),na.rm = T)) ,y=c(VectorofDifferences[AFLogical == T& !is.infinite(VectorofDifferences)])),aes(x,y)) + 
  geom_point(col = 'blue') +
  ggtitle('Comparison of History Matching and Expert Annotation Closest Warning') +xlab('Patinet Index') + ylab('Difference in Annotation Time (Minutes)')
print(mean(VectorofDifferences[!is.infinite(VectorofDifferences)] , na.rm = T ))
print(var(VectorofDifferences[!is.infinite(VectorofDifferences)] , na.rm = T ))
}


#### CT closest time #####

{
VectorofDifferences <- matrix(NA , length(listAllPatients) , 1)
AFLogical <- matrix(NA , length(listAllPatients))
duds <- matrix(FALSE , length(listAllPatients))
counter <- 1
for(PatientID in listAllPatients[1:length(listAllPatients)] ){
  
  AFLogical[counter] <- DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID ))
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID )
  
  if( file.exists(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData')) ){ 
    
    load(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData'))
    if(length(outputstruct)==2){
      duds[counter] <- TRUE
      counter <- counter +1
      next}
    t <- FM_ExtractTimeFromHMOutput(outputstruct ) 
    AFAnnotationExpert <- BC_CreateAFAnnotationFomMetaData(t , MetaData = DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientID))
    
    if(DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID )) ){
      if(AFLogical[counter] >0 & !is.na(MetaData$FirstNewAF) ){
        tmp <- AllDataStructure[[which(names(AllDataStructure) == paste0(PatientID , '-1-1-1' ))]]$timeseriesvariables$time[AllDataStructure[[which(names(AllDataStructure) == paste0(PatientID , '-1-1-1' ))]]$timeseriesvariables$Rhythm == "AF or Flutter" & !is.na(AllDataStructure[[which(names(AllDataStructure) == paste0(PatientID , '-1-1-1' ))]]$timeseriesvariables$Rhythm)]
        tmp <- difftime(tmp ,DP_StripTime(MetaData$ConfirmedFirstNewAF) , units = c('mins'))
        VectorofDifferences[counter] <- tmp[which.min(abs(tmp))]
        }
    }else{
      VectorofDifferences[counter] <- 0   
    }
    if(AFLogical[counter] ==0 & !is.na(MetaData$FirstNewAF)){
      VectorofDifferences[counter] <- -Inf   
    }
    if(AFLogical[counter] >0 & is.na(MetaData$FirstNewAF)){
      VectorofDifferences[counter] <-  Inf   
    }
  }
  DP_WaitBar(counter/length(listAllPatients))
  counter <- counter + 1
}

p_results2Nurse <- ggplot(data.frame(x = c(1:sum(AFLogical & !is.infinite(VectorofDifferences),na.rm = T)) ,y=c(VectorofDifferences[AFLogical == T& !is.infinite(VectorofDifferences)])),aes(x,y)) + 
  geom_point(col = 'blue') +
  ggtitle('Comparison of TCs and Expert Annotation Closest Warning') +xlab('Patinet Index') + ylab('Difference in Annotation Time (Minutes)')
print(mean(VectorofDifferences[!is.infinite(VectorofDifferences)] , na.rm = T ))
print(var(VectorofDifferences[!is.infinite(VectorofDifferences)] , na.rm = T ))

}


x11(20,14)
grid.arrange(p_results,p_resultsNurse , p_results2, p_results2Nurse,ncol = 2, nrow=2)

xtable(ResultsStructure,digits = 3)
