{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

DP_LoadPatientIndex()
HowtoFilterops <- read.csv(choose.files(caption = "Select listofopsSH") , stringsAsFactors = FALSE)

DP_ChooseDataReps()
listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess = 'ECGI')
Patientsprocessedlogical <- lapply(listAllPatients , function(X){DP_checkRpeaksfilesprocessed(path , X ) } )
listAllPatients <- listAllPatients[which(Patientsprocessedlogical == TRUE)]

ListConfirmedNoAF <- c( 'z380' , 'z950' , 'z1086' , 'z1163' , 'z1203' , 'z827' )

RPeaksData <- list()
SummaryStatistics <- matrix(0 , length(listAllPatients) , 2)
AFPeriods <- list()
MMPeriods <- list()
PatientNames <- list()
counter <- 1
for(ii in 1:length(listAllPatients)){
  
   outputdata <- DP_LoadRpeaksfile(path , listAllPatients[ii])
  
if( (length(outputdata$RRCombined$t) > 10000) && (outputdata$Meta_Data$TotalITUTimeHRS < 100) && DP_CheckECGreducedfilesprocessed( path , listAllPatients[ii] , 'ECGI_reduced') && DP_CheckECGreducedfilesprocessed( path , listAllPatients[ii] , 'ECGII_reduced') && DP_CheckECGreducedfilesprocessed( path , listAllPatients[ii] , 'ECGIII_reduced') )
  {
  InferenceOutput <- AFD_DetectionWrapper( outputdata$RRCombined  )
  ECGs <- DP_LoadReducedECGs(path , listAllPatients[[ii]] , numberrep  , FilestoProcess = c('ECGI' , 'ECGII' , 'ECGIII'))
  if(length(InferenceOutput$StartEndTimesAF$Start) >0 ){
    InferenceOutput$StartEndTimesAF <- AFD_Checkformissingdata(StartEndTimesAF = InferenceOutput$StartEndTimesAF , AFScore = InferenceOutput$AFScore , ECGI = ECGs$ECGI , ECGII =  ECGs$ECGII , ECGIII =  ECGs$ECGIII)}else{
    InferenceOutput$StartEndTimesAF <- InferenceOutput$StartEndTimesAF
    }
  if(length(InferenceOutput$StartEndTimesMM$Start) >0 ){
    InferenceOutput$StartEndTimesMM <- AFD_Checkformissingdata(InferenceOutput$StartEndTimesMM , InferenceOutput$AFScore , ECGI = ECGs$ECGI , ECGII =  ECGs$ECGI , ECGIII =  ECGs$ECGIII)}else{
    InferenceOutput$StartEndTimesMM <- InferenceOutput$StartEndTimesMM
    }
  RPeaksData[[counter]] <- setNames(list( outputdata$RRCombined$t , outputdata$RRCombined$RR , outputdata$Meta_Data, InferenceOutput$AFScore) , c('t' , 'RR' , 'MetaData' , 'AFScore'))
  SummaryStatistics[counter , 1] <- abs(difftime(InferenceOutput$AFScore$t[1] , InferenceOutput$AFScore$t[length(InferenceOutput$AFScore$t)] , units = c("hours")))
  SummaryStatistics[counter , 2] <- length(InferenceOutput$AFScore$IHAVFScore)
    
  AFPeriods[[counter]] <- InferenceOutput$StartEndTimesAF
  MMPeriods[[counter]] <-  InferenceOutput$StartEndTimesMM 
  PatientNames[[counter]] <-  outputdata$Meta_Data$PseudoId[1] 
  DP_WaitBar(ii/length(listAllPatients))
  counter <- counter + 1
}
}

RPeaksData <- setNames(RPeaksData , as.vector(unlist(PatientNames)))
AFPeriods <- setNames(AFPeriods , as.vector(unlist(PatientNames)))
#AFPeriods <- AFPeriods[which(unlist(lapply(RPeaksData , function(X){length(X[[1]]) > 10000}) ))]
#RPeaksData <- RPeaksData[which(unlist(lapply(RPeaksData , function(X){length(X[[1]]) > 10000}) ))]

AFStatistics <- matrix(0 , length(RPeaksData) , 3)
AFResults <- matrix(0 , length(RPeaksData) , 1)
AFStatistics <- matrix(0 , length(RPeaksData) , 3)
PatientNames <-  matrix(0 , length(RPeaksData) , 1)
for( ii in 1:length(RPeaksData) )
{
  if( !is.na(RPeaksData[[ii]]$MetaData$FirstNewAF)  && !is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 1] <- 1}
  if( !is.na(RPeaksData[[ii]]$MetaData$FirstNewAF)  &&  is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 2] <- 1}
  if(  is.na(RPeaksData[[ii]]$MetaData$FirstNewAF)  &&  is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 3] <- 1}
  if(  sum(RPeaksData[[ii]]$MetaData$PseudoId[1] ==  ListConfirmedNoAF) > 0 ){AFStatistics[ii , ] <- c(0,0,1)}
  if(  length(RPeaksData[[ii]][[1]]) < 1000){AFStatistics[ii , ] <- c(0,0,1)}
  PatientNames[ii] <- RPeaksData[[ii]]$MetaData$PseudoId
  AFResults[ii , 1] <- length(AFPeriods[[ii]]$Start) 
}


# Filter out records with large time gaps or impossible heartrates
AFResults <- AFResults[((SummaryStatistics[ , 1] < 8)*(SummaryStatistics[ , 1] > 2)*(SummaryStatistics[ , 2]>8000)*(SummaryStatistics[ , 2]<60000)) == 1]
AFStatistics <- AFStatistics[((SummaryStatistics[ , 1] < 8)*(SummaryStatistics[ , 1] > 2)*(SummaryStatistics[ , 2]>8000)*(SummaryStatistics[ , 2]<60000)) == 1 , ]
PatientNames <- PatientNames[((SummaryStatistics[ , 1] < 8)*(SummaryStatistics[ , 1] > 2)*(SummaryStatistics[ , 2]>8000)*(SummaryStatistics[ , 2]<60000)) == 1 ]


Total <- length(AFResults)
N <- apply(AFStatistics , 2 , sum)[3]
P <- apply(AFStatistics , 2 , sum)[1] + apply(AFStatistics , 2 , sum)[2]
PC <- apply(AFStatistics , 2 , sum)[1]

CorrectConfirmed <- (AFStatistics[ , 1] == 1)*( AFResults > 0 )
IncorrectConfirmed <- (AFStatistics[ , 1] == 1)*(AFResults == 0)

Correctunconfirmed <-  (AFStatistics[ , 2] == 1)*( AFResults > 0 )
IncorrectunConfirmed <- (AFStatistics[ , 2] == 1)*( AFResults == 0 )

CorrectNAF <-  (AFStatistics[ , 3] == 1)*(AFResults == 0)
IncorrectNAF <-  (AFStatistics[ , 3] == 1)*(AFResults > 0)

SenistivityConfirmed <- sum(CorrectConfirmed) / PC
SenistivityUnConfirmed <- sum(Correctunconfirmed) / (apply(AFStatistics , 2 , sum)[2])
Sensitivity <- (sum(CorrectConfirmed) + sum(Correctunconfirmed)) / P

PPV <- (sum(CorrectConfirmed) + sum(Correctunconfirmed))/ (sum(IncorrectNAF) + sum(CorrectConfirmed) + sum(Correctunconfirmed))
NPV <- sum(CorrectNAF) /(sum(CorrectNAF) + sum(IncorrectunConfirmed) + sum(IncorrectConfirmed) )

Specifictity <- sum( CorrectNAF ) / N



