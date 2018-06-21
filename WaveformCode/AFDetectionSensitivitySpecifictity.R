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

for(ii in 1:length(listAllPatients)){
  
   outputdata <- DP_LoadRpeaksfile(path , listAllPatients[ii])
  
  if(length(outputdata$RRCombined$t) < 10000)
  {
    SummaryStatistics[ii , 2] <- 0
    SummaryStatistics[ii , 1] <- 0
    RPeaksData[[ii]] <- setNames(list(1 , outputdata$Meta_Data) , c('Nonsense' , 'MetaData'))
    next 
  }
  
  InferenceOutput <- AFD_DetectionWrapper( outputdata$RRCombined  )

  RPeaksData[[ii]] <- setNames(list( outputdata$RRCombined$t , outputdata$RRCombined$RR , outputdata$Meta_Data, InferenceOutput$AFScore) , c('t' , 'RR' , 'MetaData' , 'AFScore'))
  SummaryStatistics[ii , 1] <- abs(difftime(InferenceOutput$AFScore$t[1] , InferenceOutput$AFScore$t[length(InferenceOutput$AFScore$t)] , units = c("hours")))
  SummaryStatistics[ii , 2] <- length(InferenceOutput$AFScore$IHAVFScore)
    
  AFPeriods[[ii]] <- InferenceOutput$StartEndTimesAF
  MMPeriods[[ii]] <-  InferenceOutput$StartEndTimesMM 
    
  DP_WaitBar(ii/length(listAllPatients))
}

RPeaksData <- setNames(RPeaksData , as.vector(listAllPatients))

AFStatistics <- matrix(0 , length(RPeaksData) , 3)
AFResults <- matrix(0 , length(RPeaksData) , 1)
AFStatistics <- matrix(0 , length(RPeaksData) , 3)

for( ii in 1:length(RPeaksData) )
{
  if( !is.na(RPeaksData[[ii]]$MetaData$FirstNewAF)  && !is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 1] <- 1}
  if( !is.na(RPeaksData[[ii]]$MetaData$FirstNewAF)  &&  is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 2] <- 1}
  if(  is.na(RPeaksData[[ii]]$MetaData$FirstNewAF)  &&  is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 3] <- 1}
  if(  sum(RPeaksData[[ii]]$MetaData$PseudoId[1] ==  ListConfirmedNoAF) > 0 ){AFStatistics[ii , ] <- c(0,0,1)}
  AFResults[ii , 1] <- length(AFPeriods[[ii]]$Start) 
}

PatientNames <- as.matrix(names(RPeaksData))

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



