pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )


DP_LoadPatientIndex()
HowtoFilterops <- read.csv(choose.files(caption = "Select listofopsSH") , stringsAsFactors = FALSE)

DP_ChooseDataReps()
listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
Patientsprocessedlogical <- lapply(listAllPatients , function(X){DP_checkRpeaksfilesprocessed(path , X ) } )
listAllPatients <- listAllPatients[which(Patientsprocessedlogical == TRUE)]

ListConfirmedNoAF <- c( 'z380' , 'z950' , 'z1086' , 'z1163' , 'z1203' )

RPeaksData <- list()
SummaryStatistics <- matrix(0 , length(listAllPatients) , 2)
  
for(ii in 1:length(listAllPatients))
{
  
  load(paste0(path , '\\' , listAllPatients[ii] , '\\Zip_out\\' , listAllPatients[ii] , '_RPeaks' , '.RData'))
  
  AFScore <- ExtractIHVAFScore(outputdata$ECGI ,  binlims <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3))
  RPeaksData[[ii]] <- setNames(list( outputdata$ECGI$t , outputdata$ECGI$RR , outputdata$Meta_Data, AFScore) , c('t' , 'RR' , 'MetaData' , 'AFScore'))
  SummaryStatistics[ii , 1] <- abs(difftime(outputdata$ECGI$t[1] , outputdata$ECGI$t[length(outputdata$ECGI$t)] , units = c("hours")))
  SummaryStatistics[ii , 2] <- length(outputdata$ECGI$t)
    
  print(ii/length(listAllPatients))
}

RPeaksData <- setNames(RPeaksData , as.vector(listAllPatients))

AFPeriods <- list()
AFStatistics <- matrix(0 , length(RPeaksData) , 3)
for( ii in 1:length(RPeaksData) )
{
  
  AFPeriods[[ii]] <- ASWF_GetStartEndAF(RPeaksData[[ii]]$AFScore$t , logicaltimeseries = (RPeaksData[[ii]]$AFScore$IHAVFScore > 145)  , minutethreshold = 9)
  print(ii/length(RPeaksData))
}

AFResults <- matrix(0 , length(RPeaksData) , 1)
AFStatistics <- matrix(0 , length(RPeaksData) , 3)
for( ii in 1:length(RPeaksData) )
{
  if( !is.na(RPeaksData[[ii]]$MetaData$FirstNewAF) && !is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 1] <- 1}
  if( !is.na(RPeaksData[[ii]]$MetaData$FirstNewAF) &&  is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 2] <- 1}
  if(  is.na(RPeaksData[[ii]]$MetaData$FirstNewAF)  && is.na(RPeaksData[[ii]]$MetaData$ConfirmedFirstNewAF)){AFStatistics[ii , 3] <- 1}
  if(  sum(RPeaksData[[ii]]$MetaData$PseudoId ==  ListConfirmedNoAF) > 0 ){AFStatistics[ii , ] <- c(0,0,1)}
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
NPV <- sum(CorrectNAF) /(sum(CorrectNAF) + sum(IncorrectunConfirmed) + sum(IncorrectunConfirmed) )

Specifictity <- sum( CorrectNAF ) / N



