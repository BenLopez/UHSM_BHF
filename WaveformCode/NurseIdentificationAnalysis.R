{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )
PatIndex2017 <- DP_LoadPatientIndex()}


TotalUsable <- (PatIndex2017$Usable == 1)*(PatIndex2017$TotalITUTimeHRS <= 100)
TotalUsable[is.na(TotalUsable)] <- FALSE
TotalUsable <- 758

TotalAF <- (PatIndex2017$Usable == 1)*(PatIndex2017$TotalITUTimeHRS <= 100)*(!is.na(PatIndex2017$ConfirmedFirstNewAF))*(PatIndex2017$ConfirmedFirstNewAF != 'CNAF')
TotalAF[is.na(TotalAF)] <- 0
AFPatientRecords <- PatIndex2017[TotalAF == 1,]
TotalAF <- sum(TotalAF)

TotalCNAF <- (PatIndex2017$Usable == 1)*(PatIndex2017$TotalITUTimeHRS <= 100)*(PatIndex2017$ConfirmedFirstNewAF == 'CNAF')
TotalCNAF[is.na(TotalCNAF)] <- 0
TotalCNAF <- sum(TotalCNAF)

TotalNOAF <- (TotalUsable - TotalAF)

TimeDiff <- difftime(DP_StripTime(AFPatientRecords$ConfirmedFirstNewAF) , DP_StripTime(AFPatientRecords$FirstNewAF) , units = 'secs')

NurseEarlyIdentifiction <- sum(TimeDiff > (60*60))
NurseLateIdentifiction <- sum(TimeDiff < -(60*60))

MeanLateDetection <- abs(mean(TimeDiff[TimeDiff < (60*60)]))/(60*60)
RangeLateDetection <- range(TimeDiff[TimeDiff < (60*60)]/(60^2))

plot(TimeDiff, xlab = 'Patient Index' , ylab = 'Time Difference Detection of First New AF')
abline( (60*60) , 0 , col = 'red' )
abline( -(60*60) , 0 , col = 'red')
title('Nurse Detection Analysis')

NurseSpecificity <- (TotalNOAF - TotalCNAF - NurseEarlyIdentifiction) / (TotalNOAF)
NurseSensitivity <- (TotalAF - NurseEarlyIdentifiction - NurseLateIdentifiction )/TotalAF
NursePPV <- (TotalAF - NurseEarlyIdentifiction - NurseLateIdentifiction )/( TotalAF - NurseEarlyIdentifiction - NurseLateIdentifiction + TotalCNAF)
