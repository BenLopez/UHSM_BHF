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
  set.seed(1)
}

{
  PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z209'] = "07/10/2016 10:46"
  PatIndex2017$EndFirstNewAF[PatIndex2017$PseudoId == 'z209'] = "07/10/2016 11:29"
  
  PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z925'] = "10/06/2017 08:31"
  PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z401'] = "04/12/2016 19:27"
  PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z1281'] = "26/10/2017 01:04"
  PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z580'] = "04/02/2017 06:10"
}

TotalUsable <- length(listAllPatients)

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

plot1 <- ggplot(data = data.frame(y = 1:length(as.numeric(TimeDiff[-which(abs(TimeDiff)<(60^2))[1:6]]/(60*60))) , x = as.numeric(TimeDiff[-which(abs(TimeDiff)<(60^2))[1:6]]/(60*60))) , aes(y,x))+
  geom_point( color = 'blue')+
  ylab('Difference between Diganosis of FNAF (hours)') +
  xlab('AF Patient Index') +
  ggtitle('Comparisons of Diagnosis Times for FNAF') +
  geom_hline(yintercept = 1 , color = 'red')+
  geom_hline(yintercept = -1 , color = 'red')

plot(TimeDiff[-which(abs(TimeDiff)<(60^2))[1:6]]/(60*60), xlab = 'Patient Index' , ylab = 'Time Difference Detection of First New AF')
abline( (1) , 0 , col = 'red' )
abline( -(1) , 0 , col = 'red')
title('Nurse Detection Analysis')

NurseSpecificity <- (TotalNOAF - TotalCNAF - NurseEarlyIdentifiction) / (TotalNOAF)
NurseSensitivity <- (TotalAF - NurseEarlyIdentifiction - NurseLateIdentifiction )/TotalAF
NursePPV <- (TotalAF - NurseEarlyIdentifiction - NurseLateIdentifiction )/( TotalAF - NurseEarlyIdentifiction - NurseLateIdentifiction + TotalCNAF)
