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
}

# Filter patients
listAllPatients <- DP_FilterPatients(listAllPatients = listAllPatients , 
                                     PatIndex2017 = PatIndex2017,
                                     HowtoFilterops = HowtoFilterops,
                                     path = path,
                                     FilestoProcess = FilestoProcess)

PatientCDFs <- list()
PatientRR <- list()
PatinetNames <- list()
AFLogical <- matrix(0 , 728 , 1)
counter <- 1
for(ii in 1:length(listAllPatients)){
sub_pat <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = listAllPatients[[ii]])
  if(DP_checkRpeaksfilesprocessed(path , PatientsId = listAllPatients) && DP_CheckFileExists(path , PatientsID = listAllPatients[[ii]] , paste0(listAllPatients[[ii]] , '_CDFs' )) ){
  tmp<- DP_LoadRpeaksfile(path , listAllPatients[[ii]])$RRCombined
  if(length(tmp$t) < 15000){next}
  PatientRR[[counter]] = tmp
  PatientCDFs[[counter]] = DP_LoadFile(path , PatientsID = listAllPatients[[ii]] , paste0(listAllPatients[[ii]] , '_CDFs' )) 
  PatinetNames[[counter]] <- sub_pat$PseudoId[1]
  if(DP_CheckIfAFPatient(sub_pat)){
    AFLogical[counter,1] <- 1
  }
  counter <- counter + 1
  }
DP_WaitBar(ii/length(listAllPatients))
}

PatientRR <- setNames( PatientRR , unlist(PatinetNames) )
PatientCDFs <- setNames( PatientCDFs , unlist(PatinetNames) )

BaselineCDF <- lapply( PatientCDFs , function(X){AFD_CDFCalulateBaseline(X)} )


{patienttoexamine <-'z1027'
med <- as.numeric(apply( PatientCDFs[[patienttoexamine]][,] , 1 , function(X){AFD_CalulateMedianfromBinMatrix(as.numeric(colnames(PatientCDFs[[1]])) , X)}))

index1 <- 10000
index2 <- 20000
x11()
par(mfrow = c(2 , 1))
plot(PatientRR[[patienttoexamine]]$t , PatientRR[[patienttoexamine]]$RR , col = rgb(0,0,1,alpha = 0.01) , pch = 16 , ylab = 'RR' , xlab = 't' )
title('RR Times')
abline(v = as.numeric(PatientRR[[patienttoexamine]]$t[index1]) , col = 'red' )
abline(v = as.numeric(PatientRR[[patienttoexamine]]$t[index1 + 1000]) , col = 'red')
abline(v = as.numeric(PatientRR[[patienttoexamine]]$t[index2]) , col = 'green')
abline(v = as.numeric(PatientRR[[patienttoexamine]]$t[index2 + 1000]) , col = 'green')

plot(as.numeric(colnames(PatientCDFs[[1]])) , PatientCDFs[[patienttoexamine]][index1,] - BaselineCDF[[patienttoexamine]]  , col = rgb(1,0,0 , alpha = 0.01) , type ='l' , xlab = 'x' , ylab = 'P(X - median(X) < x)', ylim = c(-1,1) , xlim =c(0,2))
for( i in 1:1000){
lines(as.numeric(colnames(PatientCDFs[[1]])) , PatientCDFs[[patienttoexamine]][index1 + i,] - BaselineCDF[[patienttoexamine]] , col = rgb(1,0,0 , alpha = 0.01) , type ='l')
lines(as.numeric(colnames(PatientCDFs[[1]])) , PatientCDFs[[patienttoexamine]][index2 + i,] - BaselineCDF[[patienttoexamine]], col = rgb(0,1,0 , alpha = 0.01) , type ='l')
}
}


AFBaselineCDFMatrices <- matrix(0 , sum(AFLogical ==1) , length( PatientCDFs[[1]][1,] ) )
NAFBaselineCDFMatrices <- matrix(0 , sum(AFLogical ==0) , length( PatientCDFs[[1]][1,] ) )
counter1 <- 1 
counter2 <- 1
for( ii in 1:length(PatientCDFs) ){
  if(AFLogical[ii] == 1){
  AFBaselineCDFMatrices[counter1,] <- BaselineCDF[[ii]]
  counter1 <- counter1 + 1
  }else{
  NAFBaselineCDFMatrices[counter2,] <- BaselineCDF[[ii]]
  counter2 <- counter2 + 1
  }
}

m_BL_AF <- apply(AFBaselineCDFMatrices[ ,5:27 ] , 2 , function(X){mean(X[!is.na(X)])})
v_BL_AF <- cov(AFBaselineCDFMatrices[ ,5:27 ])

m_BL_NAF <- apply(NAFBaselineCDFMatrices[ ,5:27 ] , 2 , function(X){mean(X[!is.na(X)])})
v_BL_NAF <- cov(NAFBaselineCDFMatrices[ ,5:27 ])


diffindit <- (m_BL_AF - m_BL_NAF)%*%solve(v_BL_NAF + v_BL_AF )%*%(m_BL_AF - m_BL_NAF)


plot( m_BL_AF , type ='l' , col = 'red' )
lines(m_BL_AF + 2*sqrt( diag(v_BL_NAF) ) , type ='l' , col = 'blue')
lines(m_BL_AF - 2*sqrt( diag(v_BL_NAF) ) , type ='l'  , col = 'blue')

lines( m_BL_NAF , type ='l' , col = 'black' )
lines(m_BL_NAF + 2*sqrt( diag(v_BL_NAF) ) , type ='l' , col = 'green')
lines(m_BL_NAF - 2*sqrt( diag(v_BL_NAF) ) , type ='l'  , col = 'green')

binlims=  c( seq(from = 0  , to = 1  , 0.1  ))
par(mfrow = c(8,2))
for( ii in 5:27 ){
  x11(20,10)
  hist(AFBaselineCDFMatrices[ , ii] , col = rgb(1 , 0 , 0 , alpha = 0.25) , breaks = binlims , axes = FALSE , freq=FALSE)
  hist(NAFBaselineCDFMatrices[ , ii] , col = rgb(0 , 0 , 1 , alpha = 0.25) , breaks = binlims, add = TRUE , axes = FALSE , freq=FALSE)
}


AFBaselineCDFMatrices[is.na(AFBaselineCDFMatrices)] <- 0
pairs(AFBaselineCDFMatrices[ ,7:20 ],col = rgb(1 , 0 , 0 , alpha = 0.1) , pch = 16)

pairs(NAFBaselineCDFMatrices[ ,7:20 ],col = rgb(0 , 1 , 0 , alpha = 0.1) , pch = 16)

x11()
PlotFun_DoublePairs(AFBaselineCDFMatrices[ ,7:20 ] , NAFBaselineCDFMatrices[ ,7:20 ])

