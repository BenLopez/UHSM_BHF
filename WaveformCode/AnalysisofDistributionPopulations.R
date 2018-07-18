{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

DP_LoadPatientIndex()
HoursBeforeAndAFter <- DP_SelectHoursBeforeandAfter()
FilestoProcess <- DP_ChooseECGstoProcess() 
HowtoFilterops <- read.csv(choose.files(caption = "Select listofopsSH") , stringsAsFactors = FALSE)

DP_ChooseDataReps()
listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
listAllPatients <- select.list(listAllPatients , graphics = TRUE , multiple = TRUE)

# Load and extract data.
{
counter <- 1
DataMatrix <- array(0 , dim = c(length(listAllPatients) , 60000 , 10) )
PatientNames <- list()
PatientDetails <- list()
PatientAF <- matrix(0 , length(listAllPatients) , 1)

for(ii in 1:length(listAllPatients)){
if(DP_CheckDistributionSummariesExists(path , listAllPatients[[ii]]) == FALSE ){next}
sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[[ii]])
if(sub_pat$PseudoId == 'z1007' || sub_pat$PseudoId == 'z1000'){next}
DistributionSummaries <- DP_LoadDistributionSummaries(path , listAllPatients[[ii]]) 
for(jj in 1:length(DistributionSummaries)){
DataMatrix[counter , 1:length(DistributionSummaries[[jj]]) , jj] <- DistributionSummaries[[jj]]
}
PatientNames[[counter]] <- listAllPatients[[ii]]
PatientDetails[[counter]] <- sub_pat
if(!is.na(PatientDetails[[counter]]$ConfirmedFirstNewAF) &  PatientDetails[[counter]]$ConfirmedFirstNewAF != 'CNAF'){
  PatientAF[counter] <- 1 
}
counter <- counter +1
}

DataMatrix[is.na(DataMatrix)] <- 0
PatientAF <- PatientAF[apply(DataMatrix  , 1 , sum) != 0 , ]
DataMatrix <- DataMatrix[ , apply(DataMatrix , 2 , sum) != 0 , ]
DataMatrix <- DataMatrix[ apply(DataMatrix  , 1 , sum) != 0 , , ]

AFMatrix <- DataMatrix[PatientAF == 1 , 1000:dim( DataMatrix )[2] , ]
NAFMatrix <- DataMatrix[PatientAF == 0 , 1000:dim( DataMatrix )[2]  , ]

AFPatients <- PatientNames[PatientAF == 1]
NAFPatients <- PatientNames[PatientAF == 0]

AFPatientDetails  <- PatientDetails[PatientAF == 1]
NAFPatientDetails <- PatientDetails[PatientAF == 0]

rm(DataMatrix)
}

dev.off()
variabletoview <- 10
plot(cumsum(NAFMatrix[1 ,NAFMatrix[1 , ,1] != 0 ,1]) , NAFMatrix[1 ,NAFMatrix[1 , ,1] != 0 ,variabletoview] , type = 'l' , 
     col = rgb(1 , 0 , 0 , alpha = 0.1) ,
     xlab ='t' , 
     ylab = names(DistributionSummaries)[variabletoview]   , 
     ylim = 2*range(NAFMatrix[1 , ,variabletoview])  , 
     xlim = c(0,25000))
title('Heart Rate Distribution Summary Analysis')
for( ii in 2:52){
  lines(cumsum(NAFMatrix[ii ,NAFMatrix[ii , ,1] != 0 ,1]) , NAFMatrix[ii , NAFMatrix[ii , ,1] != 0 ,variabletoview] , col = rgb(1 , 0 , 0 , alpha = 0.1))
}

for( ii in 2:size(AFMatrix)[1]){
  lines(cumsum(AFMatrix[ii ,AFMatrix[ii , ,1] != 0 ,1]) , AFMatrix[ii , AFMatrix[ii , ,1] != 0 ,variabletoview] , col = rgb(0 , 1 , 0 , alpha = 0.1))
}
abline(v = 15000)
abline(v = 17000 , col = 'blue')  

AFVector <- matrix(0 , 1 , 10)
PreAFVector <- matrix(0 , 1 , 10)
for(ii in 1:size(AFMatrix)[1]){
  tmp <- AFMatrix[ii ,AFMatrix[ii , ,1] != 0 ,]
  ctmp <- cumsum(tmp[  , 1])
  if(sum(ctmp > 17000) > 1){
    AFVector <- rbind(AFVector , tmp[ctmp > 17000 , ]) }
  if(sum(ctmp > 15000) > 1){
    PreAFVector <- rbind(PreAFVector , tmp[ctmp < 15000 , ])
  }
}

NAFVector<- matrix(0 , 1 , 10)
for(ii in 1:100){
  tmp <- NAFMatrix[ii ,NAFMatrix[ii , ,1] != 0 ,]
  ctmp <- cumsum(tmp[  , 1])
  if(sum(ctmp > 15000) > 1){
    NAFVector <- rbind(NAFVector , tmp[ctmp < 15000 , ] )
  }
}

AFVector <- matrix(0 , 1 , 10)
PreAFVector <- matrix(0 , 1 , 10)
for(ii in 1:size(AFMatrix)[1]){
  tmp <- AFMatrix[ii ,AFMatrix[ii , ,1] != 0 ,]
  ctmp <- cumsum(tmp[ , 1])
  if(nrow(tmp)< 15000){next}
  for(jj in 1:size(tmp)[2]){
    tmp[ , jj] <- tmp[ , jj] - mean(tmp[1:5000 , jj] , na.rm = TRUE)  
  }
  if(sum(ctmp > 17000) > 1){
  AFVector <- rbind(AFVector , tmp[ctmp > 17000 , ] ) }
  if(sum(ctmp > 15000) > 1){
  PreAFVector <- rbind(PreAFVector , tmp[ctmp < 15000 , ] )
  }
}

NAFVector<- matrix(0 , 1 , 10)
for(ii in 1:size(NAFMatrix)[1]){
  tmp <- NAFMatrix[ii ,NAFMatrix[ii , ,1] != 0 ,]
  ctmp <- cumsum(tmp[  , 1])
  if(nrow(tmp) < 15000){next}
  for(jj in 1:size(tmp)[2]){
    tmp[ , jj] <- tmp[ , jj] - mean(tmp[1:5000 , jj] , na.rm = TRUE)  
  }
  if(sum(ctmp > 15000) > 1){
    NAFVector <- rbind(NAFVector , tmp[ctmp < 15000 , ] )
  }
}


variabletoview <- 1

x11()
par(mfrow = c(2 , 5))

for(variabletoview in c(1:10)){
hist(AFVector[ , variabletoview], col=rgb(0,0,1,alpha = 0.1) ,
     main= paste0(names(DistributionSummaries)[variabletoview] , ' Histogram') , xlab = names(DistributionSummaries)[variabletoview] , freq = FALSE)
hist(PreAFVector[ , variabletoview], col=rgb(0,1,0,alpha =0.5), add=T , freq = FALSE)
hist(NAFVector[ , variabletoview], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE)
}


mAF <- apply(AFVector , 2 , mean)
SigmaAF <- cov(AFVector)

mPreAF <- apply(PreAFVector , 2 , mean)
SigmaPreAF <- cov(PreAFVector)

mNAF <- apply(NAFVector , 2 , mean)
SigmaNAF <- cov(NAFVector)


patientindex <- 2
tmp <- apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 , ] , 1 ,  function(X){t(X - mAF)%*%solve(SigmaAF)%*%(X - mAF)})
tmp2 <- apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 , ] , 1 ,  function(X){t(X - mPreAF)%*%solve(SigmaPreAF)%*%(X - mPreAF)})
tmp3 <- apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 , ] , 1 ,  function(X){t(X - mNAF)%*%solve(SigmaNAF)%*%(X - mNAF)})

stmp = 1/tmp + 1/tmp2 + 1/tmp3
plot((1/tmp)/stmp , type='l' , col = 'blue')
lines((1/tmp2)/stmp , type='l' , col = 'green')
lines((1/tmp3)/stmp , type='l' , col = 'red')

patientindex <- 2
tmp <- apply(NAFMatrix[patientindex ,  , ] , 1 ,  function(X){t(X - mAF)%*%solve(SigmaAF)%*%(X - mAF)})
tmp2 <- apply(NAFMatrix[patientindex ,  , ] , 1 ,  function(X){t(X - mPreAF)%*%solve(SigmaPreAF)%*%(X - mPreAF)})
tmp3 <- apply(NAFMatrix[patientindex , , ] , 1 ,  function(X){t(X - mNAF)%*%solve(SigmaNAF)%*%(X - mNAF)})

stmp = (1/tmp) + (1/tmp2) + (1/tmp3)
plot((1/tmp)/stmp , type='l' , col = 'blue' , xlab = 't' , ylab = 'Probability' , ylim = c(0,1) , xlim=c(0,30000))
lines((1/tmp2)/stmp , type='l' , col = 'green')
lines((1/tmp3)/stmp , type='l' , col = 'red')
