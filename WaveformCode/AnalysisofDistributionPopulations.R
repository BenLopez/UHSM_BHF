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
listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)


subList <- DP_choosepatients(listAllPatients )
# Load and extract data.
{
counter <- 1
DataMatrix <- array(0 , dim = c(length(listAllPatients) , 50000 , 11) )
PatientNames <- list()
PatientDetails <- list()
PatientAF <- matrix(0 , length(listAllPatients) , 1)

for(ii in 1:length(listAllPatients)){
if(DP_CheckDistributionSummariesExists(path , listAllPatients[[ii]]) == FALSE ){
  DP_WaitBar(ii/length(listAllPatients))  
  next}

sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[[ii]])

if(sub_pat$PseudoId == 'z1007' || sub_pat$PseudoId == 'z1000'){
  next}

DistributionSummaries <- DP_LoadDistributionSummaries(path , listAllPatients[[ii]]) 
if(length(DistributionSummaries[[1]]) > 50000){next}
if(length(DistributionSummaries[[1]]) < 5000){next}

for(jj in 1:length(DistributionSummaries)){
DataMatrix[counter , 1:min(length(DistributionSummaries[[jj]]) , 50000) , jj] <- DistributionSummaries[[jj]][1:min(length(DistributionSummaries[[jj]]) , 50000)]
#DataMatrix[counter , 1:length(DistributionSummaries[[jj]])  , jj] <- DistributionSummaries[[jj]][1:(length(DistributionSummaries[[jj]]) )]

}
PatientNames[[counter]] <- listAllPatients[[ii]]
PatientDetails[[counter]] <- sub_pat
if(!is.na(PatientDetails[[counter]]$ConfirmedFirstNewAF) &  PatientDetails[[counter]]$ConfirmedFirstNewAF != 'CNAF'){
  PatientAF[counter] <- 1 
}
counter <- counter +1
DP_WaitBar(ii/length(listAllPatients))
}

DataMatrix[is.na(DataMatrix)] <- 0
PatientAF <-  PatientAF[apply(DataMatrix  , 1 , sum) != 0 , ]

NAFMatrix <- DataMatrix[PatientAF == 0 , 1000:dim( DataMatrix )[2]  , ]
DataMatrix <- DataMatrix[PatientAF == 1 , 1000:dim( DataMatrix )[2]  , ]
NAFMatrix <- NAFMatrix[1:250,,]
AFMatrix <- DataMatrix

#AFMatrix <- AFMatrix[ , apply(NAFMatrix , 2 , sum) != 0 , ]
#AFMatrix <- AFMatrix[ apply(NAFMatrix  , 1 , sum) != 0 , , ]
#NAFMatrix <- NAFMatrix[ , apply(NAFMatrix , 2 , sum) != 0 , ]
#NAFMatrix <- NAFMatrix[ apply(NAFMatrix  , 1 , sum) != 0 , , ]

AFPatients <- PatientNames[PatientAF == 1]
NAFPatients <- PatientNames[PatientAF == 0]
NAFPatients <- NAFPatients[1:250]

AFPatientDetails  <- PatientDetails[PatientAF == 1]
NAFPatientDetails <- PatientDetails[PatientAF == 0]

rm(DataMatrix)
}

dev.off()
variabletoview <- 1
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

# Extract data with normalising constant.
{AFVector <- matrix(0 , 1 , 10)
PreAFVector <- matrix(0 , 1 , 10)
for(ii in 1:size(AFMatrix)[1]){
  tmp <- AFMatrix[ii ,AFMatrix[ii , ,1] != 0 ,]
  AFVector <- rbind(AFVector , tmp[tmp[ , 11] > 0 ,1:10]  )
  PreAFVector <- rbind(PreAFVector , tmp[tmp[ , 11] < 0 ,1:10])
  DP_WaitBar(ii/size(AFMatrix)[1])
}

NAFVector<- matrix(0 , 1 , 10)
for(ii in 1:size(NAFMatrix)[1]){
    tmp <- NAFMatrix[ii ,NAFMatrix[ii , ,1] != 0 ,]
    NAFVector <- rbind(NAFVector , tmp[,1:10] )
    DP_WaitBar(ii/size(NAFMatrix)[1])
}
}

# 
AFVector <- matrix(0 , 1 , 10)
PreAFVector <- matrix(0 , 1 , 10)
mmAF <- matrix(0 ,  size(AFMatrix)[1],10)
mmNAF <- matrix(0 ,  size(NAFMatrix)[1],10)
for(ii in 1:size(AFMatrix)[1]){
  tmp <- AFMatrix[ii ,AFMatrix[ii , ,1] != 0 ,]
  if(nrow(tmp)< 5000){
    mmAF[ii,] <- apply( tmp[,1:10],2,function(X){mean(X , na.rm = TRUE)} )
    next}
    mmAF[ii,] <- apply(tmp[1:5000,1:10],2,function(X){mean(X , na.rm = TRUE)})
  for(jj in 1:(size(tmp)[2] -1)){
    tmp[ , jj] <- tmp[ , jj] - mmAF[ii,jj] 
  }
    AFVector <- rbind(AFVector , tmp[tmp[ , 11] > 0 ,1:10]  )
    PreAFVector <- rbind(PreAFVector , tmp[tmp[ , 11] < 0 ,1:10])
    DP_WaitBar(ii/size(AFMatrix)[1])
}

NAFVector<- matrix(0 , 1 , 10)
patientindexsample <- sample(1:size(NAFMatrix)[1] , size(AFMatrix)[1] , replace = FALSE)
for(ii in 1:length(patientindexsample)){
  tmp <- NAFMatrix[patientindexsample[ii] ,NAFMatrix[patientindexsample[ii] , ,1] != 0 ,]
if(nrow(tmp) < 5000){
  mmNAF[ii,] <- apply(tmp[,1:10],2,function(X){mean(X , na.rm = TRUE)})
  next}
  mmNAF[ii,] <- apply(tmp[1:5000,1:10],2,function(X){mean(X , na.rm = TRUE)})
    for(jj in 1:(size(tmp)[2] - 1)){
    tmp[ , jj] <- tmp[ , jj] - mmNAF[ii,jj]
  }
  NAFVector <- rbind(NAFVector , tmp[,1:10] )
  DP_WaitBar(ii/size(NAFMatrix)[1])
}

variabletoview <- 1

{x11(20,14)
par(mfrow = c(2 , 5))

for(variabletoview in c(1:10)){
tmp <- hist(AFVector[ , variabletoview], col=rgb(0,0,1,alpha = 0.5) ,
     main= paste0(names(DistributionSummaries)[variabletoview] , ' Histogram') , xlab = names(DistributionSummaries)[variabletoview] , freq = FALSE)
hist(PreAFVector[ , variabletoview], col=rgb(0,1,0,alpha =0.5), add=T , freq = FALSE , breaks = c(min(PreAFVector[ , variabletoview]) , tmp$breaks, max(PreAFVector[ , variabletoview])))
hist(NAFVector[ , variabletoview], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(min(NAFVector[ , variabletoview]) , tmp$breaks, max(NAFVector[ , variabletoview]) ))
}
}

{x11(20,14)
  par(mfrow = c(2 , 5))
  
  for(variabletoview in c(1:10)){
    tmp <- hist(NAFVector[ , variabletoview], col=rgb(1,0,0,alpha = 0.5) ,
                main= paste0(names(DistributionSummaries)[variabletoview] , ' Histogram') , xlab = names(DistributionSummaries)[variabletoview] , freq = FALSE)
    hist(PreAFVector[ , variabletoview], col=rgb(0,1,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(min(PreAFVector[ , variabletoview]) , tmp$breaks, max(PreAFVector[ , variabletoview]) ))
  }
}

{x11()
  par(mfrow = c(2 , 5))
  for(variabletoview in c(1:10)){
    tmp <- hist(mmAF[ , variabletoview], col=rgb(0,0,1,alpha = 0.1) ,
                main= paste0(names(DistributionSummaries)[variabletoview] , ' Histogram') , xlab = names(DistributionSummaries)[variabletoview] , freq = FALSE)
    hist(mmNAF[ , variabletoview], col=rgb(0,1,0,alpha =0.5), add=T , freq = FALSE , breaks = c(min(mmNAF[!is.na(mmNAF[ , variabletoview]) , variabletoview] ) , tmp$breaks, max(mmNAF[!is.na(mmNAF[ , variabletoview]) , variabletoview])))
  }
}


Im1 <- mahalanobis(AFVector, mAF, SigmaAF)
Im2 <- mahalanobis(NAFVector, mNAF, SigmaNAF)
Im3 <- mahalanobis(PreAFVector, mPreAF, SigmaPreAF)


mAF <- apply(AFVector[,  ] , 2 , mean)
SigmaAF <- cov(AFVector[,  ])

mPreAF <- apply(PreAFVector[,  ] , 2 , mean)
SigmaPreAF <- cov(PreAFVector[,  ])

mNAF <- apply(NAFVector[,  ] , 2 , mean)
SigmaNAF <- cov(NAFVector[,  ])

x11()
SAMPLENUM <- 10000
variable = c(1:10)
variable = variable[ ]
tmp <- rbind(AFVector[sample(1:size(AFVector)[1] , SAMPLENUM ),variable],NAFVector[sample(1:size(NAFVector)[1] , SAMPLENUM ),variable] , PreAFVector[sample(1:size(PreAFVector)[1] , 1000 ),variable])
al = 0.01
colvector <- c(rep(rgb(1 , 0 , 0, alpha = al  ) , SAMPLENUM) , rep(rgb(0 ,0 , 1, alpha = al)   , SAMPLENUM) ,  rep(rgb(0 , 1 , 0, alpha = al)   , SAMPLENUM))
pairs(tmp , pch = 16 , col = colvector , labels = names(DistributionSummaries)[variable])

x11(30,20)
SAMPLENUM <- 10000
variable = c(1:10)
variable = variable[ ]
tmp <- rbind(NAFVector[sample(1:size(NAFVector)[1] , SAMPLENUM ),variable] , PreAFVector[sample(1:size(PreAFVector)[1] , 1000 ),variable])
al = 0.01
colvector <- c(rep(rgb(0 ,0 , 1, alpha = al)   , SAMPLENUM) ,  rep(rgb(0 , 1 , 0, alpha = al)   , SAMPLENUM))
pairs(tmp , pch = 16 , col = colvector , labels = names(DistributionSummaries)[variable])


{
patientindex <- 2
n <- 10
x11()
outputData <- DP_LoadRpeaksfile(path , AFPatients[patientindex])
tmp <- apply( t(apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmAF[patientindex,  ]} ))  , 1 ,  function(X){t(X - mAF)%*%solve(SigmaAF)%*%(X - mAF)})
tmp2 <- apply(t(apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmAF[patientindex,  ]} )) , 1 ,  function(X){t(X - mPreAF)%*%solve(SigmaPreAF)%*%(X - mPreAF)})
tmp3 <- apply(t(apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmAF[patientindex,  ]} )), 1 ,  function(X){t(X - mNAF)%*%solve(SigmaNAF)%*%(X - mNAF)})

stmp = 1/tmp + 1/tmp2 + 1/tmp3

culmulativeprobability <- function(X , Y , Z , p1= 1/3 , p2= 1/3 , p3 = 1/3 , n = 100){
  P1 <- exp(n*smth( log(X) , method = 'sma' , n=n ))
  P2 <- exp(n*smth( log(Y) , method = 'sma' , n=n ))
  P3 <- exp(n*smth( log(Z) , method = 'sma' , n=n ))
  return( (p1*P1)/ ((p1*P1)+(p2*P2)+(p3*P3)) )
  }

# P1 <- pchisq( tmp  ,  10 , lower.tail = FALSE )
# P2 <- pchisq( tmp2 ,  10 , lower.tail = FALSE )
# P3 <- pchisq( tmp3 ,  10 , lower.tail = FALSE )

x <- t(apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmAF[patientindex,  ]} ))

P1 <- exp(mvnpdf( x , mAF ,  SigmaAF ))
P3 <- exp(mvnpdf( x , mNAF  ,Sigma = SigmaNAF   ))
P2 <- exp(mvnpdf( x , mPreAF , Sigma = SigmaPreAF  ))

probabilities <- c(0.00003 , 0.00017 , (1 -0.00003 - 0.00017))
plot(culmulativeprobability(X=P1 , Y=P2, Z=P3 , p1 = probabilities[1] , p2 = probabilities[2] , p3 = probabilities[3] , n = n) , type ='l' , col = 'blue' , ylim = c(0,1))
lines(culmulativeprobability(X=P2 , Y=P1, Z=P3, p1 = probabilities[2] , p2 = probabilities[1] , p3 = probabilities[3] , n = n) , type ='l' , col = 'green')
lines(culmulativeprobability(X=P3 , Y=P1, Z=P2, p1 = probabilities[3] , p2 = probabilities[1] , p3 = probabilities[2] , n = n) , type ='l' , col = 'red')
title(AFPatients[patientindex])
}



{
  x11()
  n <- 11
  patientindex <- 8
  outputData <- DP_LoadRpeaksfile(path , NAFPatients[patientindex])
  tmp <- apply( t(apply(NAFMatrix[patientindex , NAFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmNAF[patientindex,  ]} ))  , 1 ,  function(X){t(X - mAF)%*%solve(SigmaAF)%*%(X - mAF)})
  tmp2 <- apply(t(apply(NAFMatrix[patientindex , NAFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmNAF[patientindex,  ]} )) , 1 ,  function(X){t(X - mPreAF)%*%solve(SigmaPreAF)%*%(X - mPreAF)})
  tmp3 <- apply(t(apply(NAFMatrix[patientindex , NAFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmNAF[patientindex,  ]} )), 1 ,  function(X){t(X - mNAF)%*%solve(SigmaNAF)%*%(X - mNAF)})
  stmp = 1/tmp + 1/tmp2 + 1/tmp3
  
  #P1 <- pchisq( tmp  ,  10 , lower.tail = FALSE )
  #P2 <- pchisq( tmp2 ,  10 , lower.tail = FALSE )
  #P3 <- pchisq( tmp3 ,  10 , lower.tail = FALSE )
  
  x <- t(apply(NAFMatrix[patientindex , NAFMatrix[patientindex , , 1] !=0 , 1:10 ] , 1 , function(X){X - mmNAF[patientindex,  ]} ))
  
  P1 <- exp(mvnpdf( x , mAF ,  SigmaAF ))
  P3 <- exp(mvnpdf( x , mNAF  ,Sigma = SigmaNAF   ))
  P2 <- exp(mvnpdf( x , mPreAF , Sigma = SigmaPreAF  ))
  
  
  probabilities <- c(0.00003 , 0.00017 , (1 -0.00003 - 0.00017))
  plot(culmulativeprobability(X=P1 , Y=P2, Z=P3 , p1 = probabilities[1] , p2 = probabilities[2] , p3 = probabilities[3] , n = n) , type ='l' , col = 'blue' , ylim = c(0,1))
  lines(culmulativeprobability(X=P2 , Y=P1, Z=P3, p1 = probabilities[2] , p2 = probabilities[1] , p3 = probabilities[3] , n = n) , type ='l' , col = 'green')
  lines(culmulativeprobability(X=P3 , Y=P1, Z=P2, p1 = probabilities[3] , p2 = probabilities[1] , p3 = probabilities[2] , n = n) , type ='l' , col = 'red')
  title(NAFPatients[patientindex])
}


f1 <- mvnpdf( x , mAF ,  SigmaAF )
f2 <- mvnpdf( x , mNAF  ,Sigma = SigmaNAF   )
f3 <- mvnpdf( x , mPreAF , Sigma = SigmaPreAF  )


mAF <- apply(AFVector[,  ] , 2 , mean)
SigmaAF <- cov(AFVector[,  ])

mNotAF <- apply(rbind(PreAFVector[,  ] , NAFVector[,  ]) , 2 , mean)
SigmaNotAF <- cov(rbind(PreAFVector[,  ] , NAFVector[,  ]))

{
  x11()
  patientindex <- 11
  tmp <- apply( t(apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 ,  ] , 1 , function(X){X - mmAF[patientindex,  ]} ))  , 1 ,  function(X){t(X - mAF)%*%solve(SigmaAF)%*%(X - mAF)})
  tmp2 <- apply(t(apply(AFMatrix[patientindex , AFMatrix[patientindex , , 1] !=0 ,  ] , 1 , function(X){X - mmAF[patientindex,  ]} )) , 1 ,  function(X){t(X - mNotAF)%*%solve(SigmaNotAF)%*%(X - mNotAF)})
  
  stmp = 1/tmp + 1/tmp2  
    plot((1/tmp)/stmp , type='l' , col = 'blue' , ylim = c(0,1) , xlab=c('time') , ylab=c('Probability'))
  lines((1/tmp2)/stmp , type='l' , col = 'green')
  title(AFPatients[patientindex])
}

{
  x11()
  patientindex <- 1
  tmp <- apply( t(apply(NAFMatrix[patientindex , NAFMatrix[patientindex , , 1] !=0 ,  ] , 1 , function(X){X - mmNAF[patientindex,  ]} ))  , 1 ,  function(X){t(X - mAF)%*%solve(SigmaAF)%*%(X - mAF)})
  tmp2 <- apply(t(apply(NAFMatrix[patientindex , NAFMatrix[patientindex , , 1] !=0 ,  ] , 1 , function(X){X - mmNAF[patientindex,  ]} )) , 1 ,  function(X){t(X - mNotAF)%*%solve(SigmaNotAF)%*%(X - mNotAF)})
  
  stmp = 1/tmp + 1/tmp2 + 
  plot((1/tmp)/stmp , type='l' , col = 'blue' , ylim = c(0,1) , xlab=c('time') , ylab=c('Probability'))
  lines((1/tmp2)/stmp , type='l' , col = 'green')
  title(NAFPatients[patientindex])
}



#### GMM denisty ####

densAF <- densityMclust( AFVector[ , 1:10] , G = 5)
densNAF <- densityMclust( NAFVector[ , 1:10] , G = 5)
densPreAF <- densityMclust( PreAFVector[ , 1:10] , G = 5)

summary(dens)
tmp1 <- predict(densAF , PreAFVector[10000:11000 , 1:10] , what = c('dens'))
tmp2 <- predict(densNAF , PreAFVector[10000:11000 , 1:10] , what = c('dens'))
tmp3 <- predict(densPreAF , PreAFVector[10000:11000 , 1:10] , what = c('dens'))


plot(tmp1/(tmp1 + tmp2 + tmp3))
