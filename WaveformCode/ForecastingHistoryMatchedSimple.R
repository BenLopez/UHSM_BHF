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

{For_FillEmptyRecords<-function(PatientData){
  for( i in 2:dim(PatientData)[1] ){
    PatientData[i,is.na(PatientData[i,])] <- PatientData[i-1,is.na(PatientData[i,])]
    
  }
  return(PatientData)
}
  For_BuildSOSFromData <- function(DataStructureForBayesLinear , m =10 , n=20){
    numberofvariables <- dim(DataStructureForBayesLinear)[2]
    numberofpatients <- dim(DataStructureForBayesLinear)[1]
    tmpmatrix <- matrix( DataStructureForBayesLinear[ , , (n+1):(n+m)], numberofpatients , numberofvariables*m )
    tmpmatrix <- tmpmatrix[DP_FindNARows(tmpmatrix) , ]
    
    E_D <- apply(tmpmatrix , 2, function(X){mean(X,na.rm=T)})
    V_D <- cov(tmpmatrix)
    V_D <- DP_AddNugget(V_D , 1e-8*diag(diag(V_D)) )
    
    tmpmatrix <- matrix( DataStructureForBayesLinear[ , , 1:(n)], numberofpatients , numberofvariables*n )
    tmpmatrix <- tmpmatrix[DP_FindNARows(tmpmatrix) , ]
    
    E_F <- apply(tmpmatrix , 2, function(X){mean(X,na.rm=T)})
    V_F <- cov(tmpmatrix)
    
    tmpmatrix <- matrix( DataStructureForBayesLinear[ , , 1:(n+m)], numberofpatients , numberofvariables*(n+m) )
    tmpmatrix <- tmpmatrix[DP_FindNARows(tmpmatrix) , ]
    
    cov_FD <- cov( tmpmatrix )[1:(n*numberofvariables),(n*numberofvariables + 1):((m +n)*numberofvariables)]
    
    weightmatrix <- cov_FD%*%solve(V_D)
    
    output <- setNames(list(E_D , V_D, E_F , V_F , cov_FD,weightmatrix) , c( 'E_D' , 'V_D', 'E_F' , 'V_F' , 'cov_FD' , 'weightmatrix' ) )
    return(output)
  }
  For_BayesLinearPredictionCalculateAdjustedExpectation <- function(SOSPrediction,DataInput , numberofvariables =18){
    
    m <- length(SOSPrediction$E_D)/numberofvariables
    n <- length(SOSPrediction$E_F)/numberofvariables
    
    tmpmatrix <- t(matrix( DataInput[ , (n+1):(n+m)], 1 , dim(DataInput)[1]*m ))
    
    E_D_F <- SOSPrediction$E_F + SOSPrediction$weightmatrix%*%(tmpmatrix - SOSPrediction$E_D)
    
    return(E_D_F)
  }
  For_ReshapetoMatrix <- function(E_D_F , n = 20 , numberofvariables = 18){
    output <- matrix(0 , n , numberofvariables)
    for(i in 1:numberofvariables){
      output[ , i] <- E_D_F[seq(i , numberofvariables*n , numberofvariables) ]  
    }
    return(output)
  }
  For_ReshapetoVector <- function(E_D_F){
    n <- dim(E_D_F)[1]
    numberofvariables <- dim(E_D_F)[2]
    output <- rep(0 , n*numberofvariables)
    for(i in 1:numberofvariables){
      output[ seq(i , numberofvariables*n , numberofvariables)] <- E_D_F[,i ]  
    }
    return(output)
  }
  
  For_PlotPrediction <- function(E_D_F , Ft){
    plot(E_D_F , ylim = c(min(c(E_D_F,Ft)),max(c(E_D_F,Ft))))
    points(Ft , col ='red')
    
  }
  
  
}

PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z209'] = "07/10/2016 10:46"
PatIndex2017$EndFirstNewAF[PatIndex2017$PseudoId == 'z209'] = "07/10/2016 11:29"

DataStructure <- array(NA , c(length(listAllPatients) , 25 , 400))

for(ii in 1:length(listAllPatients) ){
  
  PatientID <- listAllPatients[[ii]]
  MetaData <- DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID )
  
  
  if( file.exists(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData')) ){ 
    
    load(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData'))
    
    if(length(outputstruct) == 2){next}
    t <- FM_ExtractTimeFromHMOutput(outputstruct ) 
    SamsAnnotation <- BC_CreateAFAnnotationFomMetaData(t , MetaData)
    
    for(jj in 1:length(t)){
      if(!is.null(outputstruct[[jj]])){
        DataStructure[ii , 1,jj] <-  t[jj]
        if(nrow(outputstruct[[jj]][[1]])>=2 ){
          DataStructure[ii , 2:7,jj] <- apply(outputstruct[[jj]][[1]] , 2, median)}
          #DataStructure[ii , 2:7,jj] <- outputstruct[[jj]][[1]][sample(1:dim(outputstruct[[jj]][[1]])[1] , 1),] }
        if(!is.null(outputstruct[[jj]][[2]]$NonImplausibleSets) | !is.null(outputstruct[[jj]][[3]]$NonImplausibleSets) | !is.null(outputstruct[[jj]][[4]]$NonImplausibleSets) ){
          #tmpmatrix <- rbind(outputstruct[[jj]][[2]]$NonImplausibleSets  , outputstruct[[jj]][[3]]$NonImplausibleSets ,  outputstruct[[jj]][[4]]$NonImplausibleSets )
          #DataStructure[ii , 8:17,jj] <- tmpmatrix[sample(1:dim(tmpmatrix)[1] , 1) , ]
          DataStructure[ii , 8:17,jj] <- apply(rbind(outputstruct[[jj]][[2]]$NonImplausibleSets  , outputstruct[[jj]][[3]]$NonImplausibleSets ,  outputstruct[[jj]][[4]]$NonImplausibleSets ),2,median)
        }else{
          if(outputstruct[[length(outputstruct)]][jj,5] == 0 & outputstruct[[length(outputstruct)]][jj,2] == 1){
            DataStructure[ii , 8:17,jj] <- outputstruct[[jj]][[2]]$MinImplausiblePoint
          }
        }  
      }
    } 
  }
  
  
  DataStructure[ii , 18:22,1:dim(outputstruct[[length(outputstruct)]])[1]] <- outputstruct[[length(outputstruct)]]
  
  NonImplausiblePwaves <- DP_FindNARows(t(DataStructure[ii , 2:7, !is.na(DataStructure[ii , 18, ]) ]))
  NonImplausibleHeartRhythm <- DP_FindNARows(t(DataStructure[ii , 8:17, !is.na(DataStructure[ii , 18, ]) ]))
  
  DataStructure[ii , 23,1:length(t)] <-  NonImplausiblePwaves[1:length(t)]
  DataStructure[ii , 24,1:length(t)] <-  NonImplausibleHeartRhythm[1:length(t)]
  
  DataStructure[ii , 25,1:length(SamsAnnotation)] <- SamsAnnotation 
  
  
  if(sum(SamsAnnotation)>1){
    DataStructure[ii, ,c(1:which(SamsAnnotation ==1)[1]) ] <-  DataStructure[ii,  ,flipud(as.matrix(1:which(SamsAnnotation ==1)[1])) ]
    DataStructure[ii, ,c( (which(SamsAnnotation ==1)[1] +1):dim(DataStructure)[3] ) ] <- NA
    #DataStructure[ii, ,c(1:length(SamsAnnotation)) ] <-  DataStructure[ii,  ,flipud(as.matrix(1:length(SamsAnnotation))) ]
  }else{
    DataStructure[ii, ,c(1:length(SamsAnnotation)) ] <-  DataStructure[ii,  ,flipud(as.matrix(1:length(SamsAnnotation))) ]
  }
  
  DP_WaitBar(ii/length(listAllPatients))
}

DataStructure <- DataStructure[-which(listAllPatients == 'z209'),,]
listAllPatients <- listAllPatients[-which(listAllPatients == 'z209')]

AFLogical <- apply(as.matrix(listAllPatients) , 1 , function(X){ DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017 ,X ))})


DataStructureFilled <- DataStructure

for(ii in 1:length(listAllPatients)){
  DataStructureFilled[ii,,1:50] <- t(For_FillEmptyRecords(t(DataStructureFilled[ii , , 1:50]) ))
}


{
  load('AFDiscreteModelOutput.RData')
  
  DataForLogistic <- POM_SampledImputation(MasterData[MasterData$PseudoId %in% listAllPatients , ])
  model <- glm(formula = PreOpStepOutput[[2]]
               , family = binomial(link = "logit"),  data = DataForLogistic)
  summary(model)
  
  PreOpLogisticProbility <- predict(model , DataForLogistic , type = c('response'))
  
  
  DataForLogistic <- POM_SampledImputation(MasterData[MasterData$PseudoId %in% listAllPatients , ])
  model <- glm(formula = PostOpStepOutput[[2]]
               , family = binomial(link = "logit"),  data = DataForLogistic)
  summary(model)
  
  PostOpLogisticProbility <- predict(model , DataForLogistic , type = c('response'))}

### Get Patient data out of MasterData ###


### Analyse non-Implausible ####
rangeofpoints <- 2:30
BC_PlotCompareTwoHists(abs(cbind(c(DataStructure[AFLogical ,6 , rangeofpoints]) ,c(DataStructure[AFLogical ,7 , rangeofpoints]),
                                 c(DataStructure[AFLogical ,2 , rangeofpoints]) ,c(DataStructure[AFLogical ,3 , rangeofpoints]),
                                 c(DataStructure[AFLogical ,4 , rangeofpoints]) ,c(DataStructure[AFLogical ,5 , rangeofpoints]) )  )
                       ,abs(cbind(c(DataStructure[!AFLogical ,6 , rangeofpoints]) ,c(DataStructure[!AFLogical ,7 , rangeofpoints]) , 
                                  c(DataStructure[!AFLogical ,2 , rangeofpoints]) ,c(DataStructure[!AFLogical ,3 , rangeofpoints]) ,
                                  c(DataStructure[!AFLogical ,4 , rangeofpoints]) ,c(DataStructure[!AFLogical ,5 , rangeofpoints]) )))

BC_PlotCompareTwoHists(abs(cbind(c(DataStructure[AFLogical ,8 , rangeofpoints]) ,c(DataStructure[AFLogical ,9 , rangeofpoints]),
                                 c(DataStructure[AFLogical ,10 , rangeofpoints]) ,c(DataStructure[AFLogical ,11 , rangeofpoints]),
                                 c(DataStructure[AFLogical ,12 , rangeofpoints]) ,c(DataStructure[AFLogical ,13 , rangeofpoints]) )  )
                       ,abs(cbind(c(DataStructure[!AFLogical ,8 , rangeofpoints]) ,c(DataStructure[!AFLogical ,9 , rangeofpoints]) , 
                                  c(DataStructure[!AFLogical ,10 , rangeofpoints]) ,c(DataStructure[!AFLogical ,11 , rangeofpoints]) ,
                                  c(DataStructure[!AFLogical ,12 , rangeofpoints]) ,c(DataStructure[!AFLogical ,13 , rangeofpoints]) )))

BC_PlotCompareTwoHists(abs(cbind(c(DataStructure[AFLogical ,14 , rangeofpoints]) ,c(DataStructure[AFLogical ,15 , rangeofpoints]),
                                 c(DataStructure[AFLogical ,16 , rangeofpoints]) ,c(DataStructure[AFLogical ,17 , rangeofpoints]),
                                 c(DataStructure[AFLogical ,23 , rangeofpoints]) ,c(DataStructure[AFLogical ,24 , rangeofpoints]) )  )
                       ,abs(cbind(c(DataStructure[!AFLogical ,14 , rangeofpoints]) ,c(DataStructure[!AFLogical ,15 , rangeofpoints]) , 
                                  c(DataStructure[!AFLogical ,16 , rangeofpoints]) ,c(DataStructure[!AFLogical ,17 , rangeofpoints]) ,
                                  c(DataStructure[!AFLogical ,23 , rangeofpoints]) ,c(DataStructure[!AFLogical ,24 , rangeofpoints]) )))

{ nn <- 40
  mmatrixAF <- matrix(0,nn,18)
  CArrayAF <- array(0,c(nn,18,18))
  
  mmatrixNAF <- matrix(0,nn,18)
  CArrayNAF <- array(0,c(nn,18,18))
  
  nAF <- matrix(0,nn,1)
  nNAF <- matrix(0,nn,1) 
  
  MahalanobisMatrix <- matrix(0,nn,1) 
  NumberinAF <- matrix(0,nn,1)
  
  indicies <- 1:16
  
  ImMatrix <- array(0 , c(length(AFLogical) ,2,nn)  )
  
  for(jj in 1:nn){
    mmatrixAF[jj,] <- apply(abs(DataStructure[AFLogical ,  c(2:17 , 23 , 24) , jj]) , 2 , function(X){mean(X , na.rm = T)})  
    mmatrixNAF[jj,] <- apply(abs(DataStructure[!AFLogical , c(2:17 , 23 , 24) , jj]) , 2 , function(X){mean(X , na.rm = T)})  
    
    NumberinAF[jj,] <- sum(DataStructure[AFLogical ,25 , jj] , na.rm=T )
    
    tmp <- abs(DataStructure[AFLogical ,  c(2:17 , 23 , 24) , jj] )
    tmp <- tmp[DP_FindNARows(tmp),]
    nAF[jj,] <-dim(tmp)[1]
    CArrayAF[jj,,] <- cov(tmp)
    
    tmp <- abs(DataStructure[!AFLogical ,  c(2:17 , 23 , 24) , jj])
    tmp <- tmp[DP_FindNARows(tmp),]
    nNAF[jj,]<- dim(tmp)[1]
    CArrayNAF[jj,,] <- cov(tmp)
    
    mmatrixAF[ jj , 17:18] <- apply(abs(DataStructure[AFLogical ,  c(23 , 24) , jj]) , 2 , function(X){sum(X, na.rm = T)})/apply(DataStructure[AFLogical ,  c(23 , 24) , jj], 2,function(X){sum(!is.na(X))})
    mmatrixNAF[ jj , 17:18] <- apply(abs(DataStructure[!AFLogical ,  c(23 , 24) , jj]) , 2 , function(X){sum(X, na.rm = T)})/apply(DataStructure[!AFLogical ,  c(23 , 24) , jj], 2,function(X){sum(!is.na(X))})
    
    
    ImMatrix[, 1 ,jj ] <- apply(DataStructure[ ,  c(2:17) , jj] , 1 , function(X){mahalanobis(X , mmatrixAF[jj,indicies] , CArrayAF[jj,indicies,indicies]) })
    ImMatrix[ , 2 , jj] <- apply(DataStructure[ ,  c(2:17) , jj] , 1 , function(X){mahalanobis(X , mmatrixNAF[jj,indicies] , CArrayNAF[jj,indicies,indicies])})
    
    
    MahalanobisMatrix[jj] <- (mmatrixAF[jj,indicies] -   mmatrixNAF[jj,indicies] )%*%solve((CArrayNAF[jj,indicies,indicies]/ nNAF[jj,]) + (CArrayAF[jj,indicies,indicies]/ nAF[jj,]))%*%(mmatrixAF[jj,indicies] -   mmatrixNAF[jj,indicies] )
    
  }  
}

par(mfrow = c(3,1))
plot(sqrt(MahalanobisMatrix[2:nn]))
abline(h = sqrt(32))
title('Mahalanobis')
mmatrixAF[is.infinite(mmatrixAF[,17]),17] <- 0
plot(mmatrixAF[2:nn,17],col = 'red', ylim = c(0,1))
points(mmatrixNAF[2:nn,17],col = 'blue')
title('Proportion Non-Imaplusible P-waves')
mmatrixAF[is.infinite(mmatrixAF[,18]),18] <- 0
plot(mmatrixAF[2:nn,18],col = 'red', ylim = c(0,1))
points(mmatrixNAF[2:nn,18],col = 'blue')
title('Proportion Non-Imaplusible Heart Rhythms')


##### Full Model #####

# reshape matrix
{datasetfortest <-  matrix(0 , 50*sum(AFLogical) , length(c(2:17,23,24)) + 1 )
counter <- 1  

for( i in 1:dim(DataStructureFilled[AFLogical , c(2:17,23,24) , 1:50])[1] ){
  for( j in 3:dim(DataStructureFilled[AFLogical , c(2:17,23,24) , 1:50])[3] ){
    datasetfortest[counter , 1] <- 1     
    datasetfortest[counter , 2:dim(datasetfortest)[2] ] <- DataStructureFilled[which(AFLogical)[i] , c(2:17,23,24) , j]     
    counter <- counter + 1
  }
}  

datasetfortest2 <-  matrix(0 , 50*sum(!AFLogical) , length(c(2:17,23,24)) + 1 )
counter <- 1  

for( i in 1:dim(DataStructureFilled[!AFLogical , c(2:17,23,24) , 1:50])[1] ){
  for( j in 3:dim(DataStructureFilled[!AFLogical , c(2:17,23,24) , 1:50])[3] ){
    datasetfortest2[counter , 1] <- 0     
    datasetfortest2[counter , 2:dim(datasetfortest2)[2] ] <- DataStructureFilled[which(!AFLogical)[i] , c(2:17,23,24) , j]     
    counter <- counter + 1
  }
}

datasetfortest <- rbind(datasetfortest,datasetfortest2)
rm(datasetfortest2)


colnames(datasetfortest ) = c( "AFLogical",  "V2",  "V3",  "V4",  "V5",  "V6",  "V7",  "V8",  "V9",  "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17","V18","V19")

model <- glm(formula = AFLogical ~ V2 + V3 + V4 +V5+ V6 + V7 + V8 + V9  + V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 +V19 +
               I(V2^2)+I(V3^2)+I(V4^2)+ I(V5^2)  +I(V6^2)  +I(V8^2) +I(V9^2)+I(V11^2)+I(V12^2)+I(V13^2)+I(V14^2)+I(V15^2)+I(V16^2)+I(V17^2)
             , family = binomial(link = "logit"),  data =  data.frame(datasetfortest) )
stepoutput <- stepAIC(model)
}


model <- glm(formula = stepoutput$formula
             , family = binomial(link = "logit"),  data =  data.frame(datasetfortest))
summary(model)

AUCMatrix <- matrix(0 , nn , 1)
logisticProbabilityMatrix <- matrix(0 ,nn , length(listAllPatients) )
for(jj in 1:length(listAllPatients)){
  index <- jj
  
  DataForLogistic <-t( DataStructureFilled[index ,c(2:17,23,24) , 1:nn] )
  colnames(DataForLogistic ) <- c( "V2",  "V3",  "V4",  "V5",  "V6",  "V7",  "V8",  "V9",  "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17","V18","V19")
  DataForLogistic <- data.frame(DataForLogistic)
  
  logisticProbabilityMatrix[,jj] <- predict(model , DataForLogistic , type = c('response'))
  }
for(kk in 1:nn){
  AUCMatrix[kk] <- auc(AFLogical, logisticProbabilityMatrix[kk,] )
}

{ 
  x11()
  i <- 1
  tmp <- logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]
  tmp[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]])] <-   cummean(logisticProbabilityMatrix[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]),which(AFLogical)[i]] )
  plot(tmp , type = 'l' , ylim = c(0,1) , col = rgb(1,0,0,alpha = 0.25),xlab = 'time' , ylab = 'Output Logistic Regression')
  for(i in 2:length(which(!AFLogical)) ){
    tmp <- logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]
    tmp[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]])] <-   cummean(logisticProbabilityMatrix[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]),which(AFLogical)[i]] )
    
    lines(tmp , type = 'l' , ylim = c(0,1) , col = rgb(1,0,0,alpha = 0.25))
  }
  for(i in 1:length(which(!AFLogical)) ){
    tmp <- logisticProbabilityMatrix[nn:1,which(!AFLogical)[i]]
    tmp[!is.na(logisticProbabilityMatrix[nn:1,which(!AFLogical)[i]])] <-   cummean(logisticProbabilityMatrix[!is.na(logisticProbabilityMatrix[nn:1,which(!AFLogical)[i]]),which(!AFLogical)[i]] )
    
    lines(tmp , type = 'l' , ylim = c(0,1) , col = rgb(0,0,1,alpha = 0.1))}
  }


plot(apply(logisticProbabilityMatrix[1:nn,AFLogical] , 1, function(X){mean(X , na.rm = T)}) , ylim = c(0,0.5),xlab ='time', ylab = 'Mean Output Logistic Regression')
points(apply(logisticProbabilityMatrix[1:nn,!AFLogical] , 1, function(X){mean(X , na.rm = T)}) , col = 'red')

threshold = 0.1
{ 
  x11()
  i <- 1
  tmp <- logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]
  tmp[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]])] <-   cumsum(logisticProbabilityMatrix[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]),which(AFLogical)[i]] > 0.2 )
  plot(tmp , type = 'l' , ylim = c(0,40) , col = rgb(1,0,0,alpha = 0.25),xlab = 'time' , ylab = 'Output Logistic Regression')
  for(i in 2:length(which(!AFLogical)) ){
    tmp <- logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]
    tmp[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]])] <-   cumsum(logisticProbabilityMatrix[!is.na(logisticProbabilityMatrix[nn:1,which(AFLogical)[i]]),which(AFLogical)[i]] > 0.2 )
    
    lines(tmp , type = 'l' , ylim = c(0,1) , col = rgb(1,0,0,alpha = 0.25))
  }
  
  for(i in 1:length(which(!AFLogical)) ){
    tmp <- logisticProbabilityMatrix[nn:1,which(!AFLogical)[i]]
    tmp[!is.na(logisticProbabilityMatrix[nn:1,which(!AFLogical)[i]])] <-   cumsum(logisticProbabilityMatrix[!is.na(logisticProbabilityMatrix[nn:1,which(!AFLogical)[i]]),which(!AFLogical)[i]] > 0.2 )
    
    lines(tmp , type = 'l' , ylim = c(0,1) , col = rgb(0,0,1,alpha = 0.1))}
}

logisticProbabilityMatrix2 <- 0*logisticProbabilityMatrix
for(i in 1:length(listAllPatients)){
  tmp <- logisticProbabilityMatrix[nn:1,i]
  tmp[!is.na(logisticProbabilityMatrix[nn:1,i])] <-   cumsum(logisticProbabilityMatrix[nn:1,i][!is.na(logisticProbabilityMatrix[nn:1,i])] > 0.2 )
  
  logisticProbabilityMatrix2[ , i] <- tmp
}

maxvector <- apply(logisticProbabilityMatrix2[3:nn,] , 2 ,function(X){max(X,na.rm=T)})


Threhsold = 0

P <- sum(AFLogical)
N <- sum(!AFLogical)
TP <- sum(maxvector[AFLogical] > Threhsold)
TN <- sum(maxvector[!AFLogical] <= Threhsold)
FP <- sum(maxvector[!AFLogical & !is.infinite(AFLogical)] > Threhsold)
FN <- sum(maxvector[AFLogical] <= Threhsold)
  
Sensitivity <- TP / P
Specificity <- TN / N
PPV <- TP / (TP + FP)
NPV <- TN / (TN + FN)

