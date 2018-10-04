
###### Inference Functions ######
BC_LoadInAllDistributionSummaries <- function(path , listAllPatients , PatIndex2017 ){
  counter <- 1
  DataMatrix <- array(0 , dim = c(length(listAllPatients) , 50000 , 12) )
  PatientNames <- list()
  PatientDetails <- list()
  PatientAF <- matrix(0 , length(listAllPatients) , 1)
  
  print( 'Loading Data.' )
  for(ii in 1:length(listAllPatients)){
    if(DP_CheckDistributionSummariesExists(path , listAllPatients[[ii]]) == FALSE ){
      DP_WaitBar(ii/length(listAllPatients))  
      next}
    sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[[ii]])
    DistributionSummaries <- DP_LoadDistributionSummaries(path , listAllPatients[[ii]]) 
    if(length(DistributionSummaries[[1]]) > 50000){next}
    if(length(DistributionSummaries[[1]]) < 5000){next}
    
    PatientNames[[counter]] <- listAllPatients[[ii]]
    for(jj in 1:length(DistributionSummaries)){
      DataMatrix[counter , 1:min(length(DistributionSummaries[[jj]]) , 50000) , jj] <- DistributionSummaries[[jj]][1:min(length(DistributionSummaries[[jj]]) , 50000)]
      #DataMatrix[counter , 1:length(DistributionSummaries[[jj]])  , jj] <- DistributionSummaries[[jj]][1:(length(DistributionSummaries[[jj]]) )]
    }
    PatientDetails[[counter]] <- sub_pat
    if(DP_CheckIfAFPatient(sub_pat)){
      PatientAF[counter] <- 1 
    }
    counter <- counter +1
    DP_WaitBar(ii/length(listAllPatients))
  }
  print('Data Loaded. ')
  DataMatrix[is.na(DataMatrix)] <- 0  
  #PatientAF <-  PatientAF[apply(DataMatrix  , 1 , sum) != 0 , ]
  
  NAFMatrix <- DataMatrix[PatientAF == 0 , 1000:dim( DataMatrix )[2]  , ]
  DataMatrix <- DataMatrix[PatientAF == 1 , 1000:dim( DataMatrix )[2]  , ]
  PatientSample <- sample(1:dim(NAFMatrix)[1] ,250, replace = FALSE )
  
  NAFMatrix <- NAFMatrix[PatientSample , , ]
  AFMatrix <- DataMatrix
  
  AFPatients <- PatientNames[PatientAF == 1]
  NAFPatients <- PatientNames[PatientAF == 0]
  NAFPatients <- NAFPatients[PatientSample]
  
  AFPatientDetails  <- PatientDetails[PatientAF == 1]
  NAFPatientDetails <- PatientDetails[PatientAF == 0]
  
  output <- list(AFMatrix , NAFMatrix , AFPatients , NAFPatients)
  output <- setNames( output , c('AFPatientsDatabase' , 
                                 'NAFPatientsDatabase',
                                 'AFPatinetsNames',
                                 'NAFPatinetsNames') )
  return(output)
}
BC_LoadInAllCDFS <- function(path , listAllPatients , PatIndex2017 ){
  counter <- 1
  DataMatrix <- array(0 , dim = c(length(listAllPatients) , 30000 , 33) )
  PatientNames <- list()
  PatientDetails <- list()
  PatientAF <- matrix(0 , length(listAllPatients) , 1)

  Numberofvariables <- 32
    
  print( 'Loading Data.' )
  for(ii in 1:length(listAllPatients)){
    if(DP_CheckFileExists(path , listAllPatients[[ii]] , paste0(listAllPatients[[ii]] , '_CDFs' )) == FALSE ){
      DP_WaitBar(ii/length(listAllPatients))  
      next}
    sub_pat <- subset(PatIndex2017, PseudoId %in% listAllPatients[[ii]])
    CDFs <- DP_LoadFile(path , PatientsID = listAllPatients[[ii]] , Name = paste0(listAllPatients[[ii]] , '_CDFs' )) 
    if(dim(CDFs[[1]])[1] > 50000){next}
    if(dim(CDFs[[1]])[1] < 5000){next}
    
      PatientNames[[counter]] <- listAllPatients[[ii]]
      tmp1 <- CDFs$CDFs[(dim(CDFs[[1]])[1] - min(dim(CDFs[[1]])[1] , 30000) + 1):dim(CDFs[[1]])[1], ]
      tmp2 <-   CDFs$time[(dim(CDFs[[1]])[1] - min(dim(CDFs[[1]])[1] , 30000) + 1):dim(CDFs[[1]])[1] ]
      tmp1[is.na(tmp1)] <- 0
      
      tmp1[apply(tmp1 ==0 , 1 , sum) == Numberofvariables , ] = NaN
      tmp1 <- t(apply(tmp1 , 1 , function(X){DP_fixculmulativeprobs(X)}))
      
      DataMatrix[counter , 1:min(dim(CDFs[[1]])[1] , 30000) , 1:Numberofvariables ] <- tmp1
      DataMatrix[counter , 1:min(dim(CDFs[[1]])[1] , 30000) , (Numberofvariables + 1)] <- tmp2
      #DataMatrix[counter , 1:length(DistributionSummaries[[jj]])  , jj] <- DistributionSummaries[[jj]][1:(length(DistributionSummaries[[jj]]) )]
      PatientDetails[[counter]] <- sub_pat
    if(DP_CheckIfAFPatient(sub_pat)){
      PatientAF[counter] <- 1 
    }
    counter <- counter +1
    DP_WaitBar(ii/length(listAllPatients))
  }
  print('Data Loaded.')
  
  NAFMatrix <- DataMatrix[PatientAF == 0 , 1000:dim( DataMatrix )[2]  , ]
  DataMatrix <- DataMatrix[PatientAF == 1 , 1000:dim( DataMatrix )[2]  , ]
  PatientSample <- sample(1:dim(NAFMatrix)[1] ,250, replace = FALSE )
  NAFMatrix <- NAFMatrix[PatientSample , , ]
  AFMatrix <- DataMatrix
  
  AFPatients <- PatientNames[PatientAF == 1]
  NAFPatients <- PatientNames[PatientAF == 0]
  NAFPatients <- NAFPatients[PatientSample]
  
  AFPatientDetails  <- PatientDetails[PatientAF == 1]
  NAFPatientDetails <- PatientDetails[PatientAF == 0]
  
  output <- list(AFMatrix , NAFMatrix , AFPatients , NAFPatients)
  output <- setNames( output , c('AFPatientsDatabase' , 
                                 'NAFPatientsDatabase',
                                 'AFPatinetsNames',
                                 'NAFPatinetsNames') )
  return(output)
}
BC_CreateAFandNOAFDataStructure <- function(DataBaseMaster , PatIndex2017){
  Numberofvariables <- dim(DataBaseMaster$AFPatientsDatabase)[3] - 1
  timevariable <- dim(DataBaseMaster$AFPatientsDatabase)[3]
  AFVector <- matrix(0 , 1 , Numberofvariables)
  NAFVector<- matrix(0 , 1 , Numberofvariables)
  mmAF <- matrix(0 ,  size(DataBaseMaster$AFPatientsDatabase)[1],Numberofvariables)
  mmNAF <- matrix(0 ,  size(DataBaseMaster$NAFPatientsDatabase)[1],Numberofvariables)
  
  for(ii in 1:size(DataBaseMaster$AFPatientsDatabase)[1]){
    sub_pat <- DP_ExtractPatientRecordforIndex(PatIndex2017 , DataBaseMaster$AFPatinetsNames[[ii]])
    
    if(Numberofvariables < 20){
    tmp <- DataBaseMaster$AFPatientsDatabase[ii ,DataBaseMaster$AFPatientsDatabase[ii , ,1] != 0 ,]
    }else{
    tmp <- DataBaseMaster$AFPatientsDatabase[ii , DataBaseMaster$AFPatientsDatabase[ii , ,Numberofvariables] != 0,]
    tmp <- tmp[DP_FindNumzeroRows(tmp) < Numberofvariables , ]
    }
    
    if(nrow(tmp)< 5000){
      mmAF[ii,] <- apply( tmp[,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)} )
      next}else{
      mmAF[ii,] <- apply(tmp[1:5000,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
      }
    for(jj in 1:(Numberofvariables)){
      tmp[ , jj] <- tmp[ , jj] - mmAF[ii,jj] 
    }
    
    AFVector <- rbind(AFVector , tmp[(tmp[ , timevariable] > as.numeric(DP_StripTime(sub_pat$ConfirmedFirstNewAF)))*(tmp[ , timevariable] < as.numeric(DP_StripTime(sub_pat$EndFirstNewAF))) ==1 ,1:Numberofvariables]  )
    NAFVector <- rbind(NAFVector , tmp[tmp[ , timevariable] <  DP_StripTime(sub_pat$ConfirmedFirstNewAF) ,1:Numberofvariables])
    DP_WaitBar(ii/(size(DataBaseMaster$NAFPatientsDatabase)[1] + size(DataBaseMaster$AFPatientsDatabase)[1]))
  }
  
  for(ii in 1:size(DataBaseMaster$NAFPatientsDatabase)[1]){
    if(Numberofvariables < 20){
      tmp <- DataBaseMaster$NAFPatientsDatabase[ii ,DataBaseMaster$NAFPatientsDatabase[ii , ,1] != 0 ,]
    }else{
      tmp <- DataBaseMaster$NAFPatientsDatabase[ii , DataBaseMaster$NAFPatientsDatabase[ii , ,Numberofvariables] != 0,]
      tmp <- tmp[DP_FindNumzeroRows(tmp) < Numberofvariables , ]
    }   
    
    if(nrow(tmp) < 5000){
      mmNAF[ii,] <- apply(tmp[,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
      next}
      mmNAF[ii,] <- apply(tmp[1:5000,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
    for(jj in 1:(size(tmp)[2] - 1)){
      tmp[ , jj] <- tmp[ , jj] - mmNAF[ii,jj]
    }
    NAFVector <- rbind(NAFVector , tmp[,1:Numberofvariables] )
    DP_WaitBar((ii + + size(DataBaseMaster$AFPatientsDatabase)[1])/(size(DataBaseMaster$NAFPatientsDatabase)[1] + size(DataBaseMaster$AFPatientsDatabase)[1]))
  }
  return(setNames( list( AFVector , NAFVector , mmAF , mmNAF ) , c('AF , NAF' , 'M_AF' , 'M_NAF') ) )
}  
BC_Removezerovariancevariables <- function(DataBase){
  for(i in 1:length(DataBase)){
    DataBase[[i]] <- DataBase[[i]][ , which(DP_FindZeroVarianceRows(DataBase[[i]]) == FALSE)]
  }
  return(DataBase)
}
BC_CreateAFPreAFandNOAFDataStructure <- function(DataBaseMaster){
Numberofvariables <- dim(DataBaseMaster$AFPatientsDatabase)[2] - 1
timevariable <- dim(DataBaseMaster$AFPatientsDatabase)[2]
  
  AFVector <- matrix(0 , 1 , Numberofvariables)
  PreAFVector <- matrix(0 , 1 , Numberofvariables)
  NAFVector<- matrix(0 , 1 , Numberofvariables)
  mmAF <- matrix(0 ,  size(DataBaseMaster$AFPatientsDatabase)[1],Numberofvariables)
  mmNAF <- matrix(0 ,  size(DataBaseMaster$NAFPatientsDatabase)[1],Numberofvariables)
  sub_pat <- DP_ExtractPatientRecordforIndex(PatIndex2017 , DataBaseMaster$AFPatinetsNames[[ii]])
  
  for(ii in 1:size(DataBaseMaster$AFPatientsDatabase)[1]){
    tmp <- DataBaseMaster$AFPatientsDatabase[ii ,DataBaseMaster$AFPatientsDatabase[ii , ,1] != 0 ,]
    if(nrow(tmp)< 5000){
      mmAF[ii,] <- apply( tmp[,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)} )
      next}
    mmAF[ii,] <- apply(tmp[1:5000,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
    for(jj in 1:(size(tmp)[2] -1)){
      tmp[ , jj] <- tmp[ , jj] - mmAF[ii,jj] 
    }  
    AFVector <- rbind(AFVector , tmp[(tmp[ , timevariable] > as.numeric(DP_StripTime(sub_pat$ConfirmedFirstNewAF)))*(tmp[ , timevariable] < as.numeric(DP_StripTime(sub_pat$EndFirstNewAF))) ==1 ,1:Numberofvariables]  )
    PreAFVector <- rbind(PreAFVector , tmp[tmp[ , timevariable] <  DP_StripTime(sub_pat$ConfirmedFirstNewAF) ,1:Numberofvariables])
    DP_WaitBar(ii/(size(DataBaseMaster$NAFPatientsDatabase)[1] + size(DataBaseMaster$AFPatientsDatabase)[1]))
  }
  
  for(ii in 1:size(DataBaseMaster$NAFPatientsDatabase)[1]){
    tmp <- DataBaseMaster$NAFPatientsDatabase[ii ,DataBaseMaster$NAFPatientsDatabase[ii , ,1] != 0 ,]
    if(nrow(tmp) < 5000){
      mmNAF[ii,] <- apply(tmp[,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
      next}
    mmNAF[ii,] <- apply(tmp[1:5000,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
    for(jj in 1:(size(tmp)[2] - 1)){
      tmp[ , jj] <- tmp[ , jj] - mmNAF[ii,jj]
    }
    NAFVector <- rbind(NAFVector , tmp[,1:Numberofvariables] )
    DP_WaitBar((ii + + size(DataBaseMaster$AFPatientsDatabase)[1])/(size(DataBaseMaster$NAFPatientsDatabase)[1] + size(DataBaseMaster$AFPatientsDatabase)[1]))
  }  
  return(setNames( list( AFVector , NAFVector , PreAFVector, mmAF , mmNAF ) , c('AF , NAF' , 'PreAFVector' , 'M_AF' , 'M_NAF') ) )
} 
BC_EstimateLocalDensitiesMVN <- function(  DataBase ){
ClassMeans <- matrix(0 , length(DataBase) -2 , size(DataBase[[1]])[2] )
ClassCovarianceMatrices <- array(0 , c(length(DataBase) -2 , size(DataBase[[1]])[2] , size(DataBase[[1]])[2]) )

for(i in 1:(length(DataBase) -2) ){
  DataBase[[i]] <- DataBase[[i]][ , DP_FindZeroVarianceRows(DataBase[[i]]) == FALSE]
  ClassMeans[i , ] <- apply(DP_RemoveNaRows(DataBase[[i]]),2 , function(X){mean(X[!is.na(X)])})
  ClassCovarianceMatrices[i , ,] <- cov(DP_RemoveNaRows(DataBase[[i]]))
}  
return(setNames(list(ClassMeans , ClassCovarianceMatrices) , c('mu' , 'Sigma')))
}
BC_EstimateGlobalDensitiesMVN <- function( DataBase ){
  numberclasses <- length(DataBase) -2 
  ClassMeans <- matrix(0 , numberclasses , size(DataBase[[1]])[2] )
  ClassCovarianceMatrices <- array(0 , c(numberclasses , size(DataBase[[1]])[2] , size(DataBase[[1]])[2]) )
  
  for(i in 1:2 ){
    DataBase[[numberclasses+i]] <- DataBase[[numberclasses+i]][ , DP_FindZeroVarianceRows(DataBase[[numberclasses+i]]) == FALSE]
    ClassMeans[i , ] <- apply(DataBase[[numberclasses+i]][apply(is.na(DataBase[[numberclasses+i]]) , 1 , sum) == 0 , ],2 , function(X){mean(X[!is.na(X)])})
    ClassCovarianceMatrices[i , ,] <- cov(DataBase[[numberclasses+i]][apply(is.na(DataBase[[numberclasses+i]]) , 1 , sum) == 0 , ])
  }  
  return(setNames(list(ClassMeans , ClassCovarianceMatrices) , c('mu' , 'Sigma')))
}
BC_EstimateLocalDensitiesGMM <- function(  DataBase , numberofcomponents = 5){
  
  GMMStructureforallclasses <- list()
  for(i in 1:(length(DataBase) -2) ){
  GMMStructureforallclasses[[i]] <- densityMclust( DP_RemoveNaRows(DataBase[[i]]) , G = numberofcomponents )
  }  
  return(GMMStructureforallclasses)
}
BC_EstimateGlobalDensitiesGMM <- function( DataBase , numberofcomponents = 5){
  numberclasses <- length(DataBase) -2 
  
  GMMStructureforallclasses <- list()
  for(i in 1:2 ){
    GMMStructureforallclasses[[i]] <- densityMclust( DP_RemoveNaRows(DataBase[[numberclasses + i]]) , G = numberofcomponents )
  }  
  return(GMMStructureforallclasses)
}
BC_EstimateGlobalEffect <- function( Z ){
  if( size(Z)[1] > 5000){
  mmAF <- apply(Z[1:5000,],2,function(X){mean(X , na.rm = TRUE)})
  }else{
    mmAF <- apply(Z,2,function(X){mean(X , na.rm = TRUE)})
  }
  for(jj in 1:(size(Z)[2] -1)){
    Z[ , jj] <- Z[ , jj] - mmAF[jj] 
  }
  output <- setNames( list(mmAF , Z ) , c('M' , 'W'))
}
BC_EstimateGlobalandLocalParametersDisSum <- function( DistributionSummaries ){
  Numberofvariables <- (length(DistributionSummaries) -1)
  t <- DistributionSummaries[[length(DistributionSummaries)]]
  mm <- matrix(0 , 1 , Numberofvariables)
  W <- matrix(0 , length(DistributionSummaries[[1]]) , Numberofvariables)
  for(i in 1:Numberofvariables){
    W[,i] <- DistributionSummaries[[i]]
    }
  W[W[,8] == 1 , 9] <- 0
  W[W[,8] == 1,10] <- 0
  
  if(length(DistributionSummaries[[1]])< 5000){
    mm <- apply( W[,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)} )
    next}else{
    mm <- apply(W[1:5000,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
    }
  for(jj in 1:(size(W)[2] -1)){
    W[ , jj] <- W[ , jj] - mm[jj] 
  }
  return(setNames(list(W , mm) , c('W' , 'M') ))
}
BC_EstimateGlobalandLocalParametersCDFs <- function( CDFs ){
  Numberofvariables <- dim(CDFs$CDFs )[2]
  t <- CDFs$time
  mm <- matrix(0 , 1 , Numberofvariables)
  W <- CDFs$CDFs
  W[is.na(W)] <- 0
  W[apply(W ==0 , 1 , sum) == Numberofvariables , ] = NaN
  W <- t(apply(W , 1 , function(X){DP_fixculmulativeprobs(X)}))
  Numberofvariables <- dim(W)[2]
  
  if(length(t)< 5000){
    mm <- apply( W[,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)} )
    next}else{
    mm <- apply(W[1000:6000,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
    }
  for(jj in 1:(size(W)[2] -1)){
    W[ , jj] <- W[ , jj] - mm[jj] 
  }
  mm <- mm[2:31]
  W <- W[ , 2:31]
  
  return(setNames(list(W , mm) , c('W' , 'M') ))
}
BC_CalulateImplausabilty <- function( Z , mu , Sigma ){
  Sigma <- Sigma + 0.000000001*diag(dim(Sigma)[1])
  Implausability <- mahalanobis(x = Z , center = mu , cov = Sigma )
  return(Implausability)
}
BC_CalulateImplausabiltySecondOrderStatistics <- function(DataBase , SecondOrderStruct){
num_classes <-  length(DataBase) - 2
output <- list()
output[[1]] <- matrix(0 , num_classes , 1 )
output[[2]] <- matrix(0 , num_classes , 1 )
output <- setNames( output , c('mu' , 'Sigma') )
for(i in 1:num_classes){
  Im <- BC_CalulateImplausabilty(DP_RemoveNaRows(DataBase[[i]]) , mu = SecondOrderStruct$mu[i,], Sigma = (SecondOrderStruct$Sigma[i , ,]) )
  output[[1]][i] <- mean( Im[!is.na(Im)] )
  output[[2]][i] <- var( Im[!is.na(Im)] )
}
return(output)
}
BC_CalulateRegularisedImplausabilty <- function( Z , mu , sigma , muIm , VarIm ){
  Implausability <- (BC_CalulateImplausabilty(Z , mu , sigma) -  muIm)/sqrt(VarIm)
  return(Implausability)
}  
BC_TestingCalulateRegularisedImplausabilty<-function( TestData , SecondOrderStruct , ImSecondOrderStruct){
  nclasses <- size(SecondOrderStruct$mu)[1]
  RegularisedImplausability <- matrix(0 , dim(TestData)[1] ,  nclasses)
  for(i in 1:nclasses ){
    RegularisedImplausability[ , i] <- BC_CalulateRegularisedImplausabilty(Z = TestData ,
                                                                           mu =SecondOrderStruct$mu[i,] ,
                                                                           sigma = (SecondOrderStruct$Sigma[i,,]),
                                                                           muIm = ImSecondOrderStruct$mu[i],
                                                                           VarIm = ImSecondOrderStruct$Sigma[i] ) 
  }
return(abs(RegularisedImplausability))  
}
BC_CaculateBeatWiseProbabiltiesAFDetection<-function(Priorprobabilities){

  Priorprobabilities[[5]] <- Priorprobabilities[[1]]*Priorprobabilities[[3]] + Priorprobabilities[[2]]*Priorprobabilities[[4]] 
  Priorprobabilities[[6]] <- Priorprobabilities[[2]] + Priorprobabilities[[1]]*(1 - Priorprobabilities[[3]])
  return(Priorprobabilities)
}
BC_CalulateDenistiesGMM <- function(W , LocalDistributionStruct){
  num_classes <- length(LocalDistributionStruct)
  f_i <- matrix(0 , size(W)[1] , length(LocalDistributionStruct))
  for(i in 1:num_classes){
    f_i[ , i] <- predict( LocalDistributionStruct[[i]] , W , what = c('dens'))  
  } 
  return(f_i)  
}
BC_BayesianBeliefUpdateGMM <- function(W , LocalDistributionStruct ,  Probabilities , n =1 ){
num_classes <- length(LocalDistributionStruct)
NALogical <- DP_FindNARows(W)
if(n > 1){
f_i <- BC_CaluluateCulmulativeLikelihood(f_i  = BC_CalulateDenistiesGMM(W= DP_RemoveNaRows(W) , LocalDistributionStruct = LocalDistributionStruct) , n =n  )
}else{
f_i <- BC_CalulateDenistiesGMM(DP_RemoveNaRows(W) , LocalDistributionStruct)  
}

posteriorprobabilities <- matrix(0 , size(W)[1] , num_classes)
for(i in 1:num_classes){
  posteriorprobabilities[NALogical , i] <- Probabilities[[2*num_classes + i]]*f_i[ , i]  
}
posteriorprobabilities <- apply( posteriorprobabilities ,2, function(X){X/apply(posteriorprobabilities , 1 , sum)})
return(posteriorprobabilities)
}
BC_CalulateCulmulativeImplausability <- function(Implausability , n = 100){
  Implausability <- apply(Implausability , 2 , function(X){rollmax(x=X , k =n, na.pad = TRUE , align = c("center") )})
  return(Implausability)
}
BC_CaluluateCulmulativeLikelihood <- function( f_i , n = 5 , weight = 250 ){
for(i in 1:size(f_i)[2]){
  f_i[f_i[,i] <= 0,i] = 1e-10
  f_i[,i] <- exp(log(f_i[,i]) + (n/weight)*smth( log(f_i[,i]) , method = 'sma' , n=n ))
}
  return( f_i )  
}
BC_GlobalCalulateImplausabiltySecondOrderStatistics <- function(DataBase , SecondOrderStruct){
  num_classes <-  2
  output <- list()
  output[[1]] <- matrix(0 , num_classes , 1 )
  output[[2]] <- matrix(0 , num_classes , 1 )
  output <- setNames( output , c('mu' , 'Sigma') )
  for(i in 1:num_classes){
    Im <- BC_CalulateImplausabilty(DP_RemoveNaRows(DataBase[[2+i]]) , mu = SecondOrderStruct$mu[i,], Sigma =SecondOrderStruct$Sigma[i , ,] )
    output[[1]][i] <- mean( Im[!is.na(Im)] )
    output[[2]][i] <- var( Im[!is.na(Im)] )
  }
  return(output)
}
BC_GlobalTestingCalulateRegularisedImplausabilty<-function( M , GlobalSecondOrderStruct , GlobalImSecondOrderStruct){
  nclasses <- size(GlobalSecondOrderStruct$mu)[1]
  RegularisedImplausability <- matrix(0 , 1 ,  nclasses)
  for(i in 1:nclasses ){
    RegularisedImplausability[ , i] <- BC_CalulateRegularisedImplausabilty(Z = M ,
                                                                           mu =GlobalSecondOrderStruct$mu[i,] ,
                                                                           sigma =GlobalSecondOrderStruct$Sigma[i,,],
                                                                           muIm = GlobalImSecondOrderStruct$mu[i],
                                                                           VarIm = GlobalImSecondOrderStruct$Sigma[i] ) 
  }
  return(abs(RegularisedImplausability))  
}
BC_GlobalBayesianBeliefUpdateMVN <- function( M , GlobalSecondOrderStruct , Priorprobabilities ){
  nclasses <- size(GlobalSecondOrderStruct$mu)[1]
  output <- Priorprobabilities
  f_i <- matrix(0 , 1 , nclasses)
  for(i in 1:nclasses){
    f_i[ , i] <- exp(mvnpdf(x = t(as.matrix(M)) ,  mu = GlobalSecondOrderStruct$mu[i,] , Sigma = GlobalSecondOrderStruct$Sigma[i,,]))
  }
  output$A <- Priorprobabilities$A*f_i[ , 1] / (Priorprobabilities$A*f_i[ , 1] + Priorprobabilities$`A^c`*f_i[ , 2])
  output$`A^c` <- 1 - output$A
  output <- BC_CaculateBeatWiseProbabiltiesAFDetection(output) 
  return(output)
}
BC_GlocalBayesianBeliefUpdateGMM <- function(M , GlobalDistributionStruct ,  Priorprobabilities , n =1 ){
  num_classes <- length(GlobalDistributionStruct)
  f_i <- BC_CalulateDenistiesGMM(W = M , LocalDistributionStruct = GlobalDistributionStruct) 
  output <- Priorprobabilities
  output$A <- Priorprobabilities$A*f_i[ , 1] / (Priorprobabilities$A*f_i[ , 1] + Priorprobabilities$`A^c`*f_i[ , 2])
  output$`A^c` <- 1 - output$A
  output <- BC_CaculateBeatWiseProbabiltiesAFDetection(output) 
  return(output)
}
BC_CreateAFAnnotationFomMetaData <- function(t , MetaData){
  output <- matrix(0 , length(t) , 1)
  output[ (t>DP_StripTime(MetaData$ConfirmedFirstNewAF))*(t<DP_StripTime(MetaData$EndFirstNewAF)) == 1] <- 1
  return(output)
}
BC_CreateAnnotationFromInference <- function(t , AFLocations){
  output <- matrix(0 , length(t) , 1)
  for( i in 1:length(AFLocations$Start) ){
    output[ (t>AFLocations$Start[i])*(t<AFLocations$End[i]) == 1] <- 1 
  }
  return(output)
}
BC_CalulatePerformance <- function(AnnotatedAFMetaData ,  AnnotatedAFInference , BadDataLocations){
  P <- sum((BadDataLocations == 0)*AnnotatedAFMetaData)
  N <- sum((BadDataLocations == 0)*(AnnotatedAFMetaData ==0))
  Total <- N+P - sum(BadDataLocations)
  
  TP <- sum((BadDataLocations == 0)*(AnnotatedAFMetaData==1)*(AnnotatedAFInference == 1))
  TN <- sum((BadDataLocations == 0)*(AnnotatedAFMetaData==0 )*(AnnotatedAFInference==0))
  FP <- sum((BadDataLocations == 0)*(AnnotatedAFMetaData == 0)*(AnnotatedAFInference == 1))
  FN <- sum((BadDataLocations == 0)*(AnnotatedAFMetaData == 1)*(AnnotatedAFInference ==0))
  
  Sen  <- TP / P 
  Spec <- TN / N
  PPV  <- TP /(TP + FP) 
  NPV  <- TN /(TN + FN) 
  return( data.frame(Sensitvity = Sen , Specifictity = Spec , PPV = PPV , NPV = NPV , P = P , N = N , Total=Total , TP = TP , TN = TN , FP = FP , FN=FN) )
}
BC_CheckForTimeGaps <- function(AFLocations , BadDataLocations , t , ECGs , beatlims = c(45 , 300) , readinglim = 0.65){
  output <- AFLocations
  for(i in 1:length(AFLocations[ , 1])){
    timeinAF <- difftime( AFLocations$End[i] , AFLocations$Start[i] , units = 'secs' )
    
    minbeats <- as.numeric((beatlims[1]/60)*timeinAF)
    maxbeats <- as.numeric((beatlims[2]/60)*timeinAF)
    
    beats <- sum(((t > AFLocations$Start[i])*(t < AFLocations$End[i])) == 1 ) 
    
    if(beats < minbeats || beats > maxbeats){
      output <- output[-i , ]
      warning('Implausible Heart Rate')
      next
    }
    
    expectedmeasurments <- as.numeric(timeinAF)/as.numeric(abs(ECGs$ECGI$Date[1] - ECGs$ECGI$Date[2]))
    numberofmeasurements = c(0,0,0)
    numberofmeasurements[1] <- sum( ((ECGs$ECGI$Date > AFLocations$Start[i])*(ECGs$ECGI$Date < AFLocations$End[i])) == 1 )
    numberofmeasurements[2] <- sum( ((ECGs$ECGII$Date > AFLocations$Start[i])*(ECGs$ECGII$Date < AFLocations$End[i])) == 1 )
    numberofmeasurements[3] <- sum( ((ECGs$ECGIII$Date > AFLocations$Start[i])*(ECGs$ECGIII$Date < AFLocations$End[i])) == 1 )
    
    if(sum((numberofmeasurements / expectedmeasurments) < readinglim) > 1 ){
      BadDataLocations <<- rbind( BadDataLocations , output[i , ])
      output <- output[-i , ]
      print('Missing Data AF Period Removed')
      next
    }
    
    
  }
  
  return(output)
}
BC_ExtractValidationPriors <- function(Priorprobabilities , DataBaseMaster , DataBase){
  
  DataSetPriorProbabilities <- Priorprobabilities
  DataSetPriorProbabilities$A <- length(DataBaseMaster$AFPatinetsNames)/(length(DataBaseMaster$NAFPatinetsNames)+length(DataBaseMaster$AFPatinetsNames))
  DataSetPriorProbabilities$`A^c` <- 1 - DataSetPriorProbabilities$A
  DataSetPriorProbabilities$`B|A` <- dim(DataBase[[1]])[1] /( dim(DataBase[[1]])[1] + dim(DataBase[[2]])[1])
  DataSetPriorProbabilities <- BC_CaculateBeatWiseProbabiltiesAFDetection( DataSetPriorProbabilities )
  return(DataSetPriorProbabilities)
}
BC_EstimateGlobalandLocalParameters<-function( Z ){
  Numberofvariables <- dim(Z)[2]
  mm <- matrix(0 , 1 , Numberofvariables)
  W <- Z
  if(dim(Z)[1]< 5000){
    mm <- apply( W[,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)} )}else{
      mm <- apply(W[1:5000,1:Numberofvariables],2,function(X){mean(X , na.rm = TRUE)})
    }
  for(jj in 1:(size(W)[2] -1)){
    W[ , jj] <- W[ , jj] - mm[jj] 
  }
  return(setNames(list(W , mm) , c('W' , 'M') ))
}
BC_CreateCalibrationStructure <- function(GlobalProbCalibrationStruct ,   BinWidth = 0.1){
  ProbabililtyBins <- seq(0,1,BinWidth)
  ProbabililtyBinStruct <-  matrix(0 , (length(ProbabililtyBins) -1) , 2)
  EstimatorError <-  matrix(0 , (length(ProbabililtyBins) -1) , 1)
  
  for(i in 1:(dim(ProbabililtyBinStruct)[1])){
    tmplogical <- ((GlobalProbCalibrationStruct[,1] >= ProbabililtyBins[i])*(GlobalProbCalibrationStruct[,1] <= ProbabililtyBins[i+1])) == 1
    ProbabililtyBinStruct[i,1] <- mean( GlobalProbCalibrationStruct[tmplogical , 2] )
    ProbabililtyBinStruct[i,2] <- mean( GlobalProbCalibrationStruct[tmplogical , 1] )
    if(is.na(mean( GlobalProbCalibrationStruct[tmplogical , 1] )) ){ProbabililtyBinStruct[i,2] <- (ProbabililtyBins[i] + BinWidth/2)  }
    EstimatorError[i,] <- BC_CalulateBernoulliEstimatorUncertainty(p = ProbabililtyBinStruct[i,1], length(GlobalProbCalibrationStruct[tmplogical , 2]) )
  }
  
  return(data.frame(x = ProbabililtyBinStruct[,2], y =ProbabililtyBinStruct[,1] , sd = sqrt(EstimatorError) ) )
}
BC_CalulateBernoulliEstimatorUncertainty <- function( p , n){
  return( (p*(1 - p)/n) + 1/(n^2) )
}
BC_CreateDefaultmclustStruct <- function(){
  MclustDistributionStruct <- list( parameters = list() )  
  MclustDistributionStruct$parameters[[1]] <- 1
  MclustDistributionStruct$parameters[[2]] <- 1
  MclustDistributionStruct$parameters[[3]] <- list()
  MclustDistributionStruct$parameters <- setNames(MclustDistributionStruct$parameters , c('pro' , 'mean' , 'variance'))  
  MclustDistributionStruct$parameters$variance <- setNames(list(1) , c('sigma'))
  return(MclustDistributionStruct)
}
BC_CreateProbCalibrationStruct <- function(PosteriorProbabilities ,  alpha , numberofvalidationsamples){
  return( cbind(as.matrix(PosteriorProbabilities) , 
                rbind(matrix(1 , round(alpha[1]* numberofvalidationsamples)  , 1) , 
                      matrix(0 , round(alpha[2]* numberofvalidationsamples)  , 1) ) ) )
}
BC_PredictGMMDensity <- function(MclustDistributionStruct , x){
  N = length( MclustDistributionStruct$parameters$pro)
  f_i <- matrix(0 , N , 1)
  for(ii in 1:N){
    f_i[ii,] <- MclustDistributionStruct$parameters$pro[ii]*exp(DF_mvnpdf(x , mu = MclustDistributionStruct$parameters$mean[ , ii] , Sigma =  MclustDistributionStruct$parameters$variance$sigma[ , , ii]))   
  }
  return(sum(f_i))
}
BC_CleanProbCalibrationOutput <- function(ProbabiliticCalibrationOutput){
  ProbabiliticCalibrationOutput$sd[is.na(ProbabiliticCalibrationOutput$y)] <- 1
  ProbabiliticCalibrationOutput$y[is.na(ProbabiliticCalibrationOutput$y)] <- ProbabiliticCalibrationOutput$x[is.na(ProbabiliticCalibrationOutput$y)]
  ProbabiliticCalibrationOutput$sd[1] <- ProbabiliticCalibrationOutput$sd[1] + 0.005
  ProbabiliticCalibrationOutput$sd[length(ProbabiliticCalibrationOutput$sd)] <- ProbabiliticCalibrationOutput$sd[length(ProbabiliticCalibrationOutput$sd)] + 0.005
  return(ProbabiliticCalibrationOutput)
}

###### Plotting Functions ######

BC_PlotCreateggImplausabilities <- function(TimeVector , RegIm){
if(dim(RegIm)[2] == 2){
  Implot <- BC_PlotCreateggImplausabilitiesBase(TimeVector , RegIm ) +
    geom_line( aes( y = Im.1/max(RegIm[!is.na(RegIm)]) , colour = "blue") ) + 
    geom_line( aes( y = Im.2/max(RegIm[!is.na(RegIm)]) , colour = "red") ) +
    scale_color_discrete(name = "Normalised Implausibilities", labels = c("AFib", "Not AFib"))
}
if(dim(RegIm)[2] == 3){
    Implot <-  BC_PlotCreateggImplausabilitiesBase(TimeVector , RegIm ) +
      geom_line( aes( y = Im.1/max(RegIm[!is.na(RegIm)]) , colour = "blue") ) + 
      geom_line( aes( y = Im.2/max(RegIm[!is.na(RegIm)]) , colour = "red") )  +
      geom_line( aes( y = Im.2/max(RegIm[!is.na(RegIm)]) , colour = "black") )  +
      scale_color_discrete( name = "Implausibilities", labels = c("AFib", "Not AFib" , "Pre AF") )
}
return(Implot)  
}
BC_PlotCreateggImplausabilitiesBase <- function(TimeVector , RegIm){
  Implot <- ggplot(data.frame(t = TimeVector , Im = RegIm ), aes(t) )+
    ylab('Implausibility') +
    xlab('Time Relative to FirstNewAF')+
    ggtitle('Normalised Implausibility') + 
    geom_hline(yintercept = min(3/max(RegIm[!is.na(RegIm)])  , 1)) +
    ylim(0, 1)
  return(Implot)
}  
BC_PlotCreateggPosteriorProbabilities <- function(TimeVector , PosteriorProbabilities){
  if(dim(RegIm)[2] == 2){
    Posteriorplot <- BC_PlotCreateggPosteriorProbabilitiesBase(TimeVector = TimeVector  , PosteriorProbabilities =  PosteriorProbabilities ) +
      geom_line( aes( y = pp.1 , colour = "blue") ) + 
      geom_line( aes( y = pp.2 , colour = "red") ) +
      scale_color_discrete(name = "Probabilities", labels = c("AFib", "Not AFib"))
  }
  if(dim(RegIm)[2] == 3){
    Posteriorplot <-  BC_PlotCreateggPosteriorProbabilitiesBase(TimeVector , PosteriorProbabilities ) +
      geom_line( aes( y = pp.1 , colour = "blue") ) + 
      geom_line( aes( y = pp.2 , colour = "red") )  +
      geom_line( aes( y = pp.2 , colour = "black") )  +
      scale_color_discrete( name = "Probabilities", labels = c("AFib", "Not AFib" , "Pre AF") )
  }
  return(Posteriorplot)  
}
BC_PlotCreateggPosteriorProbabilitiesBase <- function(TimeVector , PosteriorProbabilities){
  Posteriorplot <- ggplot(data.frame(t = TimeVector , pp = PosteriorProbabilities ), aes(t) )+
    ylab('P(Class | Data)') +
    xlab('Time Relative to FirstNewAF')+
    ggtitle('Posterior probabilities') + 
    ylim(0, 1)
  return(Posteriorplot)
}  
BC_PlotCreateRRTimesPlots <- function( RPeaksStruct , color = 'blue' ,  MetaData = DP_CreateDummyMetaData() ){
  t <- RPeaksStruct$RRCombined$t[seq(1 , length(RPeaksStruct$RRCombined$t) , 3)]   
  RR <- RPeaksStruct$RRCombined$RR[seq(1 ,  length(RPeaksStruct$RRCombined$t)  , 3)] 
  
  plotstruct <- ggplot(data.frame( t = t , RR=RR ), aes(t , RR) ) +
    geom_point( colour = color ,  alpha=0.03 ) +
    ggtitle('RR times') +
    xlab( "t" ) +
    ylab( "RR" ) + 
   coord_cartesian(ylim = c(0, 1.5))
  plotstruct <- BC_PlotAddAFLines(plotstruct = plotstruct , MetaData = MetaData)  
  return(plotstruct)
}
BC_PlotAddAFLines<- function(plotstruct , MetaData , color = 'purple'){
  plotstruct <- plotstruct +
    geom_vline( xintercept = as.numeric( as.POSIXct( DP_StripTime(MetaData$EndFirstNewAF)) ) , linetype="dashed" , color = color ) +
    geom_vline( xintercept = as.numeric( as.POSIXct( DP_StripTime(MetaData$ConfirmedFirstNewAF)) )  , color = color )  
  return(plotstruct)
}
BC_PlotCreateRATimesPlots <- function(RPeaksStruct , ECGs ,  MetaData = DP_CreateDummyMetaData() ){
  t <- RPeaksStruct$RRCombined$t
  tt <- list( )
  RA <- list( )
  options( digits.secs=3 )
  for(i in 1:3){
    tt[[i]] <- ECGs[[i]]$Date[round(as.numeric(ECGs[[i]]$Date) , 3) %in% round(as.numeric(t) , 3)]
    RA[[i]] <- ECGs[[i]]$Value[round(as.numeric(ECGs[[i]]$Date) , 3) %in% round(as.numeric(t) , 3)]
    if(i<1){
      if(length(RA[[i]])/length(RA[[1]]) < 0.5)
        tt[[i]] <- ECGs[[i]]$Date[ ECGs[[2]]$Date %in% tt[[1]] ]
      RA[[i]] <- ECGs[[i]]$Value[ ECGs[[2]]$Date %in% tt[[1]] ]
    }
}
 

  plotstruct <- ggplot(data.frame(t=tt[[1]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] , RA=RA[[1]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] ) , aes(t , RA)) +
    geom_point(data = data.frame(t=tt[[1]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] , RA=RA[[1]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] ) , colour = "blue" ,  alpha=0.03 ) +
    geom_point(data = data.frame(t=tt[[2]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] , RA=RA[[2]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] ) , colour = "red"  ,  alpha=0.03 ) +
    geom_point(data = data.frame(t=tt[[3]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] , RA=RA[[3]][seq(1 , length(RPeaksStruct$RRCombined$t) , 3)] ) , colour = "green",  alpha=0.03 ) +
    ggtitle('R-amplitudes') +
    xlab("t") +
    ylab("RA") + coord_cartesian(ylim = c(-200, 200)) 
  plotstruct <- BC_PlotAddAFLines(plotstruct = plotstruct , MetaData )
  return(plotstruct)
}
BC_PlotCreateECGPlots <- function(RPeaksStruct ,  ECGStruct  , timestart = ECGStruct$Date[10000] , ECGindex = 1, timeindex = 10 , color = 'blue'){
  options( digits.secs=3 )
  t <-  ECGStruct$Date[(ECGStruct$Date > timestart )*(ECGStruct$Date < (timestart+timeindex) ) == 1]
  t1 <- RPeaksStruct$RRCombined$t[(RPeaksStruct$RRCombined$t > timestart )*(RPeaksStruct$RRCombined$t < (timestart+timeindex) ) == 1]
  t2 <- RPeaksStruct[[ECGindex]]$t[(RPeaksStruct[[ECGindex]]$t > timestart )*(RPeaksStruct[[ECGindex]]$t < (timestart+timeindex) ) == 1]
  
  f_t <- ECGStruct$Value[(ECGStruct$Date > timestart )*(ECGStruct$Date < (timestart+timeindex) ) == 1]
  RA1 <- f_t[t %in% t1]
  t1 <- t[t %in% t1]
  RA2 <- RPeaksStruct[[ECGindex]]$RA[(RPeaksStruct[[ECGindex]]$t > timestart )*(RPeaksStruct[[ECGindex]]$t < (timestart+timeindex) ) == 1]
  
  plotstruct <- ggplot(data.frame(t = t , f_t = f_t ) , aes(t , f_t))+
    geom_line(colour= color) + 
    ylab('Hz') + 
    xlab('t') +
    geom_point(data = data.frame(t = t1 , RA = RA1) , aes(t , RA) , colour = 'red' ) +  
    geom_point(data = data.frame(t = t2 , RA = RA2) , aes(t , RA) , colour = 'black'  )
  return(plotstruct)  
}
BC_PlotECGAlignAxis <- function(plotstruct , timerange){
  plotstruct <- plotstruct +   xlim(timerange[1] , timerange[2] ) 
  return(plotstruct)
}
BC_PlotECGAddTitle <- function(plotstruct , MetaData = ' ' , time = ' ' , ECGindex = 1){
  if(ECGindex == 1){
    ECGTitle <- 'ECGI'
    plotstruct <- plotstruct + ggtitle(paste0(ECGTitle , MetaData$PseudoId , '  ' ,  time))
  }
  if(ECGindex == 2){
    ECGTitle <- 'ECGII' 
    plotstruct <- plotstruct + ggtitle(paste0(ECGTitle  , '  ' ,  time))
  }
  if(ECGindex == 3){
    ECGTitle <- 'ECGIII'
    plotstruct <- plotstruct + ggtitle(paste0(ECGTitle, '  ' ,  time))
  }
  return(plotstruct)
}  
BC_PlotAddViewingRegionLines <- function(plotstruct , timerange){
  plotstruct <- plotstruct + 
    geom_vline( xintercept = as.numeric(timerange[1]) , linetype="dashed" , color = "black" ) +
    geom_vline( xintercept = as.numeric(timerange[2]) , linetype="dashed" , color = "black" )
  return(plotstruct)
}
BC_PlotCreatePosteriorBeliefPlot <- function(t , Implausability , PosteriorProbabilities){
  plotstruct <- BC_PlotCreateggImplausabilities( TimeVector = t , RegIm = Implausability )
  
  if(dim(Implausability)[2] == 2){
    plotstruct <- plotstruct +
      geom_line( data = data.frame(t = t , pp = PosteriorProbabilities), aes(x=t ,  y = pp.1 ) , color = 'red') + 
      geom_line( data = data.frame(t = t , pp = PosteriorProbabilities), aes(x=t ,  y = pp.2 ) , color = 'blue') +
      scale_color_discrete(name = " ", labels = c("AFib", "Not AFib")) + 
      ylim(0,1) +
      scale_y_continuous(sec.axis = sec_axis(~.*1, name = "P(Class | Data)")) + 
      ggtitle('Updated Beliefs')
    
  }
  
  if(dim(Implausability)[2] == 3){
    plotstruct <-  plotstruct +
      geom_line(data = data.frame(t = t , pp = PosteriorProbabilities), aes( y = pp.1 , colour = "blue"), color = 'blue') + 
      geom_line(data = data.frame(t = t , pp = PosteriorProbabilities), aes( y = pp.2 , colour = "red"), color = 'red')  +
      geom_line(data = data.frame(t = t , pp = PosteriorProbabilities), aes( y = pp.2 , colour = "black"),color = 'clack' )  +
      scale_color_discrete( name = " ", labels = c("AFib", "Not AFib" , "Pre AF") )  +
      scale_y_continuous(sec.axis = sec_axis(~.*1, name = "P(Class | Data)"))+ 
      ylim(0,1) +
      ggtitle('Updated Beliefs')
  }
  return(plotstruct)
}
BC_plotAddColouredRegions <- function(plotstruct , Locations  , fillcolor = 'pink'){
  if(length(Locations$Start) > 0){
    for( i in ( 1:length(Locations$Start) ) ){
      plotstruct <- plotstruct + annotate("rect" , xmin = Locations$Start[i], xmax = Locations$End[i], ymin = -1000, ymax= 1000 , fill = fillcolor , alpha = 0.25)
    }
  }
  return(plotstruct)  
}
BC_PlotCreateProbabilityCalibrationPlot <- function(Data){
  GlobalProbabilityCalibrationPlot <- ggplot(Data  , aes(x = x , y = y)) +
    geom_point( color = 'blue') +
    geom_errorbar(aes(ymin =  y - 2*sd , ymax = y + 2*sd ) , width = .01 ) +
    geom_line(aes(x = x , y = x))+
    xlab('Predicted Probability') +
    ylab('Estimated Probability')
  return(GlobalProbabilityCalibrationPlot)
}
BC_plotValidateDensityEstimationMarginalHistograms <- function(A , B , C , D ){
  x11(20,14)
  par(mfrow = c(2 , 5))
  for(variabletoview in c(1:10)){
    
    # Create marginal histograms  
    tmp  <- hist(A[ , variabletoview], col=rgb(0,0,1,alpha = 0.5) ,
                 main = paste0(AFD_CreateDistributionSummaryNames()[variabletoview] , ' Histogram') , freq = FALSE  , plot = FALSE)
    tmp2 <- hist(B[ , variabletoview], col=rgb(1,0,0,alpha =0.5) , freq = FALSE , 
                 breaks = c(min(B[!is.na(B[ , variabletoview]) , variabletoview] ) , tmp$breaks, max(B[!is.na(B[ , variabletoview]) , variabletoview]) ), plot = FALSE)
    
    tmp3 <- hist(C[ , variabletoview], col=rgb(1,0,0,alpha =0.5) , freq = FALSE , 
                 breaks = c(min(C[!is.na(C[ , variabletoview]) , variabletoview] ) , tmp$breaks, max(C[!is.na(C[ , variabletoview]) , variabletoview]) ), plot = FALSE)
    tmp4 <- hist(D[ , variabletoview], col=rgb(1,0,0,alpha =0.5) , freq = FALSE , 
                 breaks = c(min(D[!is.na(D[ , variabletoview]) , variabletoview] ) , tmp$breaks, max(D[!is.na(D[ , variabletoview]) , variabletoview]) ), plot = FALSE)
    
    plot(tmp2$density , tmp2$density - tmp4$density , xlab = 'index' , ylab = 'Residual Between Expected and Predicted Density' , col = 'red')
    points(tmp$density[2:(length(tmp$density) - 1)],tmp$density[2:(length(tmp$density) - 1)] - tmp3$density[3:(length(tmp2$density) - 2)] , col = rgb(0,0,1))
    title(paste0(round(sum( abs(tmp$density[2:(length(tmp$density) - 1)] - tmp3$density[3:(length(tmp2$density) - 2)]) ) , 3) , ' ',
                 round(sum( abs(tmp2$density - tmp4$density)) , 3)) )
  }
}
BC_PlotPairsFromTwoVariables <- function( X , Y , alpha = 0.01 ){
  x11(20 , 14)
  pairs( rbind(X , Y) ,  col = rbind(rep(rgb(1,0,0 , alpha = alpha) , size(X)[1]) , rep(rgb(0,0,1 , alpha = alpha) , size(Y)[1])) , pch = 16) 
}
BC_PlotPairs<- function(X , alpha = 0.01 ){
  x11(20 , 14)
  pairs( X , col = rgb(1 , 0 , 0 , alpha = 0.01) , pch =16) 
}
BC_PlotCompareTwoHists <- function( X , Y ){
  x11(20 , 14)
  par( mfrow = c(2 , ceiling(dim(X)[2]/2)) )    
  
  for(variabletoview in 1:dim(X)[2]){
    
    tmphist <- hist(rbind(X[ , variabletoview] , Y[ , variabletoview]) , breaks = 30 , plot = FALSE) 
    hist(X[ , variabletoview]
         , col=rgb(1,0,0,alpha =0.5) 
         , freq = FALSE 
         , breaks = tmphist$breaks
         , main = paste0('Variable ' , variabletoview)
         , xlabel = paste0('Variable = ' , variabletoview))
    hist(Y[ , variabletoview]
         , col = rgb(0,0,1,alpha =0.5) 
         , freq = FALSE 
         , breaks = tmphist$breaks
         , add = T)   
  }
}  
BC_CreateSecondOrderSpecificationHistImp <- function(Trainingset , Validationset ,  numberofsamples = 1000 , numberofrepitions = 100){
  Immatrix <- matrix(0 , numberofrepitions ,  size(Trainingset)[2]  )
  
  for( ii in 1:numberofrepitions ){
    
    sampleindexes <- sample( 1:dim(Validationset)[1] , numberofsamples )
    validationsample <- Validationset[sampleindexes , ]
    
    Immatrix[ii , 1:( size(Trainingset)[2] )] <- KDE_CalulateHistImplausability( Trainingset , validationsample) 
    
  }
  output <- list( mu = apply(Immatrix , 2 , mean) , Sigma = cov(Immatrix) , Inv_Sigma = solve(cov(Immatrix)) )
  return(output)
}

##### Sampling functions #####

BC_SampleGMM <- function(MclustDistributionStruct , numberofsamples){
  N = length( MclustDistributionStruct$parameters$pro)
  Sample1 <- t(rmultinom(numberofsamples, size = 1 , MclustDistributionStruct$parameters$pro)) 
  
  output <- matrix(0 , numberofsamples  , dim(MclustDistributionStruct$parameters$mean)[1])
  
  for( i in 1:N ){
    if(sum(Sample1[,i] == 1) ==0){next}
    output[Sample1[,i] == 1 , ] <- mvrnorm(n = sum(Sample1[,i]) ,
                                           mu = MclustDistributionStruct$parameters$mean[ , i], 
                                           Sigma = MclustDistributionStruct$parameters$variance$sigma[ , , i])
    
  }
  
  return(output)  
}
BC_SampleMGMM <- function( GMMmodelbyPatient , numberofsamples, alpha ){
  
  N = length( GMMmodelbyPatient)
  Sample1 <- t(rmultinom(numberofsamples, size = 1 , alpha)) 
  
  output <- matrix(0 , numberofsamples  , dim(GMMmodelbyPatient[[1]]$parameters$mean)[1])
  
  for( i in 1:N ){
    if(sum(Sample1[,i] == 1) ==0){next}
    output[Sample1[,i] == 1 , ] <- BC_SampleGMM(GMMmodelbyPatient[[i]] , numberofsamples = sum(Sample1[,i]))
    
  }
  
  return(output) 
  
}
