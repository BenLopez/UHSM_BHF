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
    plot( Ft, ylim = c(min(c(E_D_F,Ft)),max(c(E_D_F,Ft))), col ='red')
    points(E_D_F )
    
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
    #DataStructure[ii, ,c(1:which(SamsAnnotation ==1)[1]) ] <-  DataStructure[ii,  ,flipud(as.matrix(1:which(SamsAnnotation ==1)[1])) ]
    #DataStructure[ii, ,c( (which(SamsAnnotation ==1)[1] +1):dim(DataStructure)[3] ) ] <- NA
    DataStructure[ii, ,c(1:length(SamsAnnotation)) ] <-  DataStructure[ii,  ,flipud(as.matrix(1:length(SamsAnnotation))) ]
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

# load in preop and post op forecasting models.
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


##### Bayes Linear use m to predict n #####

m <- 10
n <- 30
DataStructureForBayesLinear <- abs(DataStructureFilled[,c(2:17, 23 , 24),1:(n+m) ])

##### Calculate a second order specification from the Data #####


SOSPrediction <- For_BuildSOSFromData(DataStructureForBayesLinear , m , n)
V_D_F <- SOSPrediction$V_F - SOSPrediction$weightmatrix%*%t(SOSPrediction$cov_FD)

  
# vectorise 2nd and 3rd dimension for Bayes linear


# Plots to look at prediction for AF Patients 

{
index <- 38
SOSPrediction <- For_BuildSOSFromData(DataStructureForBayesLinear[-which(AFLogical)[index],,] , m , n)
  
print(listAllPatients[which(AFLogical)[index]])
E_D_F <- For_BayesLinearPredictionCalculateAdjustedExpectation(SOSPrediction , DataInput = DataStructureForBayesLinear[which(AFLogical)[index],,] )
Ft = For_ReshapetoVector(t(DataStructureForBayesLinear[which(AFLogical)[index],,1:(n+m)]))

E_D_F <- For_ReshapetoMatrix(E_D_F)
Ft <- For_ReshapetoMatrix(Ft, n = 30 )


x11()
par(mfrow = c(2 , 3))
For_PlotPrediction(E_D_F[,1],Ft[,1])
For_PlotPrediction(E_D_F[,2],Ft[,2])
For_PlotPrediction(E_D_F[,3],Ft[,3])
For_PlotPrediction(E_D_F[,4],Ft[,4])
For_PlotPrediction(E_D_F[,5],Ft[,5])
For_PlotPrediction(E_D_F[,6],Ft[,6])


#x11()
#par(mfrow = c(2 , 5))
#For_PlotPrediction(E_D_F[,7],Ft[,7])
#For_PlotPrediction(E_D_F[,8],Ft[,8])
#For_PlotPrediction(E_D_F[,9],Ft[,9])
#For_PlotPrediction(E_D_F[,10],Ft[,10])
#For_PlotPrediction(E_D_F[,11],Ft[,11])
#For_PlotPrediction(E_D_F[,12],Ft[,12])
#For_PlotPrediction(E_D_F[,13],Ft[,13])
#For_PlotPrediction(E_D_F[,14],Ft[,14])
#For_PlotPrediction(E_D_F[,15],Ft[,15])
#For_PlotPrediction(E_D_F[,16],Ft[,16])

BC_PlotPairsFromTwoVariables( PriorNonImplausibleSet[1:1000,], E_D_F[1:10,1:6] , alpha = 0.1  ) 

print(min(apply(abs(PriorNonImplausibleSet - E_D_F[10,1:6]) / sqrt(For_ReshapetoMatrix(diag(V_D_F))[10,1:6] ) , 1 , max) ))
print(min(apply(abs(matrix(0,10,2) - E_D_F[1:10,5:6]) / sqrt(For_ReshapetoMatrix(diag(V_D_F))[1:10,5:6] ) , 1 , max) ))

}


# Plots to look at prediction for NAF Patients 


{
  index <- 10
  SOSPrediction <- For_BuildSOSFromData(DataStructureForBayesLinear[-which(!AFLogical)[index],,] , m , n)
  print(listAllPatients[which(!AFLogical)[index]])
  E_D_F <- For_BayesLinearPredictionCalculateAdjustedExpectation(SOSPrediction , DataInput = DataStructureForBayesLinear[which(!AFLogical)[index],,] )
  Ft = For_ReshapetoVector(t(DataStructureForBayesLinear[which(!AFLogical)[index],,1:(n+m)]))
  
  E_D_F <- For_ReshapetoMatrix(E_D_F)
  Ft <- For_ReshapetoMatrix(Ft, n = 30 )
  
  
  
  x11()
  par(mfrow = c(2 , 3))
  For_PlotPrediction(E_D_F[,1],Ft[,1])
  For_PlotPrediction(E_D_F[,2],Ft[,2])
  For_PlotPrediction(E_D_F[,3],Ft[,3])
  For_PlotPrediction(E_D_F[,4],Ft[,4])
  For_PlotPrediction(E_D_F[,5],Ft[,5])
  For_PlotPrediction(E_D_F[,6],Ft[,6])
  
  
  #x11()
  #par(mfrow = c(2 , 5))
  #For_PlotPrediction(E_D_F[,7],Ft[,7])
  #For_PlotPrediction(E_D_F[,8],Ft[,8])
  #For_PlotPrediction(E_D_F[,9],Ft[,9])
  #For_PlotPrediction(E_D_F[,10],Ft[,10])
  #For_PlotPrediction(E_D_F[,11],Ft[,11])
  #For_PlotPrediction(E_D_F[,12],Ft[,12])
  #For_PlotPrediction(E_D_F[,13],Ft[,13])
  #For_PlotPrediction(E_D_F[,14],Ft[,14])
  #For_PlotPrediction(E_D_F[,15],Ft[,15])
  #For_PlotPrediction(E_D_F[,16],Ft[,16])
  
  BC_PlotPairsFromTwoVariables( PriorNonImplausibleSet[1:1000,], E_D_F[,1:6] , alpha = 0.1  ) 
  print(min(apply(abs(PriorNonImplausibleSet - E_D_F[1,1:6]) / sqrt(For_ReshapetoMatrix(diag(V_D_F))[1,1:6] ) , 1 , max) ))
  print(min(apply(abs(matrix(0,10,2) - E_D_F[1:10,5:6]) / sqrt(For_ReshapetoMatrix(diag(V_D_F))[1:10,5:6] ) , 1 , max) ))
  
  }
