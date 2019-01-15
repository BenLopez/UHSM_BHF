CliqueList <- c(1 , 100 , 200 ,  400 , 600 , 1000)
ProbabilityList <- c(0.5 , 0.6 , 0.8 , 0.9 , 0.95 , 0.99)
MinutesThesholdList <- c( 2 , 4 , 6 , 8 , 10 )

PPV <- array(0 , c(length(CliqueList) , length(ProbabilityList) , length(MinutesThesholdList)))
NPV <- array(0 , c(length(CliqueList) , length(ProbabilityList) , length(MinutesThesholdList)))
Beat_Sen <- array(0 , c(length(CliqueList) , length(ProbabilityList) , length(MinutesThesholdList)))
Beat_Spec <- array(0 , c(length(CliqueList) , length(ProbabilityList) , length(MinutesThesholdList)))
Patient_Sen <- array(0 , c(length(CliqueList) , length(ProbabilityList) , length(MinutesThesholdList)))
Patient_Spec <- array(0 , c(length(CliqueList) , length(ProbabilityList) , length(MinutesThesholdList)))

if(!exists('Xstar')){
  Xstar = seq(0.5 ,1 , 0.01)}
if(!exists('PriorNonImplausibleSet')){
  PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95, 0.001 ) ,b = c(0.6  , 0.05 ) , numbersamples = 10000)
}

for(ii in 1){
  for(jj in 1){  
    for(kk in 1){

#BCParameters$TS_Likelihood_clique <- CliqueList[ii]
#BCParameters$ProbabilityThreshold <- ProbabilityList[jj]
#BCParameters$`Minute Threshold` <- MinutesThesholdList[kk]

{
  
Performancelist <- list( )
StartTimesAF <- list( )
PatientNames <- list( )
counter <- 1

for(Patinettotest in listAllPatients[1:length(listAllPatients)]){
  if(DP_checkRpeaksfilesprocessed(path , PatientsId = Patinettotest)){
    tmp <- DP_LoadRpeaksfile(path , Patinettotest )
  }
  if(DP_checkRpeaksfilesprocessed(path , PatientsId = Patinettotest) ){
    if((length(tmp$RRCombined$t) > 10000) ){
      source('BC_LoadDataandTestSinglePatient.R', echo = FALSE)
      StartTimesAF[[counter]] <- AFLocations
      if( Performance$P < 800 ){
          Performance$P <- 0
      }
      Performancelist[[counter]] <- Performance
      PatientNames[[counter]] <- Patinettotest
      DP_WaitBar(which(listAllPatients == Patinettotest)/length(listAllPatients))
      counter <- counter+1
    }
  }  
}
}

N <- 0
TN <- 0
P <- 0
TP <- 0

BeatSensitivities <- matrix(0 , length(Performancelist) , 1)
BeatSpecificties <-  matrix(0 , length(Performancelist) , 1)
CorrectIdentifiction <- matrix(0 , length(Performancelist) , 1)

CorrectIdentifiction <- matrix(0 , length(Performancelist) , 2)

for( i in 1:length(Performancelist) ){
  N  <- N + Performancelist[[i]]$N
  P  <- P + Performancelist[[i]]$P
  TN <- TN + Performancelist[[i]]$TN
  TP <- TP + Performancelist[[i]]$TP
  BeatSpecificties[i,1] <- Performancelist[[i]]$Specifictity
  BeatSensitivities[i,1] <- Performancelist[[i]]$Sensitvity
  
  if(Performancelist[[i]]$P == 0 & Performancelist[[i]]$FP == 0){
    CorrectIdentifiction[i,1] <- 1
  }
  if(Performancelist[[i]]$P > 0 & Performancelist[[i]]$TP > 0){
    CorrectIdentifiction[i,1] <- 1
    CorrectIdentifiction[i,2] <- 1
  }
  if(Performancelist[[i]]$P > 0 & Performancelist[[i]]$TP ==0 ){
    CorrectIdentifiction[i,1] <- 0
    CorrectIdentifiction[i,2] <- 1
  }
  
}

Beat_Sen <- TP/P
Beat_Spec <- TN/N

Patient_N   <-  sum( CorrectIdentifiction[, 2] == 0 )
Patient_P   <-  sum( CorrectIdentifiction[, 2] == 1 )
Patient_TN  <-  sum( CorrectIdentifiction[CorrectIdentifiction[, 2] == 0, 1] == 1 )
Patient_TP  <-  sum( CorrectIdentifiction[CorrectIdentifiction[, 2] == 1, 1] == 1 )

Patient_Sen  <- Patient_TP/Patient_P
Patient_Spec <- Patient_TN/Patient_N

print(Patient_Sen )
print(Patient_Spec)

FalseNegatives <-  PatientNames[which( ((CorrectIdentifiction[, 2] == 1)*(CorrectIdentifiction[, 1] == 0))==1)]
FalsePostives  <-  PatientNames[which( ((CorrectIdentifiction[, 2] == 0)*(CorrectIdentifiction[, 1] == 0))==1)]

#Postives <- PatientNames[which( CorrectIdentifiction[, 2] == 1 )]

PPV <- ( Priorprobabilities$A*Patient_Sen ) / (Priorprobabilities$A*Patient_Sen + Priorprobabilities$`A^c`*(1-Patient_Spec))
NPV <- ( Priorprobabilities$`A^c`*Patient_Spec ) / (Priorprobabilities$A*(1-Patient_Sen) + Priorprobabilities$`A^c`*Patient_Spec)
    }
  }
}



plot(PPV ,NPV )
title('NPV Against PPV')
abline(0.99,0)
abline( v = 0.99)
abline(0.98,0  , col = 'red')
abline(v = 0.98 , col = 'red')

plot(Patient_Sen ,Patient_Spec , xlab = 'Patient Wise Sensitivity' , ylab = 'Patient Wise Specificity' )
title('Sensitivity against Specifictity')
abline(0.99,0)
abline( v = 0.99)
abline(v = 0.98 , col = 'red')

plot(Beat_Sen , Beat_Spec , xlab = 'Beat Wise Sensitivity' , ylab = 'Beat Wise Specificity' , xlim = c(0,1))
title('Sensitivity against Specifictity')
abline(0.99,0)
abline( v = 0.99)
abline(0.999,0 , col = 'red')
abline(v = 0.98 , col = 'red')
