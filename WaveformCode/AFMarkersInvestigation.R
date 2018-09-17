{if(file.exists('CheckforDefaultsScript.R')){
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

listAllPatients <- DP_FilterPatients(listAllPatients,PatIndex2017, HowtoFilterops , path , FilestoProcess)

{ClinicalMarkers <- matrix(0,length(listAllPatients) , 5)
Opptypes <- matrix(0,length(listAllPatients) , 1)
AFLogical <- matrix(0,length(listAllPatients) , 1)

  
for( i in 1:length(listAllPatients) ){
sub_pat <- DP_ExtractPatientRecordforIndex(PatIndex2017 ,listAllPatients[i] )

ClinicalMarkers[i,1] <- sub_pat$Age[1]
ClinicalMarkers[i,2] <- sub_pat$Weight[1]
ClinicalMarkers[i,3] <- sub_pat$Height[1]

if(sub_pat$Gender[1] == "Female"){
ClinicalMarkers[i,4] <- 1
}
if(sub_pat$PreopRRT[1] == 'Yes'){
ClinicalMarkers[i , 5] <- 1
}
Opptypes[i , 1] <- sub_pat$ProcDetails[1]

if( DP_CheckIfAFPatient(sub_pat) ){
  AFLogical[i , 1] <- 1  
}
}
}

color <- matrix(rgb(1,0,0 , alpha = 0.1) , length(listAllPatients) , 1)
color[AFLogical == 1] <- rgb(0,0,1 , alpha = 0.1)

pairs(cbind(ClinicalMarkers , AFLogical) , col = color , pch = 16 , labels = c('Age' , 'Weight' , 'Height' , 'Gender' , 'PreopRRT' , 'AFLogical'))


m1 <- mean(ClinicalMarkers[AFLogical == 0 , 1])
sigma1 <- var(ClinicalMarkers[AFLogical == 0 , 1])


m2 <- mean(ClinicalMarkers[AFLogical == 1 , 1])
sigma2 <- var(ClinicalMarkers[AFLogical == 1 , 1])


dnorm(ClinicalMarkers[1,1], mean=m1, sd=sqrt(sigma1)) 
dnorm(ClinicalMarkers[1,1], mean=m2, sd=sqrt(sigma2)) 

pAF = 0.1

index <- 3
pAFgivenD <- (pAF*dnorm(ClinicalMarkers[,1], mean=m1, sd=sqrt(sigma1))) /( (((1-pAF)*dnorm(ClinicalMarkers[,1], mean=m2, sd=sqrt(sigma2)))) +  (((pAF)*dnorm(ClinicalMarkers[,1], mean=m1, sd=sqrt(sigma1)))))


