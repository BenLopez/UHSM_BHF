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
  #listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

plot(HaemIndex2017$`z1007-1-1-1`$TimeSeriesData$time , HaemIndex2017$`z1007-1-1-1`$TimeSeriesData$Platelets , xlab = 'time' , ylab = 'WBC')
abline(v = as.numeric(DP_StripTime(PatIndex2017$ConfirmedFirstNewAF[PatIndex2017$PseudoId == 'z1007']))  )
abline(v = as.numeric(DP_StripTime(PatIndex2017$FirstNewAF[PatIndex2017$PseudoId == 'z1007']))   , col ='red')
title('Platelets Time Series')



###### Pre Operation Analysis ######

HaemDataMatrix <- matrix(0 , length(HaemIndex2017) , 5)
AFLogical <- matrix(0 , length(HaemIndex2017) , 1)
PatientNames <- matrix(0 , length(HaemIndex2017) , 1)

for( i in 1:dim(HaemDataMatrix)[1] ){

  sub_pat <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = DP_ExtractPatientIDFromNewPatinetID(names(HaemIndex2017[i])))
    
  if(nrow(sub_pat) == 0){next}
  if(!is.na(sub_pat$Pre_OperativeHeartRhythm[1])){
  if(sub_pat$Pre_OperativeHeartRhythm[1] == "Atrial fibrillation/flutter" ){next}
  }
  PatientNames[i] = sub_pat$PseudoId[1]
  HaemDataMatrix[i,] = as.numeric(HaemIndex2017[[i]]$MetaData[ , c(4,5,6,8,11 ) ])
  
  if(!is.na(sub_pat$ConfirmedFirstNewAF) || !is.na(sub_pat$FirstNewAF) ){
    AFLogical[i] = 1
  }
  if(!is.na(sub_pat$ConfirmedFirstNewAF) ){
  if(sub_pat$ConfirmedFirstNewAF[1] == "CNAF"){
    AFLogical[i] = 0
  }  
  }
}

rownames(HaemDataMatrix) = PatientNames
colnames(HaemDataMatrix) = names(HaemIndex2017[[1]]$MetaData[ , c(4,5,6,8,11 ) ])

AFLogical = AFLogical[PatientNames != 0]
HaemDataMatrix = HaemDataMatrix[PatientNames != 0, ]
PatientNames = PatientNames[PatientNames !=0]

BC_PlotCompareTwoHists(HaemDataMatrix[AFLogical == 1, ],HaemDataMatrix[AFLogical == 0, ])
BC_PlotPairsFromTwoVariables(HaemDataMatrix[AFLogical == 1, ],HaemDataMatrix[sample(which(AFLogical == 0 ) , sum(AFLogical)) , ] , alpha = 0.1 )

PreOpHaem <- data.frame(cbind(HaemDataMatrix , AFLogical))
PreOpHaem$AFLogical = as.factor(PreOpHaem$AFLogical)

model <- (glm(formula = AFLogical ~ PreopHb + PreopWBC + PreopPLT    ,family=binomial(link='logit') , data=PreOpHaem ))
summary( model )

BC_PlotCompareSingleHists(model$fitted.values[model$y == 0 ],model$fitted.values[model$y == 1 ] )


###### Imediatley post Operation Analysis ######

{PostopHaemDataMatrix <- matrix(0 , length(HaemIndex2017) , 7)
AFLogical <- matrix(0 , length(HaemIndex2017) , 1)
PatientNames <- matrix(0 , length(HaemIndex2017) , 1)

for( i in 1:dim(PostopHaemDataMatrix)[1] ){
  
  sub_pat <- DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017 , PatientCode = DP_ExtractPatientIDFromNewPatinetID(names(HaemIndex2017[i])))
  
  if(nrow(sub_pat) == 0){next}
  if(!is.na(sub_pat$Pre_OperativeHeartRhythm[1])){
    if(sub_pat$Pre_OperativeHeartRhythm[1] == "Atrial fibrillation/flutter" ){next}
  }
  PatientNames[i] = sub_pat$PseudoId[1]
  if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[1 , 2])){
  PostopHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(5,6,8,9) ])
  PostopHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[1 , c(2,3 , 4) ])
  }
  if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[2 , 2])){
    PostopHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(2,3 , 4) ])
    PostopHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[1 , c(5,6,8,9) ])
  }
  if(is.na(HaemIndex2017[[i]]$TimeSeriesData[2 , 2]) & is.na(HaemIndex2017[[i]]$TimeSeriesData[1 , 2])){
    PostopHaemDataMatrix[i, ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(2,3 , 4 , 5,6,8,9) ])
    
  }
  
  if(!is.na(sub_pat$ConfirmedFirstNewAF) || !is.na(sub_pat$FirstNewAF) ){
    AFLogical[i] = 1
  }
  if(!is.na(sub_pat$ConfirmedFirstNewAF) ){
    if(sub_pat$ConfirmedFirstNewAF[1] == "CNAF"){
      AFLogical[i] = 0
    }  
  }
}


rownames(PostopHaemDataMatrix) = PatientNames
colnames(PostopHaemDataMatrix) = names(HaemIndex2017[[1]]$TimeSeriesData[1 , c(2,3 , 4, 5,6,8,9) ])

AFLogical = AFLogical[PatientNames != 0]
PostopHaemDataMatrix = PostopHaemDataMatrix[PatientNames != 0, ]
PatientNames = PatientNames[PatientNames !=0]

# rmove outlier
PostopHaemDataMatrix[725 ,1] <- median(PostopHaemDataMatrix[ ,1] , na.rm=T)
PostopHaemDataMatrix[182 ,2] <- median(PostopHaemDataMatrix[ ,2] , na.rm=T)

PostopHaemDataMatrix[ , 1] = log(PostopHaemDataMatrix[,1])
PostopHaemDataMatrix[ , 4] = log(PostopHaemDataMatrix[,4])
PostopHaemDataMatrix[ , 5] = log(PostopHaemDataMatrix[,5])
PostopHaemDataMatrix[ , 6] = log(PostopHaemDataMatrix[,6])

}


x11()
par(mfrow = c(1 , 6))
BC_PlotCompareSingleHists(PostopHaemDataMatrix[AFLogical == 0 ,1],PostopHaemDataMatrix[AFLogical == 1 ,1])
BC_PlotCompareSingleHists(PostopHaemDataMatrix[AFLogical == 0 ,2],PostopHaemDataMatrix[AFLogical == 1 ,2])
BC_PlotCompareSingleHists(PostopHaemDataMatrix[AFLogical == 0 ,3],PostopHaemDataMatrix[AFLogical == 1 ,3])
BC_PlotCompareSingleHists(PostopHaemDataMatrix[AFLogical == 0 ,4],PostopHaemDataMatrix[AFLogical == 1 ,4])
BC_PlotCompareSingleHists(PostopHaemDataMatrix[AFLogical == 0 ,5],PostopHaemDataMatrix[AFLogical == 1 ,5])
BC_PlotCompareSingleHists(PostopHaemDataMatrix[AFLogical == 0 ,6],PostopHaemDataMatrix[AFLogical == 1 ,6])


BC_PlotPairsFromTwoVariables(PostopHaemDataMatrix[AFLogical == 1, 1:3],PostopHaemDataMatrix[sample(which(AFLogical == 0 ) , sum(AFLogical)) , 1:3] , alpha = 0.5 )


PostOpHaem <- data.frame(cbind(PostopHaemDataMatrix , AFLogical))
PostOpHaem$AFLogical = as.factor(PostOpHaem$AFLogical)

model <- (glm(formula = AFLogical ~  WBC + Hb  + Platelets  ,family=binomial(link='logit') , data=PostOpHaem ))
summary( model )
BC_PlotCompareSingleHists(model$fitted.values[model$y == 0 ],model$fitted.values[model$y == 1 ] )

##### Pre AFib #####


{
PreAFHaemDataMatrix <- matrix(0 , length(HaemIndex2017) , 8)
AFLogical <- matrix(0 , length(HaemIndex2017) , 1)
PatientNames <- matrix(0 , length(HaemIndex2017) , 1)

for( i in 1:dim(HaemDataMatrix)[1] ){
  
  sub_pat <-  subset(PatIndex2017, NewPseudoId %in% names(HaemIndex2017[i]))
    
  if(nrow(sub_pat) == 0){next}
  if(!is.na(sub_pat$Pre_OperativeHeartRhythm[1])){
    if(sub_pat$Pre_OperativeHeartRhythm[1] == "Atrial fibrillation/flutter" ){next}
  }
    PatientNames[i] = sub_pat$PseudoId[1]
    HaemIndex2017[[i]]$TimeSeriesData[,1] <- HaemIndex2017[[i]]$TimeSeriesData[,1] - DP_StripTime(sub_pat$OpEnd)
    index = sample(1:max((dim(HaemIndex2017[[i]]$TimeSeriesData)[1] -1) , 1) , 1)
    if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[index , 2])){
    PreAFHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index+1 , c(5,6,8,9,1) ])
    PreAFHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index , c(2,3 , 4) ])

      }
  if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[index+1 , 2])){
    PreAFHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index+1 , c(2,3 , 4) ])
    PreAFHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index , c(5,6,8,9 , 1) ])

    }
  if(is.na(HaemIndex2017[[i]]$TimeSeriesData[2 , 2]) & is.na(HaemIndex2017[[i]]$TimeSeriesData[1 , 2])){
    PreAFHaemDataMatrix[i, ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(2,3 , 4 , 5,6,8,9,1) ])

  }
  
  if(!is.na(sub_pat$ConfirmedFirstNewAF) || !is.na(sub_pat$FirstNewAF) ){
    AFLogical[i] = 1
    
    if(!is.na(sub_pat$ConfirmedFirstNewAF)){
      if(sub_pat$ConfirmedFirstNewAF != 'CNAF'){
  
    tmp <-  HaemIndex2017[[i]]$TimeSeriesData$time + DP_StripTime(sub_pat$OpEnd)  -DP_StripTime(sub_pat$ConfirmedFirstNewAF)    
    tmp[tmp > 0] <- 1e10 
    tmp <- abs(tmp)
    index <- which.min(tmp)
      }
    }
    if(is.na(sub_pat$ConfirmedFirstNewAF) & !is.na(sub_pat$FirstNewAF)){
      tmp <-  HaemIndex2017[[i]]$TimeSeriesData$time + DP_StripTime(sub_pat$OpEnd) -DP_StripTime(sub_pat$FirstNewAF)    
      tmp[tmp > 0] <- 1e10 
      tmp <- abs(tmp)
      index <- which.min(tmp)
    }
  }
    if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[index , 2])){
      PreAFHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index + 1 , c(5,6,8,9,1) ])
      PreAFHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index , c(2,3 , 4) ])
    }
    if(!is.na(HaemIndex2017[[i]]$TimeSeriesData[index + 1 , 2])){
      PreAFHaemDataMatrix[i,c(1,2,3) ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index + 1 , c(2,3 , 4) ])
      PreAFHaemDataMatrix[i,-c(1,2,3)] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[index , c(5,6,8,9,1) ])
    }
    if(is.na(HaemIndex2017[[i]]$TimeSeriesData[2 , 2]) & is.na(HaemIndex2017[[i]]$TimeSeriesData[1 , 2])){
      PreAFHaemDataMatrix[i, ] = as.numeric(HaemIndex2017[[i]]$TimeSeriesData[2 , c(2,3 , 4 , 5,6,8,9,1) ])
      
    }
    
  if(!is.na(sub_pat$ConfirmedFirstNewAF) ){
    if(sub_pat$ConfirmedFirstNewAF[1] == "CNAF"){
      AFLogical[i] = 0
    }
  }
    HaemIndex2017[[i]]$TimeSeriesData[,1] <- HaemIndex2017[[i]]$TimeSeriesData[,1] + DP_StripTime(sub_pat$OpEnd)
  }


rownames(PreAFHaemDataMatrix) = PatientNames
colnames(PreAFHaemDataMatrix) = names(HaemIndex2017[[1]]$TimeSeriesData[1 , c(2,3 , 4, 5,6,8,9 , 1) ])

AFLogical = AFLogical[PatientNames != 0 & apply(PreAFHaemDataMatrix , 1 , function(X){sum(X, na.rm =T) }) != 0 ]
PreAFHaemDataMatrix = PreAFHaemDataMatrix[PatientNames != 0 & apply(PreAFHaemDataMatrix , 1 , function(X){sum(X, na.rm =T) }) != 0, ]
PatientNames = PatientNames[PatientNames != 0 ]
PatientNames = PatientNames[apply(PreAFHaemDataMatrix , 1 , function(X){sum(X, na.rm =T) }) != 0]

PreAFHaemDataMatrix[ , 1] = log(PreAFHaemDataMatrix[,1])
PreAFHaemDataMatrix[ , 4] = log(PreAFHaemDataMatrix[,4])
PreAFHaemDataMatrix[ , 5] = log(PreAFHaemDataMatrix[,5])
PreAFHaemDataMatrix[ , 6] = log(PreAFHaemDataMatrix[,6])

}


x11()
par(mfrow = c(1 , 7))
BC_PlotCompareSingleHists(PreAFHaemDataMatrix[AFLogical == 0 ,1],PreAFHaemDataMatrix[AFLogical == 1 ,1])
BC_PlotCompareSingleHists(PreAFHaemDataMatrix[AFLogical == 0 ,2],PreAFHaemDataMatrix[AFLogical == 1 ,2])
BC_PlotCompareSingleHists(PreAFHaemDataMatrix[AFLogical == 0 ,3],PreAFHaemDataMatrix[AFLogical == 1 ,3])
BC_PlotCompareSingleHists(PreAFHaemDataMatrix[AFLogical == 0 ,4],PreAFHaemDataMatrix[AFLogical == 1 ,4])
BC_PlotCompareSingleHists(PreAFHaemDataMatrix[AFLogical == 0 ,5],PreAFHaemDataMatrix[AFLogical == 1 ,5])
BC_PlotCompareSingleHists(PreAFHaemDataMatrix[AFLogical == 0 ,6],PreAFHaemDataMatrix[AFLogical == 1 ,6])
BC_PlotCompareSingleHists(PreAFHaemDataMatrix[AFLogical == 0 & PreAFHaemDataMatrix[,8] < 100 ,8],PreAFHaemDataMatrix[AFLogical == 1 & PreAFHaemDataMatrix[,8] < 100 ,8])


PreAFHaem <- data.frame(cbind(PreAFHaemDataMatrix , AFLogical))
PreAFHaem$AFLogical = as.factor(PreAFHaem$AFLogical)

model <- (glm(formula = AFLogical ~  WBC + Hb  + Platelets + PT   ,family=binomial(link='logit') , data=PreAFHaem ))
summary( model )
BC_PlotCompareSingleHists(model$fitted.values[model$y == 0 ],model$fitted.values[model$y == 1 ] )
