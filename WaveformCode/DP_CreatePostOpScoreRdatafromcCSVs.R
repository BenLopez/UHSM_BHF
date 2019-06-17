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


PostOpScores <-read.csv( choose.files( multi = FALSE ) , stringsAsFactors = FALSE )

PostOpScores$Day <- apply(as.matrix(PostOpScores[ , c(2,3) ]) ,1 , function(X){gsub(paste0(X[2], '-') , '' , X[1])}) 

PostOpScores$Day  <- DP_StripTime(PostOpScores$Day )
  
uniqupatients <- unique(PostOpScores$PseudoId)

PostOpData <- list()
for(ii in 1:length(uniqupatients)){

  PostOpData[[ii]] <- data.frame(Day = PostOpScores$Day[ PostOpScores$PseudoId == uniqupatients[ii] ] ,
                                 logCASUS = PostOpScores$logCASUS[ PostOpScores$PseudoId == uniqupatients[ii] ] ,
                                 RACE = PostOpScores$RACE[ PostOpScores$PseudoId == uniqupatients[ii] ] ,
                                 SOFA = PostOpScores$SOFA[ PostOpScores$PseudoId == uniqupatients[ii] ])

  PostOpData[[ii]] <- PostOpData[[ii]][order(PostOpData[[ii]]$Day) , ]
}

PostOpData <- setNames(PostOpData , uniqupatients)
PostOpScoresIndex2017 <- PostOpData
  
save(PostOpScoresIndex2017 , file = "C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\PostOpScoresMaster.RData")

  