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

numberofCSVs <- select.list( as.character(c(1:10)) , graphics = TRUE  , preselect = '2' )
HaemIndex2017 <- read.csv( choose.files( multi = FALSE ) , stringsAsFactors = FALSE )

if(as.numeric(numberofCSVs) > 1){
  for(i in 2:as.numeric(numberofCSVs) )
  {
    tmp <-  read.csv(choose.files(multi = FALSE), stringsAsFactors = FALSE)
    HaemIndex2017 <- rbind(HaemIndex2017,tmp)
    rm(tmp)
  }
}

DP_RestructureHaem <- function(HaemIndex2017){
  listoftsvariables <- names(HaemIndex2017)[19:26]
  
  uniquenames <- unique(HaemIndex2017$NewPseudoId)
  uniquenames <- uniquenames[!is.na(uniquenames)]
  
  NewData <- list()
  for( i in 1:length(uniquenames) ){
    NewData[[i]] <- setNames(list(1 , 1) , c('TimeSeriesData' , 'MetaData'))
    NewData[[i]]$TimeSeriesData <-  data.frame(time = as.POSIXct(DP_StripTime(HaemIndex2017$PostOpSampleTime[grepl(uniquenames[i] , HaemIndex2017$NewPseudoId) ]) ) , tsdata <- HaemIndex2017[grepl(uniquenames[i] , HaemIndex2017$NewPseudoId)  , 19:26] ) 
    NewData[[i]]$MetaData <- data.frame(HaemIndex2017[which(grepl(uniquenames[i] , HaemIndex2017$NewPseudoId))[1]  , -c(19:26)])
  }
  NewData <- setNames(NewData , uniquenames)
  return(NewData)
}

HaemIndex2017 <- DP_RestructureHaem( HaemIndex2017 )

save(HaemIndex2017 , file = "C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\HaemMaster.RData")
