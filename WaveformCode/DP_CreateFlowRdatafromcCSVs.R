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
FlowIndex2017 <- read.csv( choose.files( multi = FALSE ) , stringsAsFactors = FALSE )

if(as.numeric(numberofCSVs) > 1){
  for(i in 2:as.numeric(numberofCSVs) )
  {
    FlowIndex2017 <- rbind(FlowIndex2017, read.csv(choose.files(multi = FALSE), stringsAsFactors = FALSE))
  }
}

FlowIndex2017 <- DP_RestructureFlow( FlowIndex2017 )
save(FlowIndex2017 , file = "C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\FlowMaster.RData")
