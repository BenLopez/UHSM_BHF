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
BioChemIndex2017 <- read.csv( choose.files( multi = FALSE ) , stringsAsFactors = FALSE )

if(as.numeric(numberofCSVs) > 1){
  for(i in 2:as.numeric(numberofCSVs) )
  {
    BioChemIndex2017 <- rbind(BioChemIndex2017, read.csv(choose.files(multi = FALSE), stringsAsFactors = FALSE))
  }
}

BioChemIndex2017 <- DP_RestructureBioChem( BioChemIndex2017 )
save(BioChemIndex2017 , file = "C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\BIOChemMaster.RData")
