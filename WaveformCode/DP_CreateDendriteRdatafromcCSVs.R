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


numberofCSVs <- select.list( as.character(c(1:10)) , graphics = TRUE  , preselect = '2' )
DetIndex2017 <- read.csv( choose.files( multi = FALSE ) , stringsAsFactors = FALSE )


if(as.numeric(numberofCSVs) > 1){
  for(i in 2:as.numeric(numberofCSVs) )
  {
    DetIndex2017 <- rbind(DetIndex2017, read.csv(choose.files(multi = FALSE), stringsAsFactors = FALSE))
  }
}
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery == " Aortic valve"] = "Aortic Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery == "Aortic Valve "] = "Aortic Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery == "Aortic valve"] = "Aortic Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery == "Tricuspid valve"] = "Tricuspid Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery == "Tricuspid And Mitral  valve"] = "Tricuspid And Mitral Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Mitral valve "  ] =  "Mitral Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Mitral valve"  ] =  "Mitral Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Tricuspid valve "  ] =  "Tricuspid Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Tricuspid Valve "  ] =  "Tricuspid Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Tricuspid valve"  ] =  "Tricuspid Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Aortic And Mitral valve"  ] =  "Aortic And Mitral Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Aortic And Mitral And Tricuspid valve"  ] =  "Aortic And Mitral And Tricuspid Valve"
DetIndex2017$Planned.Valve.Surgery[DetIndex2017$Planned.Valve.Surgery ==  "Aortic And Mitral"  ] =  "Aortic And Mitral Valve"



save(DetIndex2017 , file = "C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\DendriteMaster.RData")
