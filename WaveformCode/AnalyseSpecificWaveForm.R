pathFiles <- choose.dir(caption="Select folder with source code.")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

source("LibrariesAndSettings.R" , print.eval  = TRUE )

load(choose.files(caption = "Select patientindexmaster.RData"))

setofrepositorynumbers <- c('1' , '2' , '3' , '4' , '5' , '6' , '7' , '8')
numberrep =  as.numeric(select.list(setofrepositorynumbers, preselect = setofrepositorynumbers[1], multiple = FALSE,
                       title = 'Select number of reporsitory locations', graphics = TRUE ))

# Select reporsitories
path = list()
for(i in 1:numberrep)
{
path[[i]] = choose.dir( caption  <- paste0( "Select folder " ,i,  " containing data repository" ))
if(i == 1){  listAllPatients <- as.matrix(list.dirs(path = path[[i]], full.names = FALSE, recursive = FALSE))  }
if(i > 1){   listAllPatients <- rbind( listAllPatients , as.matrix(list.dirs(path = path[[i]], full.names = FALSE, recursive = FALSE)) ) }
}
  
source( 'ASWF_ChooseLoadandProcessPatient.R' )

source( 'ASWF_InteractivePatientAnalysis.R' )
