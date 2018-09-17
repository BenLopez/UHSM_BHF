{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

# Optionsforprocessing

DP_LoadPatientIndex()
HoursBeforeAndAFter <- DP_SelectHoursBeforeandAfter()
FilestoProcess <- DP_ChooseECGstoProcess() 
HowtoFilterops <- read.csv(choose.files(caption = "Select listofopsSH") , stringsAsFactors = FALSE)

DP_ChooseDataReps()
listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
sublist <- select.list(listAllPatients , graphics = TRUE , multiple = FALSE)

ECGs <- DP_LoadReducedECGs( path , sublist , numberrep  , FilestoProcess)

outputdata <- DP_LoadRpeaksfile( path , sublist )

test <- CleanRpeaks( PE_MultipleECGRPeaks( outputdata , ECGs ,  thresh = 0.02 ) )
plot(outputdata$RRCombined$t , outputdata$RRCombined$RR , col = rgb(0 , 1 , 0 , alpha = 0.1) )
points(test$t , test$RR, col = rgb(1 , 0 , 0 , alpha = 0.1) )
