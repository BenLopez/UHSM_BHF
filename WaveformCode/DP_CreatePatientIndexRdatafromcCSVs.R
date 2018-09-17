{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}


numberofCSVs <- select.list(as.character(c(1:10)) , graphics = TRUE  , preselect = '2' )
PatIndex2017 <- read.csv(choose.files(multi = FALSE), stringsAsFactors = FALSE)


if(as.numeric(numberofCSVs) > 1)
{
for(i in 2:as.numeric(numberofCSVs) )
{
PatIndex2017 <- rbind(PatIndex2017, read.csv(choose.files(multi = FALSE), stringsAsFactors = FALSE))
}
}

save(PatIndex2017 , file = "C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\PatientIndexMaster.RData")
