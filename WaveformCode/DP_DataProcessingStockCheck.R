# Script to check which files have been processed
pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )

files = list.dirs( repositorylocation <- choose.dir(caption = 'Select files to check') , full.names = FALSE , recursive = FALSE)
DataTypes = c("Discrete", "ECGI", "ECGII", "ECGIII", "CVP", "ART", "SPO2", "Flow", "Paw")

load(choose.files(caption = "Select patientindexmaster.RData"))

Stocklist <- matrix(0 , length(files) ,  length(DataTypes) + 1)

for( i in 1:length(files))
{
  for(j in 1:length(DataTypes))
{
  Stocklist[i , j] <- file.exists(paste0(repositorylocation ,'\\' 
                                          , files[i] , '\\Zip_out\\' 
                                          , DataTypes[j] ,  '_' 
                                          , files[i] , '.RData' )) 
}
sub_pat = subset(PatIndex2017, PseudoId %in% files[i])
if(nrow( sub_pat ) > 0)  
{
Stocklist[i , length(DataTypes) + 1] <- sub_pat$TotalITUTimeHRS[1]
}
if(nrow( sub_pat ) == 0)  
{
Stocklist[i , length(DataTypes) + 1] <- NA
}
}

colnames( Stocklist ) <- c(DataTypes , 'TotalHours')
rownames( Stocklist ) <- files
Stocklist <- data.frame(Stocklist)
write.csv( Stocklist , file =  paste0('C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_13022018\\StockList' , Sys.Date() , '.csv'))
