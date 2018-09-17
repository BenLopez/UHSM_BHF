# Script to check which files have been processed
pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )

files = list.dirs( repositorylocation <- choose.dir(caption = 'Select files to check') , full.names = FALSE , recursive = FALSE)
DataTypes = c("Discrete", "ECGI", "ECGII", "ECGIII", "CVP", "ART", "SPO2", "Flow", "Paw")

DP_LoadPatientIndex()
HowtoFilterops <- read.csv(choose.files(caption = "Select listofopsSH") , stringsAsFactors = FALSE)
Processv2 <- read.csv(choose.files(caption = 'Select Downloadedfiles') , stringsAsFactors = FALSE)
HowtoFilterops[is.na(HowtoFilterops)] = 0

Stocklist <- matrix(0 , length(files) ,  length(DataTypes) + 2)
colnames( Stocklist ) <- c(DataTypes , 'TotalHours' , 'v1Logical_ProcessbyCamila')
rownames( Stocklist ) <- files


for( i in 1:length(files)){
for(j in 1:length(DataTypes)){
Stocklist[i , j] <- file.exists(paste0(repositorylocation ,'\\' 
                                          , files[i] , '\\Zip_out\\' 
                                          , DataTypes[j] ,  '_' 
                                          , files[i] , '.RData' )) 
}
sub_pat = subset(PatIndex2017, PseudoId %in% files[i])
if(nrow( sub_pat ) > 0){
Stocklist[i , length(DataTypes) + 1] <- sub_pat$TotalITUTimeHRS[1]
}
if(nrow( sub_pat ) == 0){
Stocklist[i , length(DataTypes) + 1] <- NA
next
}

# Filter out ops and pre op AF
index <- which(sub_pat$ProcDetails[1] == HowtoFilterops$Optype)
#if(length(index) == 0){
#Stocklist[i , length(DataTypes) + 1] <- NA
#next
#}

if( sum(Processv2$X == files[[i]])  ==1  ){
  Stocklist[ i , 11 ] <- 1
  Stocklist[ i , 10 ] <- sub_pat$TotalITUTimeHRS[1]
  next
}

if( sub_pat$Pre_OperativeHeartRhythm[1] !=  "Sinus Rhythm"){
  Stocklist[i , length(DataTypes) + 1] <- NA}
if( HowtoFilterops$Removal.all.data.if.they.ever.have.this.op[index] == 1 || HowtoFilterops$Remove.for.this.op[index] == 1  ){
  Stocklist[ i , length(DataTypes) + 1 ] <- NA
  Stocklist[ which( rownames(Stocklist) == rownames(Stocklist)[i] ) ,  length(DataTypes) + 1] <- NA  
next
}
}

Stocklist <- data.frame(Stocklist)
write.csv( Stocklist , file =  paste0('C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_02072018\\StockList' , Sys.Date() , '.csv'))
