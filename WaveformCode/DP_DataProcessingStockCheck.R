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
DP_FilterPatients<-function(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess){
    # Remove early records.
    Patientnumber <- lapply(listAllPatients , DP_ReturnPatientNumber)
    Patientnumberlogical <- lapply(Patientnumber , function(X){X > 124})
    Patientnumberlogical[is.na(Patientnumberlogical)] <- FALSE
    listAllPatients <- listAllPatients[which(as.matrix(Patientnumberlogical) == TRUE)]
    # Remove files not in patient index
    InpatIndexlogical <-    lapply(listAllPatients , function(X){DP_existsinpatientindex(PatIndex2017 = PatIndex2017 , X)})
    listAllPatients <- listAllPatients[which(as.matrix(InpatIndexlogical) == TRUE)]
    # Remove files which are marked unsuable
    Usablelogical     <-    lapply(listAllPatients , function(X){DP_isusable(PatIndex2017 = PatIndex2017 , X)})
    listAllPatients <- listAllPatients[which(as.matrix(Usablelogical) == TRUE)]
    # Remove files with unusual operations
    POSRlogical <- lapply(listAllPatients, function(X){DP_FilterbySinusRhythum(PatIndex2017 = PatIndex2017 , X)})
    listAllPatients <- listAllPatients[which(as.matrix(POSRlogical) == TRUE)]
    
    POSRlogical <- lapply(listAllPatients, function(X){DP_FilterbyOps(PatIndex2017 = PatIndex2017 , X , HowtoFilterops)})
    listAllPatients <- listAllPatients[which(as.matrix(POSRlogical) == TRUE)]
    return(listAllPatients)
  }
  
  
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}
FilestoProcess = DP_ChooseECGstoProcess()
HowtoFilterops[is.na(HowtoFilterops)] <- 0

Stocklist <- matrix(0 , length(listAllPatients) ,  length(FilestoProcess) + 2)
colnames( Stocklist ) <- c(FilestoProcess , 'TotalHours' , 'v1Logical_ProcessbyCamila')
rownames( Stocklist ) <- listAllPatients


for( i in 1:length(listAllPatients)){
for(j in 1:length(FilestoProcess)){
Stocklist[i , j] <- DP_CheckFileExists(path , listAllPatients[[i]] , paste0(FilestoProcess[j], '_' , listAllPatients[[i]]) ) 
}
sub_pat = subset(PatIndex2017, PseudoId %in% listAllPatients[i])
if(nrow( sub_pat ) > 0){
Stocklist[i , length(FilestoProcess) + 1] <- sub_pat$TotalITUTimeHRS[1]
}
if(nrow( sub_pat ) == 0){
Stocklist[i , length(FilestoProcess) + 1] <- NA
next
}

# Filter out ops and pre op AF
index <- which(sub_pat$ProcDetails[1] == HowtoFilterops$Optype)


if( sub_pat$Pre_OperativeHeartRhythm[1] !=  "Sinus Rhythm"){
  Stocklist[i , length(FilestoProcess) + 1] <- NA}
if( HowtoFilterops$Removal.all.data.if.they.ever.have.this.op[index] == 1 || HowtoFilterops$Remove.for.this.op[index] == 1  ){
  Stocklist[ i , length(FilestoProcess) + 1 ] <- NA
  Stocklist[ which( rownames(Stocklist) == rownames(Stocklist)[i] ) ,  length(FilestoProcess) + 1] <- NA  
next
}
}

Stocklist <- data.frame(Stocklist[ , 1:4])
write.csv( Stocklist , file =  paste0('C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\StockList' , Sys.Date() , '.csv'))
