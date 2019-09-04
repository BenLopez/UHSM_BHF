{
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
    set.seed( 1 )
  }
}

VectorofDifferences <- matrix(NA , length(listAllPatients) , 1)
AFLogical <- matrix(NA , length(listAllPatients))

counter <- 1
for(PatientID in listAllPatients[1:length(listAllPatients)] ){
AFLogical[counter] <- DP_CheckIfAFPatient(DP_ExtractPatientRecordforIndex(PatIndex2017 ,PatientID ))
if( file.exists(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData')) ){ 
  
load(paste0(path ,'\\',PatientID,'\\Zip_out\\', "HMOutput" , PatientID , '.RData'))
t <- FM_ExtractTimeFromHMOutput( outputstruct )

AFAnnotationHM <- FM_CreateAFAnnoation( outputstruct[[length(outputstruct)]] )[1:(length(outputstruct)-2)]
AFAnnotationExpert <- BC_CreateAFAnnotationFomMetaData(t , DP_ExtractPatientRecordforIndex(PatIndex2017 = PatIndex2017,PatientID))

if(sum(AFAnnotationExpert)>0){
VectorofDifferences[counter] <- difftime(t[which(AFAnnotationHM > 0)[1]] ,t[which(AFAnnotationExpert > 0)[1]] , units = c('hours'))
}else{
VectorofDifferences[counter] <- 0   
}
if(sum(AFAnnotationExpert)==0 & sum(AFAnnotationHM)>0){
  VectorofDifferences[counter] <- -Inf   
}
if(sum(AFAnnotationExpert)>0 & sum(AFAnnotationHM)==0){
  VectorofDifferences[counter] <-  Inf   
}
}
DP_WaitBar(counter/length(listAllPatients))
counter <- counter + 1
}

x11()
p_results <- ggplot(data.frame(x = c(1:sum(AFLogical)) ,y=c(VectorofDifferences[AFLogical == T])),aes(x,y)) + 
  geom_point(col = 'blue') +
  ggtitle('Comparison of Algorithm and Expert Annotation') +xlab('Patinet Index') + ylab('Difference in Annotation Time (hours)')
p_results