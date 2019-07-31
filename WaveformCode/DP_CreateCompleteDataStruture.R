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

{NamesTimeSeriesVariables  <- c( names( BioChemIndex2017[[1]]$TimeSeriesData ) , 
                                 names( FlowIndex2017[[1]]$TimeSeriesData ) ,
                                 names( FluidsIndex2017[[1]]$TimeSeriesData ) ,
                                 names( HaemIndex2017[[1]]$TimeSeriesData ) ,
                                 names( VentIndex2017[[1]]$TimeSeriesData ))

NamesTimeSeriesVariables <- sort(unique(NamesTimeSeriesVariables))
NumberofTimeSeriesVariables <- length(NamesTimeSeriesVariables)

NamesDiscreteVariables  <- c( names( DetIndex2017 ),
                              names( PatIndex2017 ),
                              names(VentEpisodesIndex2017),
                              names( BioChemIndex2017[[1]]$MetaData ), 
                              names( FlowIndex2017[[1]]$MetaData ),
                              names( FluidsIndex2017[[1]]$MetaData ),
                              names( HaemIndex2017[[1]]$MetaData ),
                              names( VentIndex2017[[1]]$MetaData ) )
NamesDiscreteVariables  <- sort(unique(NamesDiscreteVariables))
NumberDiscreteVariables <- length(NamesDiscreteVariables)
}

{CCD_CreateDefaultMatrix <- function(NamesX , nX , nc = 0){
  X = matrix(NA , nc , nX)
  colnames(X) <- NamesX
  return(data.frame(X))
}
CCD_PopulateDiscreteMatrix <- function(NewPuesdoId,
                                       NamesDiscreteVariables, 
                                       NumberDiscreteVariables,
                                       DetIndex2017,
                                       PatIndex2017,
                                       VentEpisodesIndex2017,
                                       BioChemIndex2017,
                                       FlowIndex2017,
                                       HaemIndex2017,
                                       VentIndex2017){
  output <- CCD_CreateDefaultMatrix(NamesDiscreteVariables , NumberDiscreteVariables , 1)
  
  if( sum(PatIndex2017$NewPseudoId == NewPuesdoId) > 0 ){
    output[1 , NamesDiscreteVariables %in% names(PatIndex2017)] <- PatIndex2017[PatIndex2017$NewPseudoId == NewPuesdoId, sort(names(PatIndex2017))]
  }
  if( sum(DetIndex2017$NewPseudoId == NewPuesdoId) > 0 ){
    output[1 , NamesDiscreteVariables %in% names(DetIndex2017)] <- DetIndex2017[DetIndex2017$NewPseudoId == NewPuesdoId, sort(names(DetIndex2017))]
  }
  if( sum(VentEpisodesIndex2017$RelevantAdmission == NewPuesdoId) > 0 ){
    output[1 , NamesDiscreteVariables %in% names(VentEpisodesIndex2017)] <- VentEpisodesIndex2017[VentEpisodesIndex2017$RelevantAdmission == NewPuesdoId, sort(names(VentEpisodesIndex2017))]
  }
  if( sum(names(BioChemIndex2017) == NewPuesdoId) > 0 ){
    output[1 , NamesDiscreteVariables %in% names(BioChemIndex2017[[which(names(BioChemIndex2017) == NewPuesdoId)]]$MetaData)] <- BioChemIndex2017[[which(names(BioChemIndex2017) == NewPuesdoId)]]$MetaData[ , sort(names(BioChemIndex2017[[which(names(BioChemIndex2017) == NewPuesdoId)]]$MetaData))]
  }
  if( sum(names(HaemIndex2017) == NewPuesdoId) > 0 ){
    output[1 , NamesDiscreteVariables %in% names(HaemIndex2017[[which(names(HaemIndex2017) == NewPuesdoId)]]$MetaData)] <- HaemIndex2017[[which(names(HaemIndex2017) == NewPuesdoId)]]$MetaData[ , sort(names(HaemIndex2017[[which(names(HaemIndex2017) == NewPuesdoId)]]$MetaData))]
  }
  
  
  return(output)
  
}
CCD_ExtractTimeSeriesVariables <- function(TimeSeriesStruct ,NamesTimeSeriesVariables ,NumberofTimeSeriesVariables ){
  
  output <- CCD_CreateDefaultMatrix(NamesTimeSeriesVariables , NumberofTimeSeriesVariables , dim(TimeSeriesStruct)[1])
  output[ , NamesTimeSeriesVariables %in% names(TimeSeriesStruct) ] <- TimeSeriesStruct[ , sort(names(TimeSeriesStruct)[names(TimeSeriesStruct)%in%NamesTimeSeriesVariables]) ] 
  return(output)
  
}
CCD_PopulateTimeSeriesMatrix <- function(NewPuesdoId,
                                         NamesTimeSeriesVariables, 
                                         NumberofTimeSeriesVariables,
                                         DetIndex2017,
                                         PatIndex2017,
                                         VentEpisodesIndex2017,
                                         BioChemIndex2017,
                                         FlowIndex2017,
                                         HaemIndex2017,
                                         VentIndex2017,
                                         FluidsIndex2017){
  
  output <- CCD_CreateDefaultMatrix(NamesTimeSeriesVariables , NumberofTimeSeriesVariables , 0)

  if(sum(which(names(BioChemIndex2017) == NewPuesdoId))>0){
    output <- rbind(output , CCD_ExtractTimeSeriesVariables(BioChemIndex2017[[which(names(BioChemIndex2017) == NewPuesdoId)]]$TimeSeriesData ,NamesTimeSeriesVariables, NumberofTimeSeriesVariables))  
  }
  if(sum(which(names(FlowIndex2017) == NewPuesdoId))>0){
    output <- rbind(output , CCD_ExtractTimeSeriesVariables(FlowIndex2017[[which(names(FlowIndex2017) == NewPuesdoId)]]$TimeSeriesData ,NamesTimeSeriesVariables, NumberofTimeSeriesVariables))  
  }
  if(sum(which(names(HaemIndex2017) == NewPuesdoId))>0){
    output <- rbind(output , CCD_ExtractTimeSeriesVariables(HaemIndex2017[[which(names(HaemIndex2017) == NewPuesdoId)]]$TimeSeriesData ,NamesTimeSeriesVariables, NumberofTimeSeriesVariables))  
  }
  if(sum(which(names(FluidsIndex2017) == NewPuesdoId))>0){
    output <- rbind(output , CCD_ExtractTimeSeriesVariables(FluidsIndex2017[[which(names(FluidsIndex2017) == NewPuesdoId)]]$TimeSeriesData ,NamesTimeSeriesVariables, NumberofTimeSeriesVariables))  
  }
  output
  return(output[order(output$time),])
}
}


AllDataStructure <- list()
for(ii in 1:length(PatIndex2017$NewPseudoId)){
  
  tmpdisdata <-   CCD_PopulateDiscreteMatrix( PatIndex2017$NewPseudoId[ii],
                                              NamesDiscreteVariables, 
                                              NumberDiscreteVariables,
                                              DetIndex2017,
                                              PatIndex2017,
                                              VentEpisodesIndex2017,
                                              BioChemIndex2017,
                                              FlowIndex2017,
                                              HaemIndex2017,
                                              VentIndex2017)
  tmptsdata <- CCD_PopulateTimeSeriesMatrix(PatIndex2017$NewPseudoId[ii],
                                           NamesTimeSeriesVariables, 
                                           NumberofTimeSeriesVariables,
                                           DetIndex2017,
                                           PatIndex2017,
                                           VentEpisodesIndex2017,
                                           BioChemIndex2017,
                                           FlowIndex2017,
                                           HaemIndex2017,
                                           VentIndex2017,
                                           FluidsIndex2017)
  
  AllDataStructure[[ii]] <- list( DiscreteData = tmpdisdata, timeseriesvariables = tmptsdata )
  DP_WaitBar(ii/length(PatIndex2017$NewPseudoId))
}
rm(tmpdisdata , tmptsdata)
names(AllDataStructure) <- PatIndex2017$NewPseudoId

save(AllDataStructure , file = "C:\\Users\\Ben\\Desktop\\UHSM_Cardiac_06082018\\AllDataStructure.RData")

