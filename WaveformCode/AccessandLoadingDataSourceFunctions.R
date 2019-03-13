# Script with WaveForm data processing and access functions.

###### Loading functions ######


DP_LoadPatientIndex <- function(){
  filetype = select.list(c('csv' , 'RData'), preselect = NULL, multiple = TRUE,
                         title = 'Choose File Type For Patient Index', graphics = TRUE )
if(filetype == 'csv')
{
     PatIndex2017 <<- read.csv(file=choose.files(caption="Select 2017 PatientIndex.csv file"), stringsAsFactors = FALSE)
}    
if(filetype == 'RData')
{
  load(choose.files( caption = 'Select PatientIndexMaster.RData'))
  PatIndex2017 <<- PatIndex2017
}
  return(PatIndex2017)
}
DP_LoadPatientDendrite <- function(){
  filetype = select.list(c('csv' , 'RData'), preselect = NULL, multiple = TRUE,
                         title = 'Choose File Type For Patient Index', graphics = TRUE )
  if(filetype == 'csv')
  {
    DetIndex2017 <<- read.csv(file=choose.files(caption="Select 2017 DendriteMaster.csv file"), stringsAsFactors = FALSE)
  }    
  if(filetype == 'RData')
  {
    load(choose.files( caption = 'Select DendriteMaster.RData'))
    DetIndex2017 <<- DetIndex2017
  }
  return(DetIndex2017)
}

DP_LoadPatientBioChem <- function(){
  filetype = select.list(c('csv' , 'RData'), preselect = NULL, multiple = TRUE,
                         title = 'Choose File Type For Patient Index', graphics = TRUE )
  if(filetype == 'csv')
  {
    BioChemIndex2017 <<- read.csv(file=choose.files(caption="Select 2017 BioChemistry.csv file"), stringsAsFactors = FALSE)
  }    
  if(filetype == 'RData')
  {
    load(choose.files( caption = 'Select BIOChemMaster.RData'))
    BioChemIndex2017 <<- BioChemIndex2017
  }
  return(BioChemIndex2017)
}

DP_LoadPatientVent <- function(){
  filetype = select.list(c('csv' , 'RData'), preselect = NULL, multiple = TRUE,
                         title = 'Choose File Type For Patient Index', graphics = TRUE )
  if(filetype == 'csv')
  {
    VentIndex2017 <<- read.csv(file=choose.files(caption="Select 2017 Ventilator.csv file"), stringsAsFactors = FALSE)
  }    
  if(filetype == 'RData'){
    load(choose.files( caption = 'Select VentMaster.RData'))
    VentIndex2017 <<- VentIndex2017
  }
  return(VentIndex2017)
}

DP_LoadECG <- function(path , subList , numberrep=1 , ECGNum = 1 ){
  # Function to load ECG data. ECGNum = (1 , 2 , 3) or ('I' , 'II' , 'III')
  if(ECGNum  == 1 || ECGNum  == 'I' || ECGNum  == 'ECGI')
  {
    # Load wave form data 
    for(i in 1:(numberrep+1))
    {
      if(i > (numberrep))
      {
        stop('Error: No ECGI data processed.') 
        break
      }
      if(file.exists(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGI_' , subList , '.RData' )))
      {
        load(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGI_' , subList , '.RData' ))
        break 
      }  
    }
  }
  
  if(ECGNum  == 2 || ECGNum  == 'II'|| ECGNum  == 'ECGII')
  {
    # Load wave form data 
    for(i in 1:(numberrep+1))
    {
      if(i > (numberrep))
      {
        stop('Error: No ECGII data processed.') 
        break
      }
      if(file.exists(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGII_' , subList , '.RData' )))
      {
        load(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGII_' , subList , '.RData' ))
        break 
      }  
    }
  }
  
  if(ECGNum  == 3 || ECGNum  == 'III'|| ECGNum  == 'ECGIII')
  {
    # Load wave form data 
    for(i in 1:(numberrep+1))
    {
      if(i > (numberrep))
      {
        stop('Error: No ECGIII data processed.') 
        break
      }
      if(file.exists(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGIII_' , subList , '.RData' )))
      {
        load(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'ECGIII_' , subList , '.RData' ))
        break 
      }  
    }
  }
  
  return(WaveData)
}
DP_LoadECGReduced <- function(path , subList , numberrep=1 , ECGNum = 1 ){
  # Function to load ECG data. ECGNum = (1 , 2 , 3) or ('I' , 'II' , 'III')
  if(ECGNum  == 1 || ECGNum  == 'I' || ECGNum  == 'ECGI')
  {
    # Load wave form data 
    for(i in 1:(numberrep+1))
    {
      if(i > (numberrep))
      {
        stop('Error: No ECGI data processed.') 
        break
      }
      if(file.exists(paste0(path , '\\' , subList , '\\Zip_out\\' , subList , '_', "ECGI_reduced" , '.RData') ))
      {
        load(paste0(path , '\\' , subList , '\\Zip_out\\' , subList , '_', "ECGI_reduced" , '.RData'))
        break 
      }  
    }
  }
  
  if(ECGNum  == 2 || ECGNum  == 'II'|| ECGNum  == 'ECGII')
  {
    # Load wave form data 
    for(i in 1:(numberrep+1))
    {
      if(i > (numberrep))
      {
        stop('Error: No ECGII data processed.') 
        break
      }
      if(file.exists(paste0(path , '\\' , subList , '\\Zip_out\\' , subList , '_', "ECGII_reduced" , '.RData') ))
      {
        load(paste0(path , '\\' , subList , '\\Zip_out\\' , subList , '_', "ECGII_reduced" , '.RData'))
        break 
      }  
    }
  }
  
  if(ECGNum  == 3 || ECGNum  == 'III'|| ECGNum  == 'ECGIII')
  {
    # Load wave form data 
    for(i in 1:(numberrep+1))
    {
      if(i > (numberrep))
      {
        stop('Error: No ECGIII data processed.') 
        break
      }
      if(file.exists(paste0(path , '\\' , subList , '\\Zip_out\\' , subList , '_', "ECGIII_reduced" , '.RData')))
      {
        load(paste0(path , '\\' , subList , '\\Zip_out\\' , subList , '_', "ECGIII_reduced" , '.RData'))
        break 
      }  
    }
  }
  
  return(WaveData)  
}
DP_LoadRpeaksfileECGI <- function(path , PatientsId ){
  load(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , PatientsId , '_RPeaks' , '.RData'))
  return(outputdata$ECGI)
}  
DP_LoadRpeaksfile <- function(path , PatientsId ){
  load(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , PatientsId , '_RPeaks' , '.RData'))
  return(outputdata)
}  
DP_LoadECGs<- function(path , subList , numberrep=1 , FilestoProcess){
  output = list()
  for(i in 1:length(FilestoProcess)){
    output[[i]] <- DP_LoadECG( path , subList , numberrep , FilestoProcess[i] )
  }
  output <- setNames(output , FilestoProcess)
}
DP_LoadReducedECGs <- function(path , subList , numberrep=1 , FilestoProcess){
  output = list()
  for(i in 1:length(FilestoProcess)){
    output[[i]] <- DP_LoadECGReduced( path , subList , numberrep , FilestoProcess[i] )
  }
  output <- setNames(output , FilestoProcess)
}
DP_LoadDistributionSummaries <- function(path , PatientsID){
  Name <- paste0(PatientsID , '_DistributionSummaries' )
  return(DP_LoadFile(path , PatientsID , Name = Name ))
}
DP_LoadFile <- function(path , PatientsID , Name){
  return(DP_loadRData(paste0(path , '\\' , PatientsID , '\\Zip_out\\', Name , '.RData')))
}
DP_loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
DP_LoadPatientsandProcessRPeaks <- function(path , subList , numberrep =1 , LoadReduced = 1,PatientRecord = DP_CreateDummyMetaData(PatIndex2017  , Name = subList)){
  # Load waveforms
  FilestoProcess <- c('ECGI' ,'ECGII', 'ECGIII' )
  if(LoadReduced == 0){
    ECGs <- DP_LoadECGs(path = path , subList = subList , numberrep =numberrep , FilestoProcess = FilestoProcess)
  }else{
    ECGs <- DP_LoadReducedECGs(path = path , subList = subList , numberrep =numberrep , FilestoProcess = FilestoProcess)
  }
  outputdata <- DP_ProcessRpeaksMultipleECGs(ECGs = ECGs , PatientRecord = PatientRecord)
  # Process Rpeaks 
  return(outputdata)
}  


##### Choosing Functions #####

DP_ChooseDataReps <- function(){
  
setofrepositorynumbers = c('1' , '2' , '3' , '4' , '5' , '6' , '7' , '8')
numberrep <<-  as.numeric(select.list(setofrepositorynumbers, preselect = setofrepositorynumbers[1], multiple = FALSE,
                                    title = 'Select number of reporsitory locations', graphics = TRUE ))

# Select reporsitories
path <<- list()
for(i in 1:numberrep)
{
  path[[i]]  <<- choose.dir( caption  = paste0( "Select folder " ,i,  " containing data repository" ))
  if(i == 1){  listAllPatients <<- as.matrix(list.dirs(path = path[[i]], full.names = FALSE, recursive = FALSE))  }
  if(i > 1){   listAllPatients <<- rbind( listAllPatients , as.matrix(list.dirs(path = path[[i]], full.names = FALSE, recursive = FALSE)) ) }
}


}
DP_choosepatient <- function(listAllPatients){
  subList <<- select.list(listAllPatients
                        , preselect = NULL
                        , multiple = FALSE
                        , title = 'Select Patient to Analyse'
                        , graphics = TRUE )
  return(subList)
}
DP_choosepatients <- function(listAllPatients){
  subList <<- select.list(listAllPatients
                          , preselect = NULL
                          , multiple = TRUE
                          , title = 'Select Patient to Analyse'
                          , graphics = TRUE )
}
DP_SelectInterestingTimePoint <- function(Wavedata , sub_pat){
  subdata <-  unique(as.vector(as.character(round.POSIXt(Wavedata$Date , units = c('hours')))))   
if(!is.na(as.POSIXct(sub_pat$FirstNewAF[1])))
  {
    if(!is.na(sub_pat$ConfirmedFirstNewAF[1]) )
    {
      interestingtimepoint <-  select.list(subdata
                                           , preselect = subdata[ which.min( abs(difftime( as.POSIXct(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1]) , tz ='GMT') , as.POSIXct(subdata, tz ='GMT') , units = 'hours')))]
                                           , multiple = FALSE
                                           , title = 'Select Interest Time Point'
                                           , graphics = TRUE )
    }
    
    if(is.na(sub_pat$ConfirmedFirstNewAF[1]) )
    {
      interestingtimepoint <-  select.list(subdata
                                           , preselect = subdata[ which.min( abs(difftime( as.POSIXct(DP_StripTime(sub_pat$FirstNewAF[1]) , tz ='GMT') , as.POSIXct(subdata, tz ='GMT') , units = 'hours')))]
                                           , multiple = FALSE
                                           , title = 'Select Interest Time Point'
                                           , graphics = TRUE )
    }
    
  }
  
  if(is.na(as.POSIXct(sub_pat$FirstNewAF[1])))
  {
    interestingtimepoint <- select.list(subdata  
                                                           , preselect = NULL
                                                           , multiple = FALSE
                                                           , title = 'Select Interest Time Point'
                                                           , graphics = TRUE )
  }
  

  
return(interestingtimepoint[1])
}
DP_SelectHoursBeforeandAfter <- function(){ 
  numberhours <- c('1' , '2' , '3' , '4' , '5' , '6' , '7' , '8')
  numberhoursbefore <- as.numeric(select.list(numberhours
                                             , preselect = '5'
                                             , multiple = FALSE
                                             ,title = 'Select number of hours before.'
                                             , graphics = TRUE ))
  numberhoursafter <- as.numeric(select.list(numberhours
                                            , preselect = '1'
                                            , multiple = FALSE
                                            ,title = 'Select number of hours after.'
                                            , graphics = TRUE ))
return(setNames(list(numberhoursbefore , numberhoursafter) , c('numberhoursbefore' , 'numberhoursafter')))
}
DP_ChooseHoursBeforeandAfter <- function(){ 
  return(DP_SelectHoursBeforeandAfter())
}
DP_CropWaveData <- function(WaveData , timeindex , HoursBeforeAndAFter){
if(length(timeindex) == 2){
  return( WaveData[ ((WaveData$Date >=  timeindex[1])*(WaveData$Date <=  timeindex[2]) == 1 ) , ] )
}
  
if(length(timeindex) == 1){  
  numberhoursbefore <- HoursBeforeAndAFter[['numberhoursbefore']]
  numberhoursafter <- HoursBeforeAndAFter[['numberhoursafter']]
  indent <- as.numeric(abs(WaveData[1 , 1] - WaveData[2 , 1]))
  if(indent > 0.012){warning('Time indent greater than 0.012. Data measured at unusual resolution.')}
  WaveData <- WaveData[ max( 1 , timeindex - (numberhoursbefore*((60^2)/indent)) ) : min(length(WaveData[ , 1]) , timeindex + (numberhoursafter*((60^2)/indent)) ) , ]
  WaveData <- ReturnWaveformwithPositiveOrientation(WaveData)
  return(WaveData)
}
}
DP_ChooseECGstoProcess <- function(){
  DataTypes = c("ECGI", "ECGII", "ECGIII")
  chooseWave2Read = select.list(DataTypes, preselect = DataTypes,
                                multiple = TRUE, graphics = TRUE, title = "Choose Waves to Read")
  return(chooseWave2Read)
}
DP_ChooseWaveformstoProcess <- function( ){
  DataTypes = c("Discrete", "ECGI", "ECGII", "ECGIII", "CVP", "ART", "SPO2", "Flow", "Paw")
  chooseWave2Read = select.list(DataTypes, preselect = DataTypes,
                                multiple = TRUE, graphics = TRUE, title = "Choose Waves to Read")
  return(chooseWave2Read)
}
DP_SelectPrecomputedFolder <- function(){
  precomputedfolderpath <- choose.dir(caption = 'Select directory where you would like the precomputed folder to be placed.')
  if(substr(precomputedfolderpath , start = nchar(precomputedfolderpath) - 38 , stop = nchar(precomputedfolderpath)) != 'PrecomputedOutputsForPopulationAnalysis'){
    precomputedfolderpath <- paste0(precomputedfolderpath , '\\PrecomputedOutputsForPopulationAnalysis')  
  }
  if(dir.exists(precomputedfolderpath) == FALSE){
    dir.create(precomputedfolderpath)
  }
  return(precomputedfolderpath)
}
DP_GetDirectories <- function( ){
  DP_LoadPatientIndex()
  DP_ChooseDataReps()
}  
DP_ChooseRegionofInterest <- function(ECG ){
  timelist <- as.vector(as.character(round.POSIXt(ECG[seq(from = 1 , to = length(ECGI[ , 1]) , by = 1000), 1] , units = 'mins')))
  
  startindex <- which.min( abs( as.POSIXct( ECGI$Date ) - as.POSIXct(select.list(unique(timelist)
                                                                                 , preselect = NULL
                                                                                 , multiple = FALSE
                                                                                 , title = 'Select time to view ECGI'
                                                                                 , graphics = TRUE ) ) ) )[1]
  
  timeinterval <- as.numeric(select.list(unique(as.character(c(1:100)))
                                         , preselect = '10'
                                         , multiple = FALSE
                                         , title = 'Select number of seconds of full ECG to be viewed'
                                         , graphics = TRUE ))
  
  endindex <- startindex + (round(timeinterval / as.numeric(abs(ECG$Date[1]- ECG$Date[2]))))
  regionofinterest <- startindex:endindex
  return(regionofinterest)
}
DP_SelectTimetoview <- function(t){
  timelist <- as.vector(as.character(round.POSIXt(t[seq(from = 1 , to = length(t) , by = 1000)] , units = 'mins')))
  startindex = which.min( abs( as.POSIXct(t) - as.POSIXct(select.list(unique(timelist)
                                                                      , preselect = NULL
                                                                      , multiple = FALSE
                                                                      , title = 'Select time to view'
                                                                      , graphics = TRUE ) ) ) )[1]
  return(t[startindex])
}

##### Checking functions ######
DP_checkfilesprocessed <- function(path , PatientsId , FilestoProcess){
  output <- rep(0 , length(FilestoProcess) , 1)
  for( i in 1:length(FilestoProcess) )
  {
    output[i] <- file.exists(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , FilestoProcess[i] , '_' , PatientsId , '.RData')  )
  }
  return(output)
}
DP_existsinpatientindex <- function(PatIndex2017 , PatientsId){
  return(nrow(subset( PatIndex2017, PseudoId %in% PatientsId )) > 0)  
}
DP_isusable <- function(PatIndex2017 , PatientsId){
  sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
  return(sub_pat$Usable[1] == 1)    
}
DP_CheckFileExists <- function(path , PatientsID , Name){
  file.exists(paste0(path , '\\' , PatientsID , '\\Zip_out\\', Name , '.RData'))
}
DP_CheckDistributionSummariesExists <- function(path , PatientsID){
  DP_CheckFileExists(path , PatientsID , Name = paste0(PatientsID , '_DistributionSummaries' ))
}
DP_CheckFieldExists <- function(DF , Field){
  return(any(names(DF) == Field))
}
DP_CheckIfAFPatient <- function(MetaData){
  if(!is.na(MetaData$ConfirmedFirstNewAF[1]) & (MetaData$ConfirmedFirstNewAF[1] != 'CNAF')){
    return(TRUE)
  }else{
    return(FALSE)
  }
} 
DP_checkRpeaksfilesprocessed<- function(path , PatientsId ){
  
  output <- file.exists(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , PatientsId , '_RPeaks' , '.RData')  )
  return(output)
}
DP_CheckECGreducedfilesprocessed <- function( path , PatientsId  , Filestoprocess){
  output <- file.exists(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , PatientsId , '_', Filestoprocess , '.RData')  )
  return(output)
}
DP_CheckfileinPrecomputedfolder <- function(precomputedfolderpath , file){
  return(file.exists(paste0(precomputedfolderpath ,'\\' ,  file)))  
}


##### Filtering functions #####
DP_FilterbySinusRhythum <- function(PatIndex2017 , PatientsId){
  sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
  return(sub_pat$Pre_OperativeHeartRhythm[1] == 'Sinus Rhythm')
}
DP_FilterbyOps <- function(PatIndex2017 , PatientsId , HowtoFilterops){
  sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
  
  index <- lapply(sub_pat$ProcDetails , function(X){which(X == HowtoFilterops$Optype)} )
  keep <- as.matrix(lapply(index , function(X){DP_Filtbyopsindividualrecord(X , HowtoFilterops = HowtoFilterops)}))
  return( sum(keep == 0) == 0)
}
DP_Filtbyopsindividualrecord <- function(index, HowtoFilterops)  {
  if(!is.na(HowtoFilterops$Keep[index]) ){keep = 1}
  if(!is.na(HowtoFilterops$Remove.for.this.op[index])){keep = 1}
  if(!is.na(HowtoFilterops$Removal.all.data.if.they.ever.have.this.op[index])){keep = 0}  
  return(keep)
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
  # Remove files without appropriate number of hours
  Numberofhourslogical <- lapply(listAllPatients , function(X){DP_numberhoursgreaterthan(PatIndex2017 = PatIndex2017 , X , numberofhours = 6 ) })
  listAllPatients <- listAllPatients[which(as.matrix(Numberofhourslogical) == TRUE)]
  # Remove files with unusual operations
  POSRlogical <- lapply(listAllPatients, function(X){DP_FilterbySinusRhythum(PatIndex2017 = PatIndex2017 , X)})
  listAllPatients <- listAllPatients[which(as.matrix(POSRlogical) == TRUE)]
  
  # Remove files without processed data
  ProcessedFiles    <-    lapply(listAllPatients , function(X){DP_checkfilesprocessed(  path = path , PatientsId = X ,  FilestoProcess = FilestoProcess )})
  Processedlogical  <-    lapply(ProcessedFiles ,  function(X){ sum(X)==length(FilestoProcess) })
  listAllPatients <- listAllPatients[which(as.matrix(Processedlogical) == TRUE)]
  
  return(listAllPatients)
}


###### Other Functions ######

DP_CalculateTimelimits <- function( sub_pat ,  HoursBeforeandAfter){
  timeindex<- c(0,0)
  timeindex[1] <- DP_AddHour(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1] ) , -HoursBeforeandAfter$numberhoursbefore)
  timeindex[2] <- DP_AddHour(DP_StripTime(sub_pat$ConfirmedFirstNewAF[1] ) , HoursBeforeandAfter$numberhoursafter)
  return(timeindex)
}

DP_FindNumberUniques <- function(X){
  # Function to find number of unique values in a column vector.
  values <- unique(X)
  n <- matrix(0 , length(values) , 1)
  nn <- matrix(0 , length(values) , 1)
  for(i in 1:length(values))
  {
    n[i] <- sum( X == values[[i]])  
    
    if(is.factor(values[[i]])){ nn[i] <- as.character(values[[i]])} else {nn[i] <- values[[i]]}
  }
  output <- data.frame( nn , n)
  output <- setNames(output , c('values' , 'n'))
  return(output) 
}
DP_StripTime <- function(X){
  if(is.POSIXct(X)){output <- X}
  if(!is.POSIXct(X))
  {
  output <- strptime(X ,  format = "%d/%m/%Y %H:%M")
  if(is.na(output)){output <- strptime(X ,  format = "%d/%m/%Y %H:%M:%S")}
  if(is.na(output)){output <- strptime(X ,  format = "%Y-%m-%d %H:%M:%S")}
  }
return(output)    
} 
DP_numberhoursgreaterthan <- function(PatIndex2017 , PatientsId , numberofhours = 6){
  sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
  output <- sub_pat$TotalITUTimeHRS[1] >= numberofhours
  if(is.na(output)){output = FALSE}
  return(output)    
}
DP_ReturnPatientNumber <- function(PatientsId){
  return(as.numeric(substr(PatientsId , start = 2 , stop = min(nchar(PatientsId)  , 6) )))
}
DP_AlignRegionofInterests <- function(Waveform1 , Waveform2 , regionofinterest){
  startindex <- which.min( abs(Waveform2[ , 1] - Waveform1[ regionofinterest[1] , 1]) )
  endindex <- startindex + length(regionofinterest)
  regionofinterest2 <- startindex:endindex
  return(regionofinterest2)
}
DP_ValidateRPeaks<-function(RPeaksOutput){
  if(sum(names(RPeaksOutput) == c("ECGI","ECGII","ECGIII","Meta_Data","RRCombined")) ==5 || sum(names(RPeaksOutput) == c("ECGI","ECGII","ECGIII","MetaData","RRCombined")) ==5){
  return(TRUE)  
  }else{return(FALSE)}
}
DP_ExtractPatientRecordforIndex <- function(PatIndex2017 , PatientCode){
  return( subset(PatIndex2017, PseudoId %in% PatientCode))
}
DP_SaveFile <- function( object , path , PatientID , Name){
  save( object , file = paste0(path , '\\' , PatientID , '\\Zip_out\\', Name , '.RData') )
}
DP_WaitBar <- function(A){
  print(paste0(A*100 , '% complete'))
}  
DP_CreateDummyMetaData <- function(PatIndex2017 , Name = NA , FirstNewAF = NA){
  output <- setNames(data.frame(NA , NA , NA , NA , NA , NA , NA , NA , NA , NA,
                                NA , NA , NA , NA , NA , NA , NA , NA , NA , NA,
                                NA , NA , NA , NA , NA , NA , NA , NA , NA , NA,
                                NA , NA , NA , NA , NA , NA , NA , NA , NA , NA,
                                NA , NA , NA , NA , NA , NA , NA) , names(PatIndex2017))
  output$TotalITUTimeHRS <- 80
  output$Pre_OperativeHeartRhythm <- "Sinus Rhythm"
  output$Usable <- 1
  output$PseudoId <- Name
  output$ConfirmedFirstNewAF <- FirstNewAF
  output$FirstNewAF <- FirstNewAF
  return(output)
}
DP_ExtractNumberofRecordsSingle <- function(Time , index  , interval = 1 ){
return(length(Time[ (Time < (Time[index] + interval ))*(Time > (Time[index] - interval) ) == 1  ]  ))  
}
DP_ExtractNumberofRecord <- function( Time , interval = 1 ){
dt <- diff(Time)
timeincriment <- round(median(dt) , digits = 6)  
output <- (2*interval/0.005)*smth( c(1,round(dt/0.005 , digits = 4) < 1.2  ) , method = 'sma' , n = (2*interval/0.005) )
output[is.na(output)] <- (2*interval/0.005)
return(output)
}
DP_AddHour <- function(X , hours){
  return(seq.POSIXt( from=X , by=paste0(as.character(hours) , " hour"), length.out=2 )[2])
}
DP_ProcessRpeaksMultipleECGs <- function(ECGs , PatientRecord = DP_CreateDummyMetaData(PatIndex2017  , Name = subList)  , Filter =  wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.5 , timethresh = 0.02){
  outputdata <- list()
  outputdata[[1]] <- CleanRpeaks(RPeakExtractionWavelet( ECGs$ECGI   , Filter  , nlevels  , ComponetsToKeep  , stdthresh ) , 2)
  outputdata[[2]] <- CleanRpeaks(RPeakExtractionWavelet( ECGs$ECGII  , Filter  , nlevels  , ComponetsToKeep  , stdthresh ) , 2)
  outputdata[[3]] <- CleanRpeaks(RPeakExtractionWavelet( ECGs$ECGIII , Filter  , nlevels  , ComponetsToKeep  , stdthresh ) , 2)
  outputdata[[4]] <- PatientRecord
  outputdata[[5]] <- 1
  outputdata <- setNames( outputdata , c('ECGI' ,'ECGII' ,'ECGIII' , 'MetaData' , 'RRCombined') )
  outputdata[[5]] <- PE_MultipleECGRPeaks(outputdata , ECGs , thresh = timethresh)
  return(outputdata)
}
DP_RemoveNaRows <- function(X){
  X <- as.matrix(X)
  return(as.matrix(X[DP_FindNARows(X) , ]))

}
DP_FindNARows <- function(X){
  return(apply(is.na(X) , 1 , function(Y){sum(Y) == 0} ))
}
DP_FindNumzeroRows <- function(X){
  return(apply(X , 1 , function(Y){sum(Y == 0)} ))
}
DP_FindZeroVarianceRows <- function(X){
  return(apply(X , 2 , function(Y){var(Y[!is.na(Y)])} ) == 0 )
}
DP_fixculmulativeprobs <- function(X){
  for(i in 1:(length(X) -1) ){
    if(is.na(X[i])){next}
    if(X[i] > X[i+1]){
      X[i+1] <- X[i]
    }
  }
  X[length(X)] <-1
  return(X)
}
DP_SampleRowsFromMatrix <- function(X , numberofsamples = 1000){
  X <- as.matrix(X)
  return(X[sample( 1:dim(X)[1] , numberofsamples ) , ])
}
DP_RescaleZeroOneToab <- function(X, a , b){
  return(X*(b-a) + a)
}
DP_pdist <- function(X , Xstar){
  return(as.matrix(pdist(X , Xstar)))
}
DP_FindMindistances <- function(X){
  dismatrix <- DP_pdist(X , X)
  diag(dismatrix) <- Inf
  return( apply(dismatrix , 1, min))  
}
DP_FindMeanMindistances <- function(X ){
  return(mean(DP_FindMindistances(X)))
  
}
DP_CalculateLimits <- function(X){
  return(c(min(X) , max(X)))
}
DP_RemoveEmptyElementsfromlist <- function(X ){
  baddataindicies <- which(unlist(lapply(X , length))==0)
  while(length(baddataindicies) > 0 ){
    X[[baddataindicies[1]]] <- NULL
    baddataindicies <- which(unlist(lapply(X , length))==0)
    
  }
  return(X)
}
DP_NormaliseData <- function(X){
  return( (X - mean(X[!is.na(X)]))/(sqrt(var(X[!is.na(X)]))) )
}
##### Quality of life functions ######
size <- function( X ){
  dim(X)}
is.POSIXct <- function( X ){ 
  inherits(X, "POSIXct")}
is.POSIXlt <- function( X ){ 
  inherits(X, "POSIXlt")}
is.POSIXt <- function( X ){ 
  inherits(X, "POSIXt")}
is.Date <- function( X ){ 
  inherits(X, "Date")}
isodd <- function( A ){
  return( (A %% 2) == 0)
}
disp <- function(A){
  print(A)
}

DP_cummean <- function(X){
  return( cumsum(X)/cumsum(rep(1  , length(X))) )
}
DP_cumvar <- function(X){
  cummu = DP_cummean(X)
  return( cumsum( (X - DP_cummean(X))^2)/(cumsum(rep(1  , length(X)) )- 1) )
}
DP_AddNugget <- function(X , nugget = 0.000000000001){
  if(is.matrix(nugget) == FALSE){
    return(X + nugget*diag(dim(X)[1])  )}
  if(is.matrix(nugget) == TRUE){
    return(X + nugget  )}
}
DP_CheckProprotionofNa <- function(X){
  return( sum(is.na(X))/length(X) )
}
DP_SortMatrix <- function(X , column = 1){
  # function to sort matrix by column i
  X <- X[order(X[,column]) , ]
}

DP_ExtractMetaDataForMultiplePatinets<- function(PatIndex2017 , listAllPatients){
  return(PatIndex2017[PatIndex2017[ , which(names(PatIndex2017) == 'PseudoId')] %in% listAllPatients , ])  
}
DP_ExtractRecordsFromDendrite <- function(DetIndex2017 , listAllPatients){
  
  LogicalVector <- matrix(0 , dim(DetIndex2017)[1] , 1)
  for( i in 1:dim(LogicalVector)[1] ){
    if( sum(apply(as.matrix(listAllPatients) ,1, function(X){ grepl( X , DetIndex2017$NewPseudoId[i]) })) >0 ){
      LogicalVector[i] <- 1
      DetIndex2017$NewPseudoId[i] <- listAllPatients[which(apply(as.matrix(listAllPatients) ,1, function(X){ grepl( X , DetIndex2017$NewPseudoId[i] , fixed = TRUE) }))][length(listAllPatients[which(apply(as.matrix(listAllPatients) ,1, function(X){ grepl( X , DetIndex2017$NewPseudoId[i] , fixed = TRUE) }))])]
    }
  }
  return(DetIndex2017[which(LogicalVector == 1) , ])
}
DP_RestructureBioChem <- function(BioChemIndex2017){
  listoftsvariables <- names(BioChemIndex2017)[15:22]
  
  uniquenames <- unique(BioChemIndex2017$NewPseudoId)
  uniquenames <- uniquenames[!is.na(uniquenames)]
  
  NewData <- list()
  for( i in 1:length(uniquenames) ){
    NewData[[i]] <- setNames(list(1 , 1) , c('TimeSeriesData' , 'MetaData'))
    NewData[[i]]$TimeSeriesData <-  data.frame(time = as.POSIXct(DP_StripTime(BioChemIndex2017$PostOpUsandEsTime[grepl(uniquenames[i] , BioChemIndex2017$NewPseudoId) ]) ) , tsdata <- BioChemIndex2017[grepl(uniquenames[i] , BioChemIndex2017$NewPseudoId)  , 15:22] ) 
    NewData[[i]]$MetaData <- data.frame(BioChemIndex2017[which(grepl(uniquenames[i] , BioChemIndex2017$NewPseudoId))[1]  , -c(15:22)])
  }
  NewData <- setNames(NewData , uniquenames)
  return(NewData)
}
DP_RestructureFluids <- function(FluidsIndex2017){
  listoftsvariables <- names(FluidsIndex2017)
  listoftsvariables <- listoftsvariables[-c(1 , 2 , 3, 4 ,5 ,6 ,7 , 8 , 9 , 10 , 187 , 188)]
  listoftsindexes  <-   which(names(FluidsIndex2017) %in% listoftsvariables)
  
  uniquenames <- unique(FluidsIndex2017$NewPseudoId)
  uniquenames <- uniquenames[!is.na(uniquenames)]
  
  NewData <- list()
  for( i in 1:length(uniquenames) ){
    NewData[[i]] <- setNames(list(1 , 1) , c('TimeSeriesData' , 'MetaData'))
    NewData[[i]]$TimeSeriesData <-  data.frame(time = as.POSIXct(DP_StripTime(FluidsIndex2017$Result.DT[grepl(uniquenames[i] , FluidsIndex2017$NewPseudoId) ]) ) , tsdata <- FluidsIndex2017[grepl(uniquenames[i] , FluidsIndex2017$NewPseudoId)  , listoftsindexes] ) 
    }
  NewData <- setNames(NewData , uniquenames)
  return(NewData)
}
DP_RestructureVent <- function(VentIndex2017){
  listoftsvariables <- names(VentIndex2017)[4:15]
  
  uniquenames <- unique(VentIndex2017$PseudoId)
  uniquenames <- uniquenames[!is.na(uniquenames)]
  
  NewData <- list()
  for( i in 1:length(uniquenames) ){
    NewData[[i]] <- setNames(list(1 ) , c('TimeSeriesData'))
    NewData[[i]]$TimeSeriesData <-  data.frame(time = as.POSIXct(DP_StripTime(VentIndex2017$Result.DT[grepl(uniquenames[i] , VentIndex2017$PseudoId) ]) ) , tsdata <- VentIndex2017[grepl(uniquenames[i] , VentIndex2017$PseudoId)  , 4:15] ) 
  }
  NewData <- setNames(NewData , uniquenames)
  return(NewData)
}
DP_ExtractPatientIDFromNewPatinetID <- function( newpatientID ){
  return(  substr(x = newpatientID , start = 1 , stop =  gregexpr( '-' , newpatientID )[[1]][1] -1 )  ) 
}