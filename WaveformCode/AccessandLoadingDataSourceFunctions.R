# Script with WaveForm data processing and access functions.


DP_LoadPatientIndex <- function()
{
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
}

DP_ChooseDataReps <- function( )
{
  
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

DP_choosepatient <- function(listAllPatients)
{
  subList <<- select.list(listAllPatients
                        , preselect = NULL
                        , multiple = FALSE
                        , title = 'Select Patient to Analyse'
                        , graphics = TRUE )
}

DP_choosepatients <- function(listAllPatients)
{
  subList <<- select.list(listAllPatients
                          , preselect = NULL
                          , multiple = TRUE
                          , title = 'Select Patient to Analyse'
                          , graphics = TRUE )
}

DP_LoadECG <- function(path , subList , numberrep=1 , ECGNum = 1 )
{
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
  
DP_SelectInterestingTimePoint <- function(Wavedata , sub_pat)
{
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

DP_SelectHoursBeforeandAfter <- function()
{ 
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

DP_CropWaveData <- function(WaveData , timeindex , HoursBeforeAndAFter)
{
  numberhoursbefore <- HoursBeforeAndAFter[['numberhoursbefore']]
  numberhoursafter <- HoursBeforeAndAFter[['numberhoursafter']]
  indent <- as.numeric(abs(WaveData[1 , 1] - WaveData[2 , 1]))
  if(indent > 0.012){warning('Time indent greater than 0.012. Data measured at unusual redolution.')}
  WaveData <- WaveData[ max( 1 , timeindex - (numberhoursbefore*((60^2)/indent)) ) : min(length(WaveData[ , 1]) , timeindex + (numberhoursafter*((60^2)/indent)) ) , ]
  WaveData <- ReturnWaveformwithPositiveOrientation(WaveData)
  return(WaveData)
}

DP_FindNumberUniques <- function(X)
{
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

DP_StripTime <- function(X)
{
  if(is.POSIXct(X)){output <- X}
  if(!is.POSIXct(X))
  {
  output <- strptime(X ,  format = "%d/%m/%Y %H:%M")
  if(is.na(output)){output <- strptime(X ,  format = "%d/%m/%Y %H:%M:%S")}
  }
return(output)    
} 

DP_ChooseECGstoProcess <- function( )
{
  DataTypes = c("ECGI", "ECGII", "ECGIII")
  chooseWave2Read = select.list(DataTypes, preselect = DataTypes,
                                multiple = TRUE, graphics = TRUE, title = "Choose Waves to Read")
  return(chooseWave2Read)
}

DP_ChooseWaveformstoProcess <- function( )
{
  DataTypes = c("Discrete", "ECGI", "ECGII", "ECGIII", "CVP", "ART", "SPO2", "Flow", "Paw")
  chooseWave2Read = select.list(DataTypes, preselect = DataTypes,
                                multiple = TRUE, graphics = TRUE, title = "Choose Waves to Read")
  return(chooseWave2Read)
}

DP_checkfilesprocessed <- function(path , PatientsId , FilestoProcess)
{
  output <- rep(0 , length(FilestoProcess) , 1)
  for( i in 1:length(FilestoProcess) )
  {
    output[i] <- file.exists(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , FilestoProcess[i] , '_' , PatientsId , '.RData')  )
  }
return(output)
}

DP_existsinpatientindex <- function(PatIndex2017 , PatientsId)
{
return(nrow(subset( PatIndex2017, PseudoId %in% PatientsId )) > 0)  
}

DP_isusable <- function(PatIndex2017 , PatientsId)
{
sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
return(sub_pat$Usable[1] == 1)    
}

DP_numberhoursgreaterthan <- function(PatIndex2017 , PatientsId , numberofhours = 6)
{
  sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
  output <- sub_pat$TotalITUTimeHRS[1] >= numberofhours
  if(is.na(output)){output = FALSE}
  return(output)    
}

DP_ReturnPatientNumber <- function(PatientsId)
{
  return(as.numeric(substr(PatientsId , start = 2 , stop = min(nchar(PatientsId)  , 6) )))
}

DP_FilterbySinusRhythum <- function(PatIndex2017 , PatientsId)
{
  sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
  return(sub_pat$Pre_OperativeHeartRhythm[1] == 'Sinus Rhythm')
}

DP_FilterbyOps <- function(PatIndex2017 , PatientsId , HowtoFilterops)
{
  sub_pat <- subset( PatIndex2017, PseudoId %in% PatientsId )
  
  index <- lapply(sub_pat$ProcDetails , function(X){which(X == HowtoFilterops$Optype)} )
  keep <- as.matrix(lapply(index , function(X){DP_Filtbyopsindividualrecord(X , HowtoFilterops = HowtoFilterops)}))
  return( sum(keep == 0) == 0)
}

DP_Filtbyopsindividualrecord <- function(index, HowtoFilterops)  
{
  if(!is.na(HowtoFilterops$Keep[index]) ){keep = 1}
  if(!is.na(HowtoFilterops$Remove.for.this.op[index])){keep = 1}
  if(!is.na(HowtoFilterops$Removal.all.data.if.they.ever.have.this.op[index])){keep = 0}  
  return(keep)
}  
  
DP_FilterPatients<-function(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
{
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

DP_checkRpeaksfilesprocessed<- function(path , PatientsId )
{
  
  output <- file.exists(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , PatientsId , '_RPeaks' , '.RData')  )
  return(output)
}

DP_CheckECGreducedfilesprocessed <- function( path , PatientsId  , Filestoprocess)
{
  output <- file.exists(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , PatientsId , '_', Filestoprocess , '.RData')  )
  return(output)
}

DP_LoadECGReduced <- function(path , subList , numberrep=1 , ECGNum = 1 )
{
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

DP_LoadRpeaksfileECGI <- function(path , PatientsId )
{
  load(paste0(path , '\\' , PatientsId , '\\Zip_out\\' , PatientsId , '_RPeaks' , '.RData'))
  return(outputdata$ECGI)
}  


is.POSIXct <- function(X){ inherits(X, "POSIXct")}
is.POSIXlt <- function(X){ inherits(X, "POSIXlt")}
is.POSIXt <- function(X){ inherits(X, "POSIXt")}
is.Date <- function(X){ inherits(X, "Date")}

