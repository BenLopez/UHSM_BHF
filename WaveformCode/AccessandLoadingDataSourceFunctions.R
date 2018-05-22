# Script with WaveForm data processing and access functions.

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

DP_LoadECG <- function(path , subList , numberrep , ECGNum = 1 )
{
  # Function to load ECG data. ECGNum = (1 , 2 , 3) or ('I' , 'II' , 'III')
  if(ECGNum  == 1 || ECGNum  == 'I')
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
  
  if(ECGNum  == 2 || ECGNum  == 'II')
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
  
  if(ECGNum  == 3 || ECGNum  == 'III')
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
  
DP_SelectInterestingTimePoint <- function(WaveData , sub_pat)
{
  subdata <-  unique(as.vector(as.character(round.POSIXt(WaveData$Date , units = c('hours')))))   
  if(!is.na(as.POSIXct(sub_pat$FirstNewAF[1])))
  {
    interestingtimepoint <-  select.list(subdata
                                                           , preselect = subdata[ which.min( abs(as.POSIXct(sub_pat$FirstNewAF[1]) - as.POSIXct(subdata)))]
                                                           , multiple = FALSE
                                                           , title = 'Select Interest Time Point'
                                                           , graphics = TRUE )
    
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