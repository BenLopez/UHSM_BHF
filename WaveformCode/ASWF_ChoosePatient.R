subList = select.list(listAllPatients
                      , preselect = NULL
                      , multiple = FALSE
                      , title = 'Select Patient to Analyse'
                      , graphics = TRUE )

if(DP_checkfilesprocessed(path , subList , 'Discrete') == 1){  
# Load discrete data
for(i in 1:(numberrep+1))
{
  if(i > (numberrep))
  {
    print('Error: No discrete data processed. Choose another patient or process data.') 
    source('ASWF_ChoosePatient.R')
    break
  }
  if(file.exists(paste0( path[[i]] , '\\' , subList , '\\Zip_out\\' , 'Discrete_' , subList , '.RData' )))
  {
    DataSet <- LoadSingleDiscreteDataHeartRates( subList , path[[i]] )
    break 
  }  
}
}
if(DP_checkfilesprocessed(path , subList , 'Discrete') == 0){warning('No discrete data processed.')}

PatientRecord <- DP_ExtractPatientRecordforIndex(PatIndex2017  , subList)

if(nrow(PatientRecord) > 0){
if( PatientRecord$TotalITUTimeHRS > 100 & DP_checkRpeaksfilesprocessed(path ,  subList ) == FALSE ){print(paste0('Total hours over 100'))
  source('ASWF_ChoosePatient.R')}}
