subList = select.list(listAllPatients
                      ,preselect = NULL
                      , multiple = FALSE
                      ,title = 'Select Patient to Analyse'
                      , graphics = TRUE )

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