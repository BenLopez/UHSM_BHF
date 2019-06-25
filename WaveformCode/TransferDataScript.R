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

pathtocopyfrom = choose.dir( caption = "Select folder to copy from")
Listofiles <- select.list(list.files(pathtocopyfrom) , graphics = TRUE , multiple = T , title = "Select patients to copy")
pathtocopyto = choose.dir(caption = "Select folder to copy to")

# copy processed files
for( ii in 1:length(Listofiles) ){
  dir.create(paste0(pathtocopyto , '\\' , Listofiles[[ii]] ) ,  recursive = TRUE)
  file.copy( from = paste0(pathtocopyfrom , '\\' , Listofiles[[ii]] , '\\Zip_out' ) , to = paste0(pathtocopyto , '\\' , Listofiles[[ii]] ) ,  recursive = TRUE)
  print(paste0('Percentage complete ' , (ii/length(Listofiles))*100 , '%'))
}

# copy files that have not been processed
for( ii in 1:length(Listofiles) ){
if(DP_checkfilesprocessed(path , Listofiles[[ii]] , 'ECGII') ==1  ){
  next}else{
    dir.create(paste0(pathtocopyto , '\\' , Listofiles[[ii]] ) ,  recursive = TRUE)
    file.copy( from = paste0(pathtocopyfrom , '\\' , Listofiles[[ii]] ) , to = paste0(pathtocopyto  ) ,  recursive = TRUE)
    print(paste0('Percentage complete ' , (ii/length(Listofiles))*100 , '%')) 
  }
}