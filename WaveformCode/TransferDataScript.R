{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

pathtocopyfrom = choose.dir( caption = "Select folder to copy from")
Listofiles <- select.list(list.files(pathtocopyfrom) , graphics = TRUE , multiple = T , caption = "Select patients to copy")
pathtocopyto = choose.dir(caption = "Select folder to copy to")

for(ii in 1:length(Listofiles)){
  dir.create(paste0(pathtocopyto , '\\' , Listofiles[[ii]] ) ,  recursive = TRUE)
  file.copy( from = paste0(pathtocopyfrom , '\\' , Listofiles[[ii]] , '\\Zip_out' ) , to = paste0(pathtocopyto , '\\' , Listofiles[[ii]] ) ,  recursive = TRUE)
  print(paste0('Percentage complete ' , (ii/length(Listofiles))*100 , '%'))
}