
pathZIPs = paste0(path,"\\",PatientCode,"\\Zip_out\\")

# Check if data has already been processed.
numberofprocessedfiles <- 0
for(i in 1:length(chooseWave2Read)){
  numberofprocessedfiles <- numberofprocessedfiles + as.numeric(file.exists(paste0( pathZIPs , chooseWave2Read[i] , '_' ,  PatientCode , '.RData') ) )
}

if(numberofprocessedfiles == length(chooseWave2Read)){
print(paste0('All ECG files processed for' , PatientCode , ' moving to next patient'))
next}

zip_files = list.files(path = paste0(path,"\\",PatientCode),pattern="*\\.zip*$")


objStr = c("discrete","ECG I.W","ECG II.W", "ECG III.W", "CVP.W", "ART.W", "SPO2.W", "Flow.W","Paw.W")
objTypes = c("Discrete","ECGI","ECGII","ECGIII","CVP","ART","SPO2","Flow","Paw")
names(objStr) = objTypes

filesInZip = unlist(lapply(zip_files,function(filename){paste0(filename,"/",unzip(paste0(path,"\\",PatientCode,"\\",filename), list = TRUE)$Name)}))
filesInZipW = filesInZip[rowSums(sapply(objStr[chooseWave2Read],function(x) {grepl(x, filesInZip)}))==1]
sizesInZip = unlist(lapply(zip_files,function(filename){unzip(paste0(path,"\\",PatientCode,"\\",filename), list = TRUE)$Length}))
sizesInZipW = sizesInZip[rowSums(sapply(objStr[chooseWave2Read],function(x) {grepl(x, filesInZip)}))==1]
dir.create(path = paste0(path,"\\",PatientCode,"\\temp_zip"), showWarnings = FALSE)

# Clear folder
do.call(file.remove, list(list.files(paste0(path,"\\",PatientCode,"\\temp_zip"), full.names = TRUE)))
sapply(list.files(paste0(path,"\\",PatientCode,"\\temp_zip"), full.names = TRUE), function(x){unlink(x,recursive = TRUE)})


# Create and execute shell script to unzip all zips in folder
for (zipF in zip_files){
  for (wave in chooseWave2Read){
    print(paste0("reading ", wave),sep="")
    unzipfileMatch(paste0(zipF,"unzip.vbs",sep=""),zipF,paste0(path,"\\",PatientCode,"\\"), objStr[wave])
  }    
}

####### Wait for files to unzip before proceeding! ##########
while (length(setdiff(filesInZipW,
                      list.files(paste0(path,"\\",PatientCode,"\\temp_zip\\"),recursive=TRUE)))>0 || 
       sum(file.info(list.files(paste0(path,"\\",PatientCode,"\\temp_zip\\"),
                                full.names = TRUE,
                                recursive = TRUE, include.dirs = FALSE))$size)<sum(sizesInZipW)) {
  Sys.sleep(1)
}
cnt = 0
for (zipF in zip_files){
  cnt = cnt + 1
  print(zipF)
  print(cnt)
  tpath = paste0(path,"\\",PatientCode,"\\temp_zip\\",zipF)
  filez = list.files(tpath)
  sapply(filez,FUN=function(eachPath){
    file.rename(from=paste0(tpath,"\\",eachPath),to=sub(pattern=".csv",replacement=paste0(cnt,".csv",sep=""),paste0(tpath,"\\",eachPath)))
  })
  filez = list.files(tpath, include.dirs = FALSE)
  sapply(filez,FUN=function(eachPath){
    file.move(paste0(tpath,"\\",eachPath),paste0(path,"\\",PatientCode,"\\temp_zip\\"))
  })
  Sys.sleep(1)
  filez = list.files(tpath, include.dirs = FALSE)
  sapply(filez,FUN=function(eachPath){
    file.move(paste0(tpath,"\\",eachPath),paste0(path,"\\",PatientCode,"\\temp_zip\\"))
  })
  Sys.sleep(1)
  unlink(tpath, recursive = TRUE)
}
print("moved")
####
check_text = list.files(path = paste0(path,"\\",PatientCode,"\\"),pattern=".txt*$", recursive = TRUE)

csv_files = list.files(path=paste0(path,"\\",PatientCode,"\\temp_zip"),pattern=".csv", recursive = TRUE)
#Discrete
Disc_files = csv_files[grepl("discrete",csv_files)]
#ECG I
ECGI_files = csv_files[grepl("ECG I.W", csv_files)]
#ECG II
ECGII_files = csv_files[grepl("ECG II.W", csv_files)]
#ECG III
ECGIII_files = csv_files[grepl("ECG III.W", csv_files)]
#CVP
CVP_files = csv_files[grepl("CVP.W", csv_files)]
#ART
ART_files = csv_files[grepl("ART.W", csv_files)]
#SPO2
SPO2_files = csv_files[grepl("SPO2.W", csv_files)]
#Flow
Flow_files = csv_files[grepl("Flow.W", csv_files)]
#Paw
Paw_files = csv_files[grepl("Paw.W", csv_files)]

#Path CSV files
pathCSV = paste0(path,"\\",PatientCode,"\\temp_zip\\")
pathZIPs = paste0(path,"\\",PatientCode,"\\Zip_out\\")

dir.create(pathZIPs, showWarnings = FALSE)

pathOutFilesExtra = paste0(path, "\\", PatientCode, "\\Extra_clean\\")
dir.create(pathOutFilesExtra, showWarnings = FALSE)

# Discrete Data -----------------------------------------------------------

if(length(Disc_files)>0 & ("Discrete" %in% chooseWave2Read) & file.exists(paste0( pathZIPs , 'Discrete_' ,  PatientCode, '.RData') ) == FALSE){
  print("Reading Discrete Data")
  source(paste0(pathFiles,"/sourceDiscrete.R"), echo = TRUE)
} else if (length(Disc_files)==0) {
  print("No discrete data")
}

# ECG Data ----------------------------------------------------------------
# some options for saving .RData files
options(digits = 15)
options(digits.secs=3)

source(paste0(pathFiles,"/sourceReadWrite.R"))
pathIn = paste0(path,"\\",PatientCode,"\\")
options(warn = -1)

#### Lead I
if (length(ECGI_files)>0 & ("ECGI" %in% chooseWave2Read) & file.exists(paste0( pathZIPs , 'ECGI_' ,  PatientCode, '.RData') ) == FALSE ){
  # Wave_files = ECGI_files
  # wavename = "ECGI"
  # path = pathIn

  readWriteWave(ECGI_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "ECGI",pathIn, pathOutFilesExtra, Use7z, UseZip)
}

# Wave_files = ECGI_files
# wavename = "ECGI"
# path = pathIn

#### Lead 2
if (length(ECGII_files)>0 & ("ECGII" %in% chooseWave2Read) & file.exists(paste0( pathZIPs , 'ECGII_' ,  PatientCode, '.RData')) == FALSE){
readWriteWave(ECGII_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "ECGII", pathIn, pathOutFilesExtra, Use7z, UseZip)
}
#### Lead 3
if (length(ECGIII_files)>0 & ("ECGIII" %in% chooseWave2Read) & file.exists(paste0( pathZIPs , 'ECGIII_' ,  PatientCode, '.RData')) == FALSE){
readWriteWave(ECGIII_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "ECGIII", pathIn, pathOutFilesExtra,Use7z,UseZip)

}
# CVP Files ---------------------------------------------------------------

if(length(CVP_files)>0 & ("CVP" %in% chooseWave2Read)){
readWriteWave(CVP_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "CVP",pathIn, pathOutFilesExtra,Use7z,UseZip)

}
# ART Files ---------------------------------------------------------------
if(length(ART_files)>0 & ("ART" %in% chooseWave2Read)){
readWriteWave(ART_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "ART",pathIn, pathOutFilesExtra,Use7z,UseZip)
}
# SPO2 Files ---------------------------------------------------------------
if(length(SPO2_files)>0 & ("SPO2" %in% chooseWave2Read)){
readWriteWave(SPO2_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "SPO2",pathIn,pathOutFilesExtra,Use7z,UseZip)
}

# Flow Files ---------------------------------------------------------------

if(length(Flow_files)>0 & ("Flow" %in% chooseWave2Read)){
readWriteWave(Flow_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "Flow",pathIn,pathOutFilesExtra,Use7z,UseZip)
}
# Paw Files ---------------------------------------------------------------

if(length(Paw_files)>0 & ("Paw" %in% chooseWave2Read)){
readWriteWave(Paw_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "Paw",pathIn,pathOutFilesExtra,Use7z, UseZip)
}
