zip_files = list.files(path = paste0(path,"\\",PatientCode),pattern="*\\.zip*$")
filesInZip = unlist(lapply(zip_files,function(filename){unzip(paste0(path,"\\",PatientCode,"\\",filename), list = TRUE)$Name}))

dir.create(path = paste0(path,"\\",PatientCode,"\\temp_zip"), showWarnings = FALSE)

# Clear folder
do.call(file.remove, list(list.files(paste0(path,"\\",PatientCode,"\\temp_zip"), full.names = TRUE)))

for (zipF in zip_files){
  unzipfile(paste0(zipF,"unzip.vbs",sep=""),zipF,paste0(path,"\\",PatientCode,"\\"))
}
####### Wait for files to unzip before proceeding! ##########
while (length(setdiff(filesInZip, list.files(paste0(path,"\\",PatientCode,"\\temp_zip\\"))))>0) {
  Sys.sleep(1)
}

####
check_text = list.files(path = paste0(path,"\\",PatientCode,"\\"),pattern=".txt*$")

csv_files = list.files(path=paste0(path,"\\",PatientCode,"\\temp_zip"),pattern=".csv")
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

if(length(Disc_files)>0 & ("Discrete" %in% chooseWave2Read)){
  print("Reading Discrete Data")
  source(paste0(pathFiles,"/sourceDiscrete.R"), echo = TRUE)
}

# ECG Data ----------------------------------------------------------------

options(digits = 15)
options(digits.secs=3)

source(paste0(pathFiles,"/sourceReadWrite.R"))
pathIn = paste0(path,"\\",PatientCode,"\\")
options(warn = -1)

#### Lead I
if (length(ECGI_files)>0 & ("ECGI" %in% chooseWave2Read)){
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
if (length(ECGII_files)>0 & ("ECGII" %in% chooseWave2Read)){
readWriteWave(ECGII_files,
              cleanWave,
              sub_pat,
              choose_outputs,
              "ECGII", pathIn, pathOutFilesExtra, Use7z, UseZip)
}
#### Lead 3
if (length(ECGIII_files)>0 & ("ECGIII" %in% chooseWave2Read)){
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