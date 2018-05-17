# Select directory containing source
pathFiles = choose.dir(caption="Select folder with source code")
pathFiles = paste0(pathFiles, "\\")

##### Select file PatientIndex.csv if it exists

filetype = select.list(c('csv' , 'RData'), preselect = NULL, multiple = TRUE,
            title = 'Choose File Type For Patient Index', graphics = TRUE )

if(filetype == 'csv')
{
path_PatIndex = choose.files(caption="Select 2017 PatientIndex.csv file")

if(length(path_PatIndex)>0){
  PatIndex2017 = read.csv(file=path_PatIndex, stringsAsFactors = FALSE)
  # PatIndex2017$FirstITUEntry=as.POSIXct(PatIndex2017$FirstITUEntry, format="%d/%m/%Y %H:%M")
  # PatIndex2017$LastITUEntry=as.POSIXct(PatIndex2017$LastITUEntry, format="%d/%m/%Y %H:%M")
} else {
  warning("No Patient Info provided")
  sub_pat = list()
}
}
if(filetype = 'RData')
{
  path_PatIndex =  choose.files()
  load(path_PatIndex)
}

choose_outputs = c(0,0,1) #csv, mat, rdata --- DO NOT USE MAT, needs testing, very slow
## Install R.matlab library and uncomment next line if you want to generate mat files
# library(R.matlab)

# Option to choose maximium number of hours to process to prevent memory bottlenecks
maxhourstoprocess <- 100

Use7z = 1
UseZip = 0
KeepCSVs = 0


DataTypes = c("Discrete", "ECGI", "ECGII", "ECGIII", "CVP", "ART", "SPO2", "Flow", "Paw")
chooseWave2Read = select.list(DataTypes, preselect = DataTypes,
                              multiple = TRUE, graphics = TRUE, title = "Choose Waves to Read")

# Run functions source file
source(paste0(pathFiles,"/sourceFunctions.R"))

# Select directory containing all patients
path = choose.dir(caption="Select folder containing data repository")
#setwd(path)
listAllPatients = list.dirs(path = path, full.names = FALSE, recursive = FALSE)

#
subList = select.list(listAllPatients, preselect = NULL, multiple = TRUE,
                      title = NULL, graphics = TRUE )

# path = paste0(path, "\\")
# findZ = regexpr("z[^z]*$", path)[1]
# PatientCode = substr(path, findZ, nchar(path)-1)
for (PatientCode in subList){
  if(length(path_PatIndex)>0)
  {
    sub_pat = subset(PatIndex2017, PseudoId %in% PatientCode)
    # warning("Patient not in PatientIndex") # not sure this warning is in the right place
    if(nrow(sub_pat) == 0){
      sub_pat = list()
      print(paste0("Patient ", PatientCode , " not in PatientIndex: skipping to next patient"))
      next }
  } 
  else 
  {
    error("No Patient Info provided")
    # sub_pat = list()
  }
  if(sub_pat[["TotalITUTimeHRS" ]] > maxhourstoprocess)
  {
    print(paste0(PatientCode , " skipped as over max hours"))
    next}
  
  # Added redundacy to stop memory overload.
  if(sub_pat[["TotalITUTimeHRS" ]] < maxhourstoprocess)
  {
    print(paste0("Processing ",PatientCode))
    source(paste0(pathFiles,"/readPatient.R"))
  }
  
  if(KeepCSVs == 0)
  {
    for(PatientCode in subList)
    {
      pathIn = path
      unlink(paste0(pathIn,"\\Disc_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ECGI_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ECGII_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ECGIII_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\CVP_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ART_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\SPO2_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\Flow_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\Paw_clean"), recursive = TRUE)
      unlink(paste0(pathIn,"\\Disc_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ECGI_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ECGII_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ECGIII_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\CVP_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\ART_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\SPO2_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\Flow_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\Paw_cleanMAT"), recursive = TRUE)
      unlink(paste0(pathIn,"\\",PatientCode,"\\temp_zip"), recursive = TRUE)
    }
  }
  
}
