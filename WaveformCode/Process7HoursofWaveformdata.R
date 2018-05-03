pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# Load in patientindex data
##### Select file PatientIndex.csv if it exists
path_PatIndex = choose.files(caption="Select 2017 PatientIndex.csv file")

if(length(path_PatIndex)>0){
  PatIndex2017 = read.csv(file=path_PatIndex, stringsAsFactors = FALSE)
  # PatIndex2017$FirstITUEntry=as.POSIXct(PatIndex2017$FirstITUEntry, format="%d/%m/%Y %H:%M")
  # PatIndex2017$LastITUEntry=as.POSIXct(PatIndex2017$LastITUEntry, format="%d/%m/%Y %H:%M")
} else {
  warning("No Patient Info provided")
  sub_pat = list()
}


path = choose.dir(caption="Select folder containing data repository")

#setwd(path)
listAllPatients = list.dirs(path = path, full.names = FALSE, recursive = FALSE)

#
subList = select.list(listAllPatients, preselect = NULL, multiple = TRUE,
                      title = NULL, graphics = TRUE )

HoursBeforeEnd <- 5
Filter <-  wt.filter(filter = "d6" , modwt=TRUE, level=1)


for( PatientCode in subList )
{
  
  # Check ECG1Rdata has been prcossed
  if(!file.exists( paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'ECGI_', PatientCode , '.RData') )){
    print(paste0(PatientCode , ' has no waveform.rData skipping to next file'))
    next
  }
  # Check the Rpeaks file has not already been produced.
  if(file.exists( paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'ECGI_', PatientCode , '_Rpeaks' , '.RData') )){
    print(paste0(PatientCode , ' already processed skipping to next file'))
    next
  }
  
  sub_pat = subset(PatIndex2017, PseudoId %in% PatientCode)
  
  # Load wavedata 
  print(paste0('Processing ' , PatientCode))
  load(paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'ECGI_', PatientCode , '.RData'))
  
  # Extract hours before end.
  if(is.na(sub_pat$FirstNewAF)  )
  {
    TimeofNAFEvent<-0
    index = round(length(WaveData$Date)/2)
    WaveData <- WaveData[  max(( index - ( HoursBeforeEnd*( (60^2) /0.005 ) ) ) , 1) : min(( index + (( (60^2) /0.005 ) ) )  , length(WaveData$Date)),  ]
  }
  
  if(is.na(sub_pat$FirstNewAF) == FALSE  )
  {
    TimeofNAFEvent <- as.POSIXct(strptime(sub_pat$FirstNewAF , "%d/%m/%Y %H:%M" ) )
    index <- which.min( abs(WaveData$Date - TimeofNAFEvent) )
    WaveData <- WaveData[  max(( index - ( HoursBeforeEnd*( (60^2) /0.005 ) ) ) , 1) : min(( index + (( (60^2) /0.005 ) ) )  , length(WaveData$Date)),  ]
  }
  
  WaveData <- ReturnWaveformwithPositiveOrientation(WaveData)
  RWaveExtractedData <- RPeakExtractionWavelet( WaveData , Filter , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8)
  
  print(paste0( PatientCode , ' processed; saving file.'))
  save('RWaveExtractedData' , file = paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'ECGI_', PatientCode , '_Rpeaks' , '.RData'))
  
}