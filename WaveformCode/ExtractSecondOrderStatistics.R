

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


E_medRR <- list()
Std_medRR <- list()
E_medRA <- list()
Std_medRA <- list()
NEWAFLogical<-list()
tt <- list()
ID <-list()
counter <- 1

for( PatientCode in subList )
{
  
  # Check ECG1Rdata has been prcossed
  if(!file.exists( paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'ECGI_', PatientCode , '.RData') )){
    print(paste0(PatientCode , ' has no waveform.rData skipping to next file'))
    next
  }
  # Check the Rpeaks file has not already been produced.
  if(!file.exists( paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'ECGI_', PatientCode , '_Rpeaks' , '.RData') )){
    print(paste0(PatientCode , ' already processed skipping to next file'))
    next
  }
  
  sub_pat = subset(PatIndex2017, PseudoId %in% PatientCode)
  
  # Load wavedata 
  print( paste0('Processing ' , PatientCode) )
  load( paste0(path , '\\' , PatientCode , '\\' , 'Zip_out\\' , 'ECGI_', PatientCode , '_Rpeaks' , '.RData'))
  
  # Extract hours before end.
  if(is.na(sub_pat$FirstNewAF)  )
  {
    tt[[counter]] <- RWaveExtractedData$t
    Filter2 = rep(1 , 1200)/1200
    RWaveExtractedData$RR[RWaveExtractedData$RR > 2] <- 0.8
    E_medRR[[counter]] <- imfilter1D(as.numeric(RWaveExtractedData$RR) , Filter2 )
    Std_medRR[[counter]] <- sqrt(imfilter1D(as.numeric(RWaveExtractedData$RR)^2 , Filter2) - E_medRR[[counter]]^2)
    E_medRA[[counter]] <- imfilter1D(as.numeric(RWaveExtractedData$RA) , Filter2 )
    Std_medRA[[counter]] <- sqrt(imfilter1D(as.numeric(RWaveExtractedData$RA)^2 , Filter2) - E_medRA[[counter]]^2)
    NEWAFLogical[[counter]] <- 0       
    ID[[counter]]<-PatientCode
  }
  
  if(is.na(sub_pat$FirstNewAF) == FALSE  )
  {
    tt[[counter]] <- RWaveExtractedData$t
    RWaveExtractedData$RR[RWaveExtractedData$RR > 2] <- 0.8
    Filter2 = rep(1 , 1200)/1200
    E_medRR[[counter]] <- imfilter1D(as.numeric(RWaveExtractedData$RR) , Filter2 )
    Std_medRR[[counter]] <- sqrt(imfilter1D(as.numeric(RWaveExtractedData$RR)^2 , Filter2) - E_medRR[[counter]]^2)
    E_medRA[[counter]] <- imfilter1D(as.numeric(RWaveExtractedData$RA) , Filter2 )
    Std_medRA[[counter]] <- sqrt(imfilter1D(as.numeric(RWaveExtractedData$RA)^2 , Filter2) - E_medRA[[counter]]^2)
    NEWAFLogical[[counter]] <- 1 
    ID[[counter]]<-PatientCode
  }
  
  counter <- counter +1  
  
}

num_points <- matrix(0,length(NEWAFLogical) , 1)
var_tt <- num_points
difftt <- num_points

for( i in 1:length(NEWAFLogical))
{
  num_points[i] = length(tt[[i]])
  var_tt[i] = var(tt[[i]])
  difftt[i] =  tt[[i]][length(tt[[i]])] - tt[[i]][1]
}

difftt[difftt<0]<-100

dev.off()
par(mfrow=c(2 , 1))
plot(tt[[1]] - tt[[1]][1] , E_medRR[[1]] , type ='l' , ylim = c(0.3,1.25) , xlim = c(2500,20000) , ylab = 'local mean RR' , xlab = 't' )
for(i in 2:length(NEWAFLogical))
{
if(difftt[i] > 7  ){next}
if(length(tt[[i]]) < 10000){next}
if(length(tt[[i]]) > 35000){next}
if(NEWAFLogical[[i]] == 0){ lines(tt[[i]] - tt[[i]][1] , E_medRR[[i]]) }  
if(NEWAFLogical[[i]] == 1){ lines(tt[[i]] - tt[[i]][1] , E_medRR[[i]] , col = 'red') }  
}  
title('Local Mean RR times for 168 Patients')


j <- 1
plot(tt[[1]] - tt[[1]][1] , Std_medRR[[1]] , type ='l' , ylim = c(0,0.5) , xlim = c(2500,20000), ylab = 'local variance RR' , xlab = 't')
for(i in 2:length(NEWAFLogical))
{
  if(difftt[i] > 7  ){next}
  if(length(tt[[i]]) < 10000){next}
  if(length(tt[[i]]) > 35000){next}
  if(NEWAFLogical[[i]] == 0){ lines(tt[[i]] - tt[[i]][1] , Std_medRR[[i]]) }  
  if(NEWAFLogical[[i]] == 1){ lines(tt[[i]] - tt[[i]][1] , Std_medRR[[i]] , col = 'red') }  
}  

title('Local Variance RR times for 168 Patients')
j<-j+1

