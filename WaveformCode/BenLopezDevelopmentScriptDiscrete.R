pathFiles <- choose.dir(caption="Select folder with source code")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

# Load settings
source("LibrariesAndSettings.R" , print.eval  = TRUE )

# load data % Tested on z1139.Rdata
FileLocation <- choose.files(caption="Select .Rdata file of cleaned ECG data")
load(FileLocation)

HeartRate <- DiscDataTrim$Value[DiscDataTrim$VarName == DiscDataTrim$VarName[3]]
tt <- DiscDataTrim$Date[DiscDataTrim$VarName == DiscDataTrim$VarName[3]]
  
plot( tt , HeartRate  , xlab = 't' , ylab = 'Heartrate')
title('z812 Discrete Heartrate')
abline(60,0) 
abline(100,0 , col = 'red') 
abline(80 , 0 , col = 'blue') 

