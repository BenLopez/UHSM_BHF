{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )
DP_GetDirectories()}

PatientID <- DP_choosepatient( listAllPatients )

ECGs <- DP_LoadReducedECGs(path ,  PatientID , FilestoProcess = DP_ChooseECGstoProcess()  )
PeakData <- DP_LoadRpeaksfile(path , PatientID)

index <- 1
QS <- AFD_ExtractSQ( ECGs$ECGI , PeakData$RRCombined , index = 1)

QS <- AFD_ExtractAllSQ( ECGs$ECGI , PeakData$RRCombined[1:1000,] )

plot(0  , 0  , xlim = c(-0,0.7) , ylim = c(-2,10) , xlab='t' , ylab ='Hz' )
for(i in 1:length(QS)){
  lines(QS[[i]]$Date -QS[[i]]$Date[1], QS[[i]]$Value - mean(QS[[i]]$Value) , col= rgb(1 , 0 , 0 , alpha = 0.01))
}
title('P and T waves Distribution')
