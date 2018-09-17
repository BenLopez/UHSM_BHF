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
}

ECGs <- DP_LoadReducedECGs( path , 'z1026' , numberrep = numberrep , FilestoProcess = FilestoProcess )
RPeakData <- DP_LoadRpeaksfile( path , 'z1026' )

QSwidth <- 12
QS <- AFD_ExtractAllSQ(ECG = ECGs$ECGII , RPeaks = RPeakData$RRCombined[10000:11000,] , QSwidth = QSwidth)
goodbeatslogical <- (QS$numvalues < as.numeric(quantile( QS$numvalues  , 0.95)))*(QS$numvalues > as.numeric(quantile( QS$numvalues  , (1-0.95)))) == 1


z <- QS$Value[1,] - median(QS$Value[1,] , na.rm = TRUE)
t <- QS$Date[1,]
plot(t , z  , type ='l')

{x <- c(rep(0,7))
x[1] <- -6  # Baseline
x[2] <- 0.46  # P_Cen  
x[3] <- 20   # P_Amp
x[4] <- 0.1  # P_Width
x[5] <- 0.15  # T_Cen
x[6] <- 10      # T_Amp
x[7] <- 0.111   # T_Width

plot(t , z  , type ='l')
f_x <- ECGSim_QSSimulator(x , t)
Var_md <- ECGSim_QSModelDiscrepancy(X = x , t = t , a = 0 , b = 0.35, d = 0.25 , e = 0.001)
Var_me <- 0.05

lines(t , f_x , col = 'red')
lines(t , f_x + 2*sqrt(Var_md+ Var_me) , col = 'blue')
lines(t , f_x - 2*sqrt(Var_md + Var_me), col = 'blue')
title('Observation, Simulation and Model Discrepancy')
}


