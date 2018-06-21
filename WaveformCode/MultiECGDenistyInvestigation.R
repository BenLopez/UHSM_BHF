{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

DP_LoadPatientIndex()
DP_ChooseDataReps()

PatientID <- DP_choosepatient(listAllPatients)
sub_pat = subset( PatIndex2017, PseudoId %in% PatientID )

outputdata <- DP_LoadRpeaksfile( path , PatientID )
SettingsAFDetection <- AFD_CreateDefaultSettings()
RRStruct <- outputdata$RRCombined
n <- 251
SummaryStats <- list()
SummaryStats[[1]] <- rollmedian( RRStruct$RR , n , na.pad = TRUE)
SummaryStats[[2]] <- rollmean(   RRStruct$RR , n, na.pad = TRUE)
SummaryStats[[3]] <- rollapply(  RRStruct$RR , width = n , FUN = var , na.pad = TRUE)
SummaryStats[[4]] <- rollapply(  RRStruct$RR , width = n , FUN = skewness , na.pad = TRUE)
SummaryStats[[5]] <- rollapply(  RRStruct$RR , width = n , FUN = kurtosis , na.pad = TRUE)
binMatrix <- AFD_CalulateBinMatrixKernelDensityEstimated(RRStruct , n =n)
SummaryStats[[6]] <- apply(binMatrix , 1 , function(X){var(X[X>0] , na.rm = TRUE)} )
SummaryStats[[7]] <- apply(binMatrix , 1 , function(X){max(X , na.rm = TRUE)} )
SummaryStats[[7]][is.infinite(SummaryStats[[7]])] <- mean(SummaryStats[[7]][(2*n):10000] , rm.na = TRUE)
SummaryStats[[8]]  <- apply(binMatrix , 1 , function(X){ length(AFD_ExtractModeStatistics(X)$densities) } )
SummaryStats[[8]][SummaryStats[[8]] == 0] <- 1
SummaryStats[[8]][is.na(SummaryStats[[8]]) == 0] <- 1
SummaryStats[[9]]  <- apply(binMatrix , 1 , function(X){ var(AFD_ExtractModeStatistics(X)$locations) } )
SummaryStats[[9]][is.na(SummaryStats[[9]])] <- mean(SummaryStats[[9]][1:min(10000,length(SummaryStats[[9]]))] , rm.na = TRUE)
SummaryStats[[10]] <- apply(binMatrix , 1 , function(X){ var(AFD_ExtractModeStatistics(X)$densities) } )
SummaryStats[[10]][is.na(SummaryStats[[10]])] <- mean(SummaryStats[[10]][1:min(10000,length(SummaryStats[[10]]))] , rm.na = TRUE)

SummaryStats <- setNames(SummaryStats , c('Median' , 
                                          'Mean' , 
                                          'Variance' , 
                                          'Skewness' , 
                                          'Kurtosis' , 
                                          'Variance densities' ,
                                          'Max Densitiy' , 
                                          'Num Modes',
                                          'Var Modes',
                                          'var Mode Density'))

Beta <- lapply( SummaryStats , function(X){mean(X[1:10000] , na.rm = T)} )

x11(20,20)
par(mfrow = c(5,2))
plot(  outputdata$RRCombined$t , outputdata$RRCombined$RR , xlab = 't', ann=FALSE , col = rgb(0,0,1,alpha = 0.1) )
plot(  outputdata$RRCombined$t , SummaryStats[['Median']] - Beta[[1]] , xlab = 't' , type ='l', ann=FALSE )
lines( outputdata$RRCombined$t , SummaryStats[['Mean']] - Beta[[2]], col = 'red' )
title('Average')
plot(  outputdata$RRCombined$t , SummaryStats[['Variance']]- Beta[[3]] , xlab = 't' , type ='l', ann=FALSE )
title('Variance')
plot(  outputdata$RRCombined$t , SummaryStats[['Skewness']]- Beta[[4]] , xlab = 't' , type ='l', ann=FALSE )
title('Skewness')
plot(  outputdata$RRCombined$t , SummaryStats[['Kurtosis']]- Beta[[5]] , xlab = 't' , type ='l', ann=FALSE )
title('Kurtosis')
plot(  outputdata$RRCombined$t , SummaryStats[['Variance densities']]- Beta[[6]] , xlab = 't' , type ='l', ann=FALSE )
title('Variance densities')
plot(  outputdata$RRCombined$t , SummaryStats[['Max Densitiy']]- Beta[[7]] , xlab = 't' , type ='l', ann=FALSE )
title('Max Densitiy')
plot(  outputdata$RRCombined$t , SummaryStats[[ 'Num Modes']]- Beta[[8]] , xlab = 't' , type ='l', ann=FALSE )
title('Num Modes')
plot(  outputdata$RRCombined$t , SummaryStats[['Var Modes']]- Beta[[9]] , xlab = 't' , type ='l', ann=FALSE )
title('Var Modes')
plot(  outputdata$RRCombined$t , SummaryStats[['var Mode Density']]- Beta[[10]] , xlab = 't' , type ='l', ann=FALSE )
title('Var Mode Density')
x11(20,20)
plot(  outputdata$RRCombined$t , outputdata$RRCombined$RR , xlab = 't', ann=FALSE , col = rgb(0,0,1,alpha = 0.1) )




