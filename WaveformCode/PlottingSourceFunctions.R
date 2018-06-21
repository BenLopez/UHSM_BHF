

PlD_DistributionSummaryStats <- function(outputdata , SummaryStats){
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
}
  