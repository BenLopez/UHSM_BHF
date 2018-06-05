ListGoodPatients <- PatientNames[CorrectNAF ==1]

maxlength <- max(unlist(lapply(RPeaksData , function(X){length(X[[1]]) })))

Moments = list( )

Moments[[ 1 ]] <- lapply(RPeaksData , function(X){  smth( as.numeric(X[[2]])   , method = 'sma' ,  n = 250) })
Moments[[ 2 ]] <- lapply(RPeaksData , function(X){  smth( as.numeric(X[[2]])^2 , method = 'sma' ,  n = 250) })
Moments[[ 3 ]] <- lapply(RPeaksData , function(X){  smth( as.numeric(X[[2]])^3 , method = 'sma' ,  n = 250) })
Moments[[ 4 ]] <- lapply(RPeaksData , function(X){  smth( as.numeric(X[[2]])^4 , method = 'sma' ,  n = 250) })

Moments <- setNames( Moments , c('First' , 'Second' , 'Thrid' , 'Fourth') )

c <- matrix(0 , length(ListGoodPatients) , 50000)
m[ii] <- matrix(0 , length(ListGoodPatients) , 1)
  
for( ii in 1:length(ListGoodPatients) )
{
  
  m[ii , 1] <- mean( RPeaksData[[ListGoodPatients[ii]]]$AFScore$IHAVFScore[(min(501 ,length(RPeaksData[[ListGoodPatients[ii]]]$AFScore$IHAVFScore) ) -1):min(5000 , length(RPeaksData[[ListGoodPatients[ii]]]$AFScore$IHAVFScore))] , rm.na = TRUE )
  DataMatrix[ii , 1:length(RPeaksData[[ListGoodPatients[ii]]]$AFScore$IHAVFScore)] <- 
    RPeaksData[[ListGoodPatients[ii]]]$AFScore$IHAVFScore - m[ii]

}

plot(apply(DataMatrix , 2 , mean) + 3*sqrt(apply(DataMatrix ,2 , function(X){var(X[X>0])}) ) , type = 'l' , ylab = 'E[X] + 3sqrt[X]' , xlab = 't' , xlim = c(0,30000))
abline(70,0)
title('Threshold Choice')