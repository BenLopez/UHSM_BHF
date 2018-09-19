# Script to test probabilitic calibration methodod

NumberofSample = 100000
alpha = c(0.2 , 0.8)

for( j in c(2:11) ){
Indexes <- c(1:j)
Sample1 <- mvrnorm(n = round(alpha[1]*NumberofSample) , 
                   mu = t(as.matrix(LocalSecondOrderStruct$mu[1,Indexes])) ,
                   Sigma = LocalSecondOrderStruct$Sigma[1,Indexes,Indexes] )
Sample2 <- mvrnorm(n = round(alpha[2]*NumberofSample) , 
                   mu = t(as.matrix(LocalSecondOrderStruct$mu[2,Indexes])) , 
                   Sigma = LocalSecondOrderStruct$Sigma[2,Indexes,Indexes] )

Z = rbind(Sample1 , Sample2)

f_i <- matrix(0, dim( Z )[ 1 ] , 2 )
for(i in 1:dim(f_i)[2]){
  f_i[ , i] <-  exp(mvnpdf( Z ,  mu = LocalSecondOrderStruct$mu[i,Indexes] , Sigma = LocalSecondOrderStruct$Sigma[i,Indexes,Indexes]))
}

LocalProbCalibrationStruct <- matrix(0, dim( Z )[ 1 ] , 2 )
LocalProbCalibrationStruct[,1] <- (alpha[1]* f_i[ , 1])/ (alpha[1]* f_i[ , 1] + alpha[2]* f_i[ , 2])
LocalProbCalibrationStruct[1:round(alpha[1]*NumberofSample),2] <- 1 

  
ProbabiliticCalibrationOutput <- BC_CreateCalibrationStructure(DP_RemoveNaRows(LocalProbCalibrationStruct) , BinWidth = 0.05)

if(i == j){
  plot(ProbabiliticCalibrationOutput$x , ( ProbabiliticCalibrationOutput$y)  , type = 'l', col = rgb(1/(length(Indexes)-1) , 0 , 0 , alpha = 0.9)  , 
       xlab = 'Estimated Probabilities' , 
       ylab = 'Predicted Probabilities' ,
       ylim = c(0 , 1))
  title('Validation of Validation Methodology')}else{
    lines(ProbabiliticCalibrationOutput$x , ( ProbabiliticCalibrationOutput$y)  , type = 'l', col = rgb(1/(length(Indexes) -1) , 0 , 0 , alpha = 0.9))
  }
}
#
#BC_PlotCreateProbabilityCalibrationPlot(ProbabiliticCalibrationOutput) + ggtitle('Local Calibration Probabilities')

# Estimation Error
for( j in c(2:11) ){
  Indexes <- c(1:j)
  Sample <- list()
  Sample[[1]] <- mvrnorm(n = round(alpha[1]*NumberofSample) , 
                     mu = t(as.matrix(LocalSecondOrderStruct$mu[1,Indexes])) ,
                     Sigma = LocalSecondOrderStruct$Sigma[1,Indexes,Indexes] )
  Sample[[2]] <- mvrnorm(n = round(alpha[2]*NumberofSample) , 
                     mu = t(as.matrix(LocalSecondOrderStruct$mu[2,Indexes])) , 
                     Sigma = LocalSecondOrderStruct$Sigma[2,Indexes,Indexes] )
  
  Z = rbind(Sample[[1]] , Sample[[2]])
  
  f_i <- matrix(0, dim( Z )[ 1 ] , 2 )
  for(i in 1:dim(f_i)[2]){
    f_i[ , i] <-  exp(mvnpdf( Z ,  mu = apply(Sample[[i]] , 2 , mean) , Sigma = cov(Sample[[i]]))) 
  }
  
  LocalProbCalibrationStruct <- matrix(0, dim( Z )[ 1 ] , 2 )
  LocalProbCalibrationStruct[,1] <- (alpha[1]* f_i[ , 1])/ (alpha[1]* f_i[ , 1] + alpha[2]* f_i[ , 2])
  LocalProbCalibrationStruct[1:round(alpha[1]*NumberofSample),2] <- 1 
  
  
  ProbabiliticCalibrationOutput <- BC_CreateCalibrationStructure(DP_RemoveNaRows(LocalProbCalibrationStruct) , BinWidth = 0.01)
  
  if(i == j){
    plot(ProbabiliticCalibrationOutput$x , ( ProbabiliticCalibrationOutput$y)  , type = 'l', col = rgb(1/(length(Indexes)-1) , 0 , 0 , alpha = 0.9)  , 
         xlab = 'Estimated Probabilities' , 
         ylab = 'Predicted Probabilities' ,
         ylim = c(0 , 1))
    title('Validation of Validation Methodology')}else{
      lines(ProbabiliticCalibrationOutput$x , ( ProbabiliticCalibrationOutput$y)  , type = 'l', col = rgb(1/(length(Indexes) -1) , 0 , 0 , alpha = 0.9))
    }
}

Sample <- list()
Sample[[1]] <- t(apply(rmvt(df = 3 , n = round(alpha[1]*NumberofSample) , sigma = LocalSecondOrderStruct$Sigma[1,,])/sqrt(3) , 1 , function(x){x + t(as.matrix(LocalSecondOrderStruct$mu[1,]))}))
Sample[[2]] <- t(apply(rmvt(df = 3 , n = round(alpha[2]*NumberofSample) , sigma = LocalSecondOrderStruct$Sigma[2,,])/sqrt(3) , 1 , function(x){x + t(as.matrix(LocalSecondOrderStruct$mu[2,]))}))


for( j in c(2:11) ){
  Indexes <- c(1:j)

  Z = rbind(Sample[[1]][, Indexes] , Sample[[2]][ , Indexes])
  
  f_i <- matrix(0, dim( Z )[ 1 ] , 2 )
  for(i in 1:dim(f_i)[2]){
    f_i[ , i] <-  exp(mvnpdf( Z ,  mu = apply(Sample[[i]][ , Indexes] , 2 , mean) , Sigma = cov(Sample[[i]][ , Indexes]))) 
  }
  
  LocalProbCalibrationStruct <- matrix(0, dim( Z )[ 1 ] , 2 )
  LocalProbCalibrationStruct[,1] <- (alpha[1]* f_i[ , 1])/ (alpha[1]* f_i[ , 1] + alpha[2]* f_i[ , 2])
  LocalProbCalibrationStruct[1:round(alpha[1]*NumberofSample),2] <- 1 
  
  
  ProbabiliticCalibrationOutput <- BC_CreateCalibrationStructure(DP_RemoveNaRows(LocalProbCalibrationStruct) , BinWidth = 0.05)
  
  if(i == j){
    plot(ProbabiliticCalibrationOutput$x , ( ProbabiliticCalibrationOutput$y)  , type = 'l', col = rgb(1/(length(Indexes)-1) , 0 , 0 , alpha = 0.9)  , 
         xlab = 'Estimated Probabilities' , 
         ylab = 'Predicted Probabilities' ,
         ylim = c(0 , 1))
    title('Validation of Validation Methodology')}else{
      lines(ProbabiliticCalibrationOutput$x , ( ProbabiliticCalibrationOutput$y)  , type = 'l', col = rgb(1/(length(Indexes) -1) , 0 , 0 , alpha = 0.9))
    }
abline(0,1  , col ='blue')
  }


# GMM


