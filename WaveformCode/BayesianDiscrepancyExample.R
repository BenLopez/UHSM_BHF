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
  set.seed(1)
}

beta = c(0.3 , 0.7)

{
MclustDistributionStruct <- BC_CreateDefaultmclustStruct()

MclustDistributionStruct$parameters$pro <- c(0.1 , 0.8 ,0.1)
MclustDistributionStruct$parameters$mean <- matrix(0, 1 , 3)
MclustDistributionStruct$parameters$mean[ , 1] <- 5
MclustDistributionStruct$parameters$mean[ , 2] <- 10
MclustDistributionStruct$parameters$mean[ , 3] <- 15

MclustDistributionStruct$parameters$variance$sigma <- array(0 , c(1 , 1 , 3))
MclustDistributionStruct$parameters$variance$sigma[ , , 1] <- c( 2.5) 
MclustDistributionStruct$parameters$variance$sigma[ , , 2] <- c( 5) 
MclustDistributionStruct$parameters$variance$sigma[ , , 3] <- c( 2.5)
x11()
SampleofPoints <- BC_SampleGMM(MclustDistributionStruct = MclustDistributionStruct , numberofsamples =    beta[1]*100000)
hist(SampleofPoints , breaks = 50 , col = rgb(1 , 0 , 0 , alpha = 0.5), freq = F , xlab = TeX('z') , ylab = TeX('Density') , main = TeX('Histograms of Datasets for Both Classes'))

MclustDistributionStruct2 <- BC_CreateDefaultmclustStruct()

MclustDistributionStruct2$parameters$pro <- c(0.1 , 0.8 ,0.1)
MclustDistributionStruct2$parameters$mean <- matrix(0, 1 , 3)
MclustDistributionStruct2$parameters$mean[ , 1] <- 0
MclustDistributionStruct2$parameters$mean[ , 2] <- 5
MclustDistributionStruct2$parameters$mean[ , 3] <- 10

MclustDistributionStruct2$parameters$variance$sigma <- array(0 , c(1 , 1 , 3))
MclustDistributionStruct2$parameters$variance$sigma[ , , 1] <- c( 2.5) 
MclustDistributionStruct2$parameters$variance$sigma[ , , 2] <- c( 5) 
MclustDistributionStruct2$parameters$variance$sigma[ , , 3] <- c( 2.5)

SampleofPoints2 <- BC_SampleGMM(MclustDistributionStruct = MclustDistributionStruct2 , numberofsamples =  beta[2]*100000)

hist(SampleofPoints2 , breaks = 50 , col = rgb(0 , 0 , 1 , alpha = 0.5) , add = TRUE , freq = F)
}



m1 <- mean(SampleofPoints)
s1 <- var(SampleofPoints)
m2 <- mean(SampleofPoints2)
s2 <- var(SampleofPoints2)

Z <- rbind(SampleofPoints , SampleofPoints2)
P = (beta[1]*dnorm(Z , m1 , sqrt(s1))) / (beta[1]*dnorm(Z , m1 , sqrt(s1)) + beta[2]*dnorm(Z , m2 , sqrt(s2)) )


ProbabiliticCalibrationOutput <- BC_CreateCalibrationStructure(BC_CreateProbCalibrationStruct(P  , alpha = beta,numberofvalidationsamples = 100000 ) , BinWidth = 0.05)
x11()
Data <- ProbabiliticCalibrationOutput
GlobalProbabilityCalibrationPlot <- ggplot(Data  , aes(x = x , y = y)) +
  geom_point( color = 'blue') +
  geom_errorbar(aes(ymin =  y - 2*sd , ymax = y + 2*sd ) , width = .01 ) +
  geom_line(aes(x = x , y = x))+
  xlab(TeX('$B( P , z )$')) +
  ylab(TeX('Estimated $P_t(X| z)$' )) +
  ggtitle(TeX('Bayesian Estimated probabilities Against Data Estimated Probabilities'))
print(GlobalProbabilityCalibrationPlot)


{Z <- rbind(as.matrix(rnorm(beta[1]*100000 , m1 , sqrt(s1))) , as.matrix(rnorm(beta[2]*100000 , m2 , sqrt(s2))))
P = (beta[1]*dnorm(Z , m1 , sqrt(s1))) / (beta[1]*dnorm(Z , m1 , sqrt(s1)) + beta[2]*dnorm(Z , m2 , sqrt(s2)) )


ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(BC_CreateProbCalibrationStruct(P  , alpha = beta,numberofvalidationsamples = 100000 ) , BinWidth = 0.05))
x11()
Data <- ProbabiliticCalibrationOutput
GlobalProbabilityCalibrationPlot <- ggplot(Data  , aes(x = x , y = y)) +
  geom_point( color = 'blue') +
  geom_errorbar(aes(ymin =  y - 2*sd , ymax = y + 2*sd ) , width = .01 ) +
  geom_line(aes(x = x , y = x))+
  xlab(TeX('$B( P , z )$')) +
  ylab(TeX('Estimated $P_t(X| z)$' )) +
  ggtitle(TeX('Analysis Given Known Densities'))
print(GlobalProbabilityCalibrationPlot)}



# Create Sen against Spec curves

prob_thresh = seq(0 , 1 , 0.01)
ProbCalibStruct <- BC_CreateProbCalibrationStruct(P  , alpha = beta,numberofvalidationsamples = 100000 )
ROCStruct <- matrix(0 , length(prob_thresh) , 4)

for(i in 1:dim(ROCStruct)[1]){
  temp <- BD_CalulateSenSpecNPVPPV(ProbCalibStruct , prob_thresh[i])
  ROCStruct[i,1] <-temp$Sen
  ROCStruct[i,2] <- temp$Spec
  ROCStruct[i,3] <- temp$NPV
  ROCStruct[i,4] <- temp$PPV
}


x11()
p1 <- ggplot( data.frame(Sensitivity = ROCStruct[ , 1] , Specificity = ROCStruct[ , 2]) , aes(Sensitivity , Specificity)) + 
  geom_point( color = 'blue') + ggtitle('ROC Curve')
p2 <- ggplot( data.frame(PPV = ROCStruct[ , 3] , NPV = ROCStruct[ , 4]) , aes(PPV , NPV)) + 
  geom_point( color = 'blue') + ggtitle('NPV vs PPV Curve')
grid.arrange(p1 , p2 , nrow = 1 , ncol = 2)

# Create a Bayes Linear Update using Data from emperical calibration curve.



plot(ProbabiliticCalibrationOutput$y , ProbabiliticCalibrationOutput$x )
abline(0,1)
