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

numbertrainingpoints = 100000

# Set up data parameters
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
SampleofPoints <- BC_SampleGMM(MclustDistributionStruct = MclustDistributionStruct , numberofsamples =    beta[1]*numbertrainingpoints)
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

SampleofPoints2 <- BC_SampleGMM(MclustDistributionStruct = MclustDistributionStruct2 , numberofsamples =  beta[2]*numbertrainingpoints)

hist(SampleofPoints2 , breaks = 50 , col = rgb(0 , 0 , 1 , alpha = 0.5) , add = TRUE , freq = F)




m1 <- mean(SampleofPoints)
s1 <- var(SampleofPoints)
m2 <- mean(SampleofPoints2)
s2 <- var(SampleofPoints2)

Z <- rbind(SampleofPoints , SampleofPoints2)
P = (beta[1]*dnorm(Z , m1 , sqrt(s1))) / (beta[1]*dnorm(Z , m1 , sqrt(s1)) + beta[2]*dnorm(Z , m2 , sqrt(s2)) )


ProbabiliticCalibrationOutput1 <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(BC_CreateProbCalibrationStruct(P  , alpha = beta,numberofvalidationsamples = numbertrainingpoints ) ,
                                                                BinWidth = 0.05))
x11()
Data <- ProbabiliticCalibrationOutput1
GlobalProbabilityCalibrationPlot1 <- ggplot(Data  , aes(x = x , y = y)) +
  geom_point( color = 'blue') +
  geom_errorbar(aes(ymin =  y - 2*sd , ymax = y + 2*sd ) , width = .01 ) +
  geom_line(aes(x = x , y = x))+
  xlab(TeX('$B( P , z )$')) +
  ylab(TeX('Estimated $P_t(X| z)$' )) +
  ggtitle(TeX('Bayesian Estimated probabilities Against Data Estimated Probabilities'))
print(GlobalProbabilityCalibrationPlot1)
}

{Z2 <- rbind(as.matrix(rnorm(beta[1]*numbertrainingpoints , m1 , sqrt(s1))) , as.matrix(rnorm(beta[2]*numbertrainingpoints , m2 , sqrt(s2))))
P2 = (beta[1]*dnorm(Z2 , m1 , sqrt(s1))) / (beta[1]*dnorm(Z2 , m1 , sqrt(s1)) + beta[2]*dnorm(Z2 , m2 , sqrt(s2)) )


ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(BC_CreateProbCalibrationStruct(P2  , alpha = beta,numberofvalidationsamples = numbertrainingpoints ) ,
                                                                                             BinWidth = 0.05))
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
{
  prob_thresh = seq(0 , 1 , 0.01)
ProbCalibStruct <- BC_CreateProbCalibrationStruct(P  , alpha = beta,numberofvalidationsamples = numbertrainingpoints )
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
}


# Inverse problem solve 

F1 <- kde( x = SampleofPoints  , h = 0.1*var(SampleofPoints)   , eval.points = Z)$estimate
F2 <- kde( x = SampleofPoints2  ,  h = 0.1*var(SampleofPoints)   , eval.points = Z)$estimate
P2 <- as.matrix(beta[1]*F1 / (beta[1]*F1 + beta[2]*F2))
ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(BC_CreateProbCalibrationStruct(PosteriorProbabilities = P2  , alpha = beta , numberofvalidationsamples = numbertrainingpoints ) ,
                                                                                             BinWidth = 0.05))

Simulator <- function(X , X2 ){
  F1 <- kde( x = SampleofPoints  , h = X*var(SampleofPoints)   , eval.points = Z)$estimate
  F2 <- kde( x = SampleofPoints2  ,  h = X2*var(SampleofPoints)   , eval.points = Z)$estimate
  P2 <- as.matrix(beta[1]*F1 / (beta[1]*F1 + beta[2]*F2))
  ProbabiliticCalibrationOutput <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(BC_CreateProbCalibrationStruct(PosteriorProbabilities = P2  , alpha = beta , numberofvalidationsamples = numbertrainingpoints ) ,
                                                                                               BinWidth = 0.05))
  
  return(ProbabiliticCalibrationOutput)
}

ImMeasure <- function( ProbabiliticCalibrationOutput ){
  return(abs( max((ProbabiliticCalibrationOutput$x - ProbabiliticCalibrationOutput$y)/ProbabiliticCalibrationOutput$sd)) )
}

PriorSample <- randomLHS(n = 100000 , k = 2 )
Implausability <- matrix(0 , dim(PriorSample)[1] , 1)


for(i in 1:dim(PriorSample)[1]){
  Implausability[i, ] <- ImMeasure(BC_CleanProbCalibrationOutput(Simulator(PriorSample[i,1] , PriorSample[i,2])))
  DP_WaitBar(i/dim(PriorSample)[1])
}



plot(PriorSample[Implausability >3 , 1] , PriorSample[Implausability >3 , 2] , pch = 16 , col = rgb(1,0,0, alpha = 0.5))
points(PriorSample[Implausability <3 , 1] , PriorSample[Implausability <3 , 2] , pch = 16 , col = rgb(0,0,1, alpha = 0.5))


HistoryMatchOutput <- BE_HistoryMatch(TrainingSet , TrainingSet2, EmulatorSettings = EmulatorSettings  , HistoryMatchSettings = HistoryMatchSettings , PriorRange = PriorRange )
EmulatorSettings = HistoryMatchOutput$EmulatorSettings
chi_star = HistoryMatchOutput$chi_star


# Create a Bayes Linear Update using Data from emperical calibration curve.
{
EmulatorSettings <- BE_CreateDefaultEmulationClass()
EmulatorSettings$MeanFunction <- function(X){
  X <- as.matrix(X)
  H = cbind(as.matrix(1 + 0*X)  )
  return(H)
}
EmulatorSettings$X <- ProbabiliticCalibrationOutput1$y
EmulatorSettings$Y <- ProbabiliticCalibrationOutput1$x - ProbabiliticCalibrationOutput1$y
EmulatorSettings$w <- function(X){
  return( mean( ProbabiliticCalibrationOutput1$sd^2 )/50 )
}

xstar <- seq(0,1,0.001)

em_output <- BE_BayesLinearEmulatorLSEstimates(xstar =  xstar , EmulatorSettings = EmulatorSettings)
BE_PlotOneDOutput(Em_output = em_output , EmulatorSettings = EmulatorSettings , Xstar = xstar)

xstar <- P
em_output <- BE_BayesLinearEmulatorLSEstimatesBatchMode(xstar =  xstar , EmulatorSettings = EmulatorSettings)



ProbabiliticCalibrationOutput2 <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(BC_CreateProbCalibrationStruct(P - em_output$E_D_fX  ,
                                                                                               alpha = beta,numberofvalidationsamples = numbertrainingpoints ) ,
                                                                BinWidth = 0.05))
Data <- ProbabiliticCalibrationOutput2
GlobalProbabilityCalibrationPlot2 <- ggplot(Data  , aes(x = x , y = y)) +
  geom_point( color = 'blue') +
  geom_errorbar(aes(ymin =  y - 2*sd , ymax = y + 2*sd ) , width = .01 ) +
  geom_line(aes(x = x , y = x))+
  xlab(TeX('Estimated $P_t(X| z)$ from adjusted expectations of P_t(X| z)')) +
  ylab(TeX('Estimated $P_t(X| z)$' )) +
  ggtitle(TeX('Adjusted Expected probabilities Against Data Estimated Probabilities'))

x11()
grid.arrange(GlobalProbabilityCalibrationPlot1 , GlobalProbabilityCalibrationPlot2 , nrow = 1 , ncol = 2)
}
