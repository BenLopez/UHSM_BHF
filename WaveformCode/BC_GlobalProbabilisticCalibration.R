
DataSetPriorProbabilities <- BC_ExtractValidationPriors(Priorprobabilities , DataBaseMaster , DataBase)

GlobalUpdateDiagnostics <- matrix(0 , dim(DataBase[[3]])[1] + dim(DataBase[[4]])[1] , 1 )
counter <-1

if(BCOptions$DensityEstimationGlobal == 'MVN'){
  for(i in 1:dim(DataBase[[3]])[1]){
    GlobalUpdateDiagnostics[counter,] <- BC_GlobalBayesianBeliefUpdateMVN(M = DataBase[[3]][i,] , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = Priorprobabilities)$A
    counter <- counter +1
  }
  for(i in 1:dim(DataBase[[4]])[1]){
    GlobalUpdateDiagnostics[counter,] <- BC_GlobalBayesianBeliefUpdateMVN(M = DataBase[[4]][i,] , GlobalSecondOrderStruct = GlobalSecondOrderStruct , Priorprobabilities = Priorprobabilities)$A
    counter <- counter +1
  }
}
if(BCOptions$DensityEstimationGlobal == 'GMM'){
  for(i in 1:dim(DataBase[[3]])[1]){
    if(sum(is.na(t(as.matrix(DataBase[[3]][i,]))))>0){next}
    GlobalUpdateDiagnostics[counter,] <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix(DataBase[[3]][i,])) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = Priorprobabilities)$A
    counter <- counter +1
  }
  for(i in 1:dim(DataBase[[4]])[1]){
    if(sum(is.na(t(as.matrix(DataBase[[4]][i,]))))>0){next}
    GlobalUpdateDiagnostics[counter,] <- BC_GlocalBayesianBeliefUpdateGMM(M = t(as.matrix(DataBase[[4]][i,])) , GlobalDistributionStruct = GlobalDistributionStruct , Priorprobabilities = Priorprobabilities)$A
    counter <- counter +1
  }
}  

GlobalProbCalibrationStruct <- cbind(GlobalUpdateDiagnostics , 0*GlobalUpdateDiagnostics)
GlobalProbCalibrationStruct[1:52 , 2] <- 1 
GlobalProbCalibrationStruct <- DP_RemoveNaRows(GlobalProbCalibrationStruct)

#plot(GlobalProbCalibrationStruct[,1] , GlobalProbCalibrationStruct[,2]  , pch = 16 , col = rgb(0,0,1,alpha = 0.1))

BinWidth <- 0.1
ProbabililtyBins <- seq(0,1,BinWidth)
ProbabililtyBinStruct <-  matrix(0 , (length(ProbabililtyBins) -1) , 1)
EstimatorError <- ProbabililtyBinStruct
  
for(i in 1:(dim(ProbabililtyBinStruct)[1])){
  tmplogical <- ((GlobalProbCalibrationStruct[,1] >= ProbabililtyBins[i])*(GlobalProbCalibrationStruct[,1] <= ProbabililtyBins[i+1])) == 1
  ProbabililtyBinStruct[i,] <- sum(GlobalProbCalibrationStruct[tmplogical , 2])/length((GlobalProbCalibrationStruct[tmplogical , 2]))
  EstimatorError[i,] <- (ProbabililtyBinStruct[i,])*(1-ProbabililtyBinStruct[i,])/length((GlobalProbCalibrationStruct[tmplogical , 2]))
}
x11()
EstimatorError[is.na(EstimatorError)] <- max(EstimatorError[!is.na(EstimatorError)] )
EstimatorError[EstimatorError < 0.0000000000001] <- max(EstimatorError[!is.na(EstimatorError)] )
GlobalProbabilityCalibrationPlot <- ggplot(data.frame(x =  ProbabililtyBins[1:dim(ProbabililtyBinStruct)[1]] + BinWidth/2, 
                                                     y = ProbabililtyBinStruct,  
                                                     sd = sqrt(EstimatorError) )  , aes(x = x , y = y)) +
                                    geom_point( color = 'blue') +
                                    geom_errorbar(aes(ymin = y - 2*sd , ymax = y + 2*sd ) , width = .01 ) +
                                    geom_line(aes(x = x , y = x))+
                                    ggtitle('Global Calibrated Expert Probability') +
                                    xlab('Predicted Probability') +
                                    ylab('Estimated Probability')
print(GlobalProbabilityCalibrationPlot)

#bound = 0.000001

#sensitivity <- sum( GlobalUpdateDiagnostics[c(1:52), ] > bound )/52
#specificity <- sum(GlobalUpdateDiagnostics[-c(1:52), ] < bound)/250 
#PPV <- ( 0.2*sensitivity )/( (0.2*sensitivity) + ( (1-specificity)*0.8 ))
#NPV <- ( 0.8*specificity )/( (0.2*(1-sensitivity)) + ( (specificity)*0.8 ))
