
# A structure for options of what type to process. 
# BCOptions[['DataType']] = 'DistributionSummaries' , 'CDFs'
# BCOptions[['Classifier Type']] = 'AFClassifier' , 'AFPreAFClassifier'
# BCOptions[[ 'Density Estimation Global']] = 'MVN' , 'GMM'
# BCOptions[['Density Estimation local']] = 'MVN' , 'GMM'


{BCOptions <- setNames( list(1 , 1 , 1 , 1 , 1 , 1) , c('DataType' ,
                                                        'ClassifierType' , 
                                                        'DensityEstimationGlobal', 
                                                        'DensityEstimationlocal' , 
                                                        'GlobalUpdate' ,
                                                        'ProbabilisticCalibration') )
BCOptions[[1]] <- 'DistributionSummaries'
BCOptions[[2]] <- 'AFClassifier' 
BCOptions[[3]] <- 'GMM' 
BCOptions[[4]] <- 'GMM'
BCOptions[[5]] <- 'Yes'
BCOptions[[6]] <- 'No'


# Parameters for Priors
BCParameters <- setNames( list(1,1,1,1 , 1) , c( 'TS_Likelihood_clique' , 'NumberComponentsforGMM' , 'ProbabilityThreshold' , 'Minute Threshold' ,'GlobalNumberComponentsforGMM' ) )
BCParameters[[1]] <- 100
BCParameters[[2]] <- 25
BCParameters[[3]] <- 0.8
BCParameters[[4]] <- 6
BCParameters[[5]] <- 1

# Prior probabilities of going into AF.
Priorprobabilities <- setNames( list(1 , 1 , 1 , 1 , 1 , 1) ,
                                c('A' ,
                                  'A^c',
                                  'B|A',
                                  'B|A^c',
                                  'B',
                                  'B^c') )
Priorprobabilities[[ 1 ]] <- 0.2
Priorprobabilities[[ 2 ]] <- 1 - Priorprobabilities[[1]]
Priorprobabilities[[ 3 ]] <- 0.0005
Priorprobabilities[[ 4 ]] <- 0
Priorprobabilities <- BC_CaculateBeatWiseProbabiltiesAFDetection(Priorprobabilities)

precomputedfolderpath <- DP_SelectPrecomputedFolder()
}