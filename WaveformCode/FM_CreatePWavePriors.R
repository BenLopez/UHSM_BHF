
PsimulatorFunction <- function( x , t_observation ){
  
  Pcen <- x[1]
  Pwidth <- x[2]
  Pcen2 <- x[3]
  Pwidth2 <- x[4]
  output <- ECGSim_Gaussian( x = t_observation , mu = Pcen , sigma = Pwidth ) + ECGSim_Gaussian( x = t_observation , mu = Pcen2 , sigma = Pwidth2 )
  
  return( output /max(output) )
}
ModelDiscrepancy <- function(x , t_observation, PsimulatorFunction){
  1/( PsimulatorFunction( x , t_observation ) + 0.5 )
}
CalculateImplausability <- function( t_observation , x ,  z  , mdfun = ModelDiscrepancy){
  H = PWaveHM_CreateDesignMatrix(t_observation , x , PsimulatorFunction)
  Beta = PWaveHM_CalculateBetas(H , z)
  Im = mean( abs(z - H%*%Beta) / sqrt(ModelDiscrepancy(x , t_observation, PsimulatorFunction)))
  return( Im )
} 
PWaveHM_CreateDesignMatrix <- function(t_observation , x , PsimulatorFunction){
  t_observation <- as.matrix(t_observation)
  return( cbind( matrix(1 , dim(t_observation)[1] , 1 ) , t_observation , t_observation^2, as.matrix(PsimulatorFunction(x , t_observation))  ) )
}

EmulatorParameters <- PWaveHM_CreateDefaultEmulationclass()

{ Xstar = seq(0.5 ,1 , 0.01)
  PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.001 , 0.95, 0.001 ) , b = c(0.6  , 0.04 , 0.6 , 0.04 ) , numbersamples = 100000 )
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.1 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.02 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]
  QSwidth = 10 
  
  {Hinvstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
    H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)
    return(solve(t(H)%*%H)%*%t(H))
  })
    Hstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction))
    })
    ModelDiscrepancyMatrix <- apply(PriorNonImplausibleSet , 1 , function(X){sqrt(ModelDiscrepancy(X , Xstar , PsimulatorFunction))} )
  }
  
}