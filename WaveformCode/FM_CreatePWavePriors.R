
PsimulatorFunction <- function( x , t_observation ){
  
  Pcen <- x[1]
  Pwidth <- x[2]
  Pcen2 <- x[3]
  Pwidth2 <- x[4]
  output <- ECGSim_Gaussian( x = t_observation , mu = Pcen , sigma = Pwidth ) + ECGSim_Gaussian( x = t_observation , mu = Pcen2 , sigma = Pwidth2 )
  
  return( output /max(output) )
}
ModelDiscrepancy <- function(x , t_observation, PsimulatorFunction){
  return(4/( PsimulatorFunction( x , t_observation ) + 0.5 ) + 0.01)
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
{E_Beta = matrix(c(0 , 0 , 0 , 14) , 4 , 1)
  V_Beta = diag(diag(matrix(0 , 4 , 4)))
  V_Beta[1 , 1] <- 1000
  V_Beta[2 , 2] <- 1000
  V_Beta[3 , 3] <- 1000
  V_Beta[4 , 4] <- 450
  
  C = diag(4)
  C[1 , 2] <- -0.95
  C[2 , 1] <- C[1 , 2]
  C[1 , 3] <- 0.94
  C[3 , 1] <- C[1 , 3]
  C[1 , 4] <- 0.2
  C[4 , 1] <- C[1 , 4]
  C[3 , 2] <- -0.8
  C[2 , 3] <- C[3 , 2] 
  C[4 , 2] <- -0.2
  C[2 , 4] <- C[4 , 2]
  C[4 , 3] <- 0.2
  C[3 , 4] <- C[4 , 3]
  
  V_Beta <- ((diag(sqrt(V_Beta) ))%*%t(diag(sqrt(V_Beta) )) )*C}


EmulatorParameters <- PWaveHM_CreateDefaultEmulationclass()

{ 
  Xstar = seq(0.5 ,1 , 0.01)
  PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.001 , 0.95, 0.001 ) , b = c(0.55  , 0.05 , 0.55 , 0.05 ) , numbersamples = 100000 )
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.1 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.02 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]
  QSwidth = 10 
  
  {
    Hinvstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
    H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)
    return(t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + 10*diag(dim(H)[1])  ) )
    #return(solve(t(H)%*%H)%*%t(H))
  })
  HinvstructP <- apply(PriorNonImplausibleSet , 1 , function(X){
      H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)
      return(solve(t(H[,1:3])%*%H[,1:3])%*%t(H[,1:3]))
    })
    Hstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction))
    })
    EHstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)%*%E_Beta)
    })
    
    HstructP <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)[ , 1:3] )
    })
    ModelDiscrepancyMatrix <- apply(PriorNonImplausibleSet , 1 , function(X){(ModelDiscrepancy(X , Xstar , PsimulatorFunction))} )
  }
  
}


