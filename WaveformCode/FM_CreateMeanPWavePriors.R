
{
PsimulatorFunction <- function( x , t_observation ){
  
  Pcen <- x[1]
  Pwidth <- x[2]
  Pcen2 <- x[3]
  Pwidth2 <- x[4]
  output <-(ECGSim_Gaussian( x = t_observation , mu = Pcen , sigma = Pwidth ) + ECGSim_Gaussian( x = t_observation , mu = Pcen2 , sigma = Pwidth2 ))
  
  return( x[5]*(output /max(output)) )
}
MDCOmpnent1 <- function(x , t_observation){
  mu <- 0.5*x[1] + 0.5*x[3]
  output <- 1/ECGSim_Gaussian(Xstar , mu , max(abs(1-mu) ,abs(0.5-mu) ) /3 )
  return(output / max(output))
}
ModelDiscrepancy <- function(x , t_observation, PsimulatorFunction){
  if(x[5] == 0){
  return( rep((1.5)/( 0.5 ) + 0.01 , length(t_observation) ) )  
  }else{
  return( (1.5)/( ((1/x[5])*PsimulatorFunction( x , t_observation )) +  0.5 ) + 5*MDCOmpnent1(x , t_observation)  + 0.01)
  }
}
CalculateImplausability <- function( t_observation , x ,  z  , mdfun = ModelDiscrepancy){
  H = PWaveHM_CreateDesignMatrix(t_observation , x , PsimulatorFunction)
  Beta = PWaveHM_CalculateBetas(H , z)
  Im = mean( abs(z - H%*%Beta) / sqrt(ModelDiscrepancy(x , t_observation, PsimulatorFunction)))
  return( Im )
} 
PWaveHM_CreateDesignMatrix <- function(t_observation , x , PsimulatorFunction){
  t_observation <- as.matrix(t_observation)
  return( cbind( matrix(1 , dim(t_observation)[1] , 1 ) , t_observation , t_observation^2)  )
}
{
{E_Beta = matrix(c(0 , 0 , 0 ) , 3 , 1)
  V_Beta = diag(3)
  V_Beta[1 , 1] <- 1000
  V_Beta[2 , 2] <- 1000
  V_Beta[3 , 3] <- 1000

  C = diag(3)
  C[1 , 2] <- -0.95
  C[2 , 1] <- C[1 , 2]
  C[1 , 3] <- 0.94
  C[3 , 1] <- C[1 , 3]
  C[3 , 2] <- -0.8
  C[2 , 3] <- C[3 , 2] 

  V_Beta <- ((diag(sqrt(V_Beta) ))%*%t(diag(sqrt(V_Beta) )) )*C}
}

EmulatorParameters <- PWaveHM_CreateDefaultEmulationclass()

{ 
  Xstar = seq(0.5 ,1 , 0.01)
  PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.001 , 0.95, 0.001 , 40 ) , b = c(0.55  , 0.05 , 0.55 , 0.05 , -40 ) , numbersamples = 300000 )
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.1 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.02 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]
  PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,5]) < 8 ,  1:4] <- cbind(rep(PriorNonImplausibleSet[which(abs(PriorNonImplausibleSet[,5]) > 8)[1],1] , sum(abs(PriorNonImplausibleSet[,5]) < 8)) , 
                                                                               rep(PriorNonImplausibleSet[which(abs(PriorNonImplausibleSet[,5]) > 8)[1],2] , sum(abs(PriorNonImplausibleSet[,5]) < 8)) , 
                                                                               rep(PriorNonImplausibleSet[which(abs(PriorNonImplausibleSet[,5]) > 8)[1],3] , sum(abs(PriorNonImplausibleSet[,5]) < 8)) , 
                                                                               rep(PriorNonImplausibleSet[which(abs(PriorNonImplausibleSet[,5]) > 8)[1],4] , sum(abs(PriorNonImplausibleSet[,5]) < 8)))
  PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,5]) < 8 ,  5] <- 0
  QSwidth = 10  
  
  BC_PlotPairsFromTwoVariables(X =PriorNonImplausibleSet[which(PriorNonImplausibleSet[,5]==0)[1:2000],] , Y = PriorNonImplausibleSet[which(PriorNonImplausibleSet[,5]!=0)[1:2000],] , alpha = 0.1 , labels = c('X1' ,'X2'  , 'X3','X4','X5' ) , main = TeX('Prior Non-implausible Sets P-waves'))

  {
    Hinvstruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      H = PWaveHM_CreateDesignMatrix(Xstar , X , PsimulatorFunction)
      return(t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + 10*diag(dim(H)[1])  ) )
      #return(solve(t(H)%*%H)%*%t(H))
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
    
    FStruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  PsimulatorFunction(X , Xstar) )
    })
    
    }
}
}

