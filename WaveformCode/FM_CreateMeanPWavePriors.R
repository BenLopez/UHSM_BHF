{
PsimulatorFunction <- function( x , t_observation ){
  
  Pcen <- x[1]
  Pwidth <- x[2]
  Pcen2 <- x[3]
  Pwidth2 <- x[4]
  output1 <- (ECGSim_Gaussian( x = t_observation , mu = Pcen , sigma = Pwidth ))
  output2 <- (ECGSim_Gaussian( x = t_observation , mu = Pcen2 , sigma = Pwidth2 ))
  return( x[5]*(output1 /max(ECGSim_Gaussian( x = seq(0.5,1,0.5/510) , mu = Pcen , sigma = Pwidth ))) + x[6]*(output1 / max(ECGSim_Gaussian( x = seq(0.5,1,0.5/510) , mu = Pcen2 , sigma = Pwidth2 ))))  
}
MDCOmpnent1 <- function(x , t_observation){
  mu <- 0.5*x[1] + 0.5*x[3]
  output <- 1/ECGSim_Gaussian(Xstar , mu , max(abs(1-mu) ,abs(0.5-mu) ) /3 )
  return(output / max(output))
}
ModelDiscrepancy <- function(x , t_observation, PsimulatorFunction , eta = 0.01 , alpha = 0.4){
  if(x[5] == 0){
    return( rep(1, length(t_observation) ) )  
  }else{
    ff <- PsimulatorFunction(x , t_observation)
    return( (alpha*abs(ff)) +  ((eta*max(abs(ff)))/(abs(ff) + 1))   )
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
  
{ 
  Xstar = seq(0.5 ,1 , 0.5/51)
  PriorNonImplausibleSet <- BE_SampleLHSinab( a = c( 0.95 , 0.01 , 0.95, 0.01 , 20 , 20 ) , b = c(0.55  , 0.05 , 0.55 , 0.05 , -20 , -20 ) , numbersamples = 5000000 )
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,1] - PriorNonImplausibleSet[,3]) < 0.2 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(PriorNonImplausibleSet[,2] - PriorNonImplausibleSet[,4]) < 0.03 ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ PriorNonImplausibleSet[,1] < PriorNonImplausibleSet[,3] ,  ]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ (abs(PriorNonImplausibleSet[,5]) + abs(PriorNonImplausibleSet[,6])) > 8 ,  1:6]
  PriorNonImplausibleSet <- PriorNonImplausibleSet[ abs(abs(PriorNonImplausibleSet[,5]) - abs(PriorNonImplausibleSet[,6])) < 10  ,  1:6]
  
  
    FStruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  PsimulatorFunction(X , Xstar) )
    })
    PriorNonImplausibleSet <- PriorNonImplausibleSet[apply(abs(FStruct) , 2 , max) > 8 & apply(abs(FStruct) , 2 , max) < 40 , ]

    PriorNonImplausibleSet <- rbind(PriorNonImplausibleSet,t(matrix(c(0.95 , 0.001 , 0.95, 0.001 , 0 , 0) , 6 , 2)) )
    FStruct <- apply(PriorNonImplausibleSet , 1 , function(X){
      return(  PsimulatorFunction(X , Xstar) )
    })
    PriorNonImplausibleSet <- PriorNonImplausibleSet[DP_FindNARows(t(FStruct) ),]
    FStruct <- FStruct[,DP_FindNARows(t(FStruct) )]
    FStruct <- t(FStruct)

    XPwave <- t(apply(FStruct[1:(dim(FStruct)[1] - 2) , ] ,1,  function(X){as.matrix(c(Xstar[which(abs(X)>0.01)[1]],Xstar[which(abs(X)>0.01)[length(which(abs(X)>0.01))]]) ) } ))
    XPwave <- rbind(XPwave , c(0.5,1),c(0.5,1))
    XPwave <- t(apply(XPwave , 1 , function(X){seq(from = max(X[1]-0.05,0.5) ,to =  min(X[2]+0.05,1) ,by =  (abs(max(X[1]-0.05,0.5)  - min(X[2]+0.05,1))/51) )[1:51] }))
  
    FStruct <- apply(cbind(PriorNonImplausibleSet ,XPwave ) , 1 , function(X){
      return(  PsimulatorFunction(X[1:6] , X[7:length(X)]) )
    })
}
    {  
    Hinvstruct <- apply(XPwave , 1 , function(X){
      H = PWaveHM_CreateDesignMatrix(X , PriorNonImplausibleSet[1,] , PsimulatorFunction)
      return(t(H%*%V_Beta)%*%solve(H%*%V_Beta%*%t(H) + 10*diag(dim(H)[1])  ) )
      #return(solve(t(H)%*%H)%*%t(H))
    })
    Hstruct <- apply(XPwave , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(X , PriorNonImplausibleSet[1,] , PsimulatorFunction))
    })
    EHstruct <- apply(XPwave , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(X , PriorNonImplausibleSet[1,] , PsimulatorFunction)%*%E_Beta)
    })
    HstructP <- apply(XPwave , 1 , function(X){
      return(  H = PWaveHM_CreateDesignMatrix(X , PriorNonImplausibleSet[1,] , PsimulatorFunction)[ , 1:3] )
    })
    
    ModelDiscrepancyMatrix <- apply(cbind(PriorNonImplausibleSet ,XPwave ) , 1 , function(X){(ModelDiscrepancy(X[1:6] , X[7:length(X)], PsimulatorFunction))} )
    }
}
}
XPwave[is.na(XPwave)] <- 0.5

source('CTEM_BuildPwaveEmulators.R')



x11(20 , 14)
pairs( PriorNonImplausibleSet[1:2000,] , col = rgb(0 , 0 , 1 , alpha = 0.05) , pch = 16 ) 