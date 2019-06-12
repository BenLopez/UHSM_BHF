{
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
    set.seed( 1 )
  }
}  

numberofsimulations = 1100
numberofrepeats = 50000
precomputedfolderpath <- DP_SelectPrecomputedFolder()

if( DP_CheckfileinPrecomputedfolder(precomputedfolderpath,'EmulatingCriticalThresholdDataCDF.RData') ){
  load(paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataCDF.RData'))
  Im_Crit <- EmulatingCriticalThresholdData[[1]]
  X <- EmulatingCriticalThresholdData[[2]]
  xx <- seq(0,5,0.01)
}else{
  xx <- seq(0,5,0.01)
  outputstruct <- CTEm_CreateTrainingSetCDF(xx,x,PriorNonImplausibleSetTotal,F_total,f_total,MD_Total  ,numberofsimulations = numberofsimulations, c = c , l = l , N = n ,numberofrepeats = numberofrepeats)
  EmulatingCriticalThresholdData <- list( Im_crit = outputstruct , X = PriorNonImplausibleSetTotal[1:numberofsimulations , ])
  save( EmulatingCriticalThresholdData , file = paste0(precomputedfolderpath , '\\EmulatingCriticalThresholdDataCDF.RData') )
}

x11()
plot(xx ,  Im_Crit[1,,1] , type ='l', col = rgb(0,0,1,alpha = 0.1) , xlim = c(0,3) , ylab = ' CDF' , xlab = 'Mean Implausibility')
for(i in 2:dim(Im_Crit)[1]){
  lines(xx ,  Im_Crit[i,,1] , col = rgb(0,0,1,alpha = 0.01) )
}
title('Simulated Mean Implausibility CDFs')
abline(h = 0.99)
x11()
plot(xx , 1 - Im_Crit[1,,2] , type ='l', col = rgb(0,0,1,alpha = 0.1) , xlim = c(0,5) , ylab = '1 - CDF' , xlab = 'Max Implausibility')
for(i in 2:dim(Im_Crit)[1]){
  lines(xx , 1 - Im_Crit[i,,2] , col = rgb(0,0,1,alpha = 0.01) )
}
title('Emulating Observed Implausibility')
abline(h = 0.01)

##### Build Emulators #####

Training_X = X[1:1000 , ]
Training_ImCrit = Im_Crit[1:1000 , , 1:2]
Validation_X = X[1001:1100 ,]
Validation_ImCrit = Im_Crit[1001:1100 , , 1:2]

{ 
  EmulatorParametersCDFMean <- BE_CreateDefaultEmulationClass()
  EmulatorParametersCDFMean$X <- Training_X[ , -3]
  EmulatorParametersCDFMean$Y <- Training_ImCrit[ , , 1]
  EmulatorParametersCDFMean$w <- function(X){
    return(0.0001*diag(dim(as.matrix(X))[1]))
  }
 # EmulatorParametersCDFMean$MeanFunction <- function(X){
 #   X <- as.matrix(X)
 #   H = cbind(as.matrix(1 + 0*X[,1]) , X ) 
#    return(H)}
  EmulatorParametersCDFMean$MeanFunction <- function(X){
   X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X ,
              X[,1]*X[ , 2:dim(X)[2]] , 
              X[,2]*X[ , 3:dim(X)[2]] , 
              X[,4]*X[ , 5:dim(X)[2]] , 
              X[,6]*X[ , 7:dim(X)[2]] ,
              X[,7]*X[ , 8:dim(X)[2]] ,
              X[,8]*X[ , 9:dim(X)[2]] ,
              X^2 )
    return(H)
  }
  EmulatorParametersCDFMean$CorrelationLength <- function(X , n){
    return(0.25*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  {
    Xstar <- Validation_X[ , -3]
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimatesMO(xstar = Xstar ,EmulatorSettings = EmulatorParametersCDFMean  )
    
    StdResiduals <- matrix(0 , dim(EmulatorOutput$E_D_fX)[1] , dim(EmulatorOutput$E_D_fX)[2])
    
    for(i in 1:dim(Validation_ImCrit)[1]){
      StdResiduals[i,] <-  ( (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ))
    }  
    
    sampleofpoints <- sample(1:dim(StdResiduals)[1] , 100)
    
    x11()
    output <- ggplot( data.frame(x_2 = Validation_X[sampleofpoints , 8] , StdResids = StdResiduals[sampleofpoints , 116]) , aes(x_2 , StdResids) ) +
      geom_point( color = 'blue') + ggtitle('Mean Implausibility Validation') + geom_hline( yintercept  = 3)+ geom_hline( yintercept  = -3) 

    print(output)
    
    #x11()
    #plot(xx,(EmulatorOutput$E_D_fX[1,] -  Validation_ImCrit[ 1 ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[1,1]*EmulatorOutput$SigmaHat) ) , type ='l'  , col = rgb(0,0,1,alpha =0.1) , ylim = c(0,5) , xlab = 'a' , ylab = 'StandardisedResiudals' )
    #for(i in 2:dim(Validation_ImCrit)[1]){
    #  lines(xx , (EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,, 1])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ), col = rgb(0,0,1,alpha =0.1))
    #}
      
    Rsq <- 1 - var( (matrix(EmulatorOutput$E_D_fX , length(Validation_ImCrit[ ,  , 1]) , 1 ) -  matrix(Validation_ImCrit[ ,  , 1] , length(Validation_ImCrit[ ,  , 1]) , 1) ) )/var(matrix(Validation_ImCrit[ ,  , 1] , length(Validation_ImCrit[ ,  , 1]) , 1))
    
  }
}

{ 
  EmulatorParametersCDFMax <- BE_CreateDefaultEmulationClass()
  EmulatorParametersCDFMax$X <- Training_X[ , -3]
  EmulatorParametersCDFMax$Y <- Training_ImCrit[ , , 2]
  EmulatorParametersCDFMax$w <- function(X){
    return(0.001*diag(dim(as.matrix(X))[1]))
  }
  EmulatorParametersCDFMax$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X ) 
    return(H)}
  EmulatorParametersCDFMax$CorrelationLength <- function(X , n){
    return(0.25*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  {
    Xstar <- Validation_X[ , -3]
    EmulatorOutput <- BE_BayesLinearEmulatorLSEstimatesMO(xstar = Xstar ,EmulatorSettings = EmulatorParametersCDFMax  )
    
    x11()
    plot(xx,abs(EmulatorOutput$E_D_fX[1,] -  Validation_ImCrit[ 1 ,, 2])/diag(sqrt(EmulatorOutput$C_D_fX[1,1]*EmulatorOutput$SigmaHat) ) , type ='l'  , col = rgb(0,0,1,alpha =0.1) , ylim = c(0,5))
    for(i in 2:dim(Validation_ImCrit)[1]){
      lines(xx , abs(EmulatorOutput$E_D_fX[i,] -  Validation_ImCrit[ i ,, 2])/diag(sqrt(EmulatorOutput$C_D_fX[i,i]*EmulatorOutput$SigmaHat) ), col = rgb(0,0,1,alpha =0.1))
    }
    
    Rsq <- 1 - var( (matrix(EmulatorOutput$E_D_fX , length(Validation_ImCrit[ ,  , 2]) , 1 ) -  matrix(Validation_ImCrit[ ,  , 2] , length(Validation_ImCrit[ ,  , 2]) , 1) ) )/var(matrix(Validation_ImCrit[ ,  , 2] , length(Validation_ImCrit[ ,  , 1]) , 1))
    
  }
}

load("C:\\Users\\Ben\\Documents\\Output Images\\AllPatientAnnotations\\HMOutputz1007.RData")

Xstar <- outputstruct[[51]][[3]]$NonImplausibleSets[, -3]
EmulatorOutput <- BE_BayesLinearEmulatorLSEstimatesMO(xstar = Xstar ,EmulatorSettings = EmulatorParametersCDFMax  )

EmulatorOutput$E_D_fX[EmulatorOutput$E_D_fX < 0] <- 0
EmulatorOutput$E_D_fX[EmulatorOutput$E_D_fX > 1] <- 1

tmpmatrix <- EmulatorOutput$E_D_fX
for(i in 1:dim(tmpmatrix)){
  tmpmatrix[i,] <-  tmpmatrix[i,] - (3*sqrt(EmulatorOutput$C_D_fX[i,i]) * diag(EmulatorOutput$SigmaHat))
}

index <- 3
plot(xx , 1-tmpmatrix[index,]  , type ='l' , xlab = 'Mean Implausblity' , ylab = '1 - CDF')
abline( v = outputstruct[[51]][[3]]$Implausability[index,1] , col = 'red')
abline( h = abs(1-tmpmatrix[index,which.min( abs(outputstruct[[51]][[3]]$Implausability[index,1] - xx) )]) , col = 'red')
title('Emulating Observed Implausibility')


dd <- 1 - diag(tmpmatrix[1:dim(Xstar)[1] , apply(as.matrix(outputstruct[[51]][[3]]$Implausability[,1])  , 1 , function(X){which.min(abs(X - xx))}  )]) 
histogram(1-diag(tmpmatrix[1:dim(Xstar)[1] , apply(as.matrix(outputstruct[[51]][[3]]$Implausability[,1])  , 1 , function(X){which.min(abs(X - xx))}  )]) , xlab = 'P(X>Im)' ,main= 'Histogram of P(X>Im)')

BC_PlotPairs(Xstar , alpha = 0.1)
BC_PlotPairs(Xstar , alpha = dd/max(dd))
