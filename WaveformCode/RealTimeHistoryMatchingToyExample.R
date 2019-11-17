#
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
#
source('RealTimeHistoryMatchingSourefunctions.R')
f = F = function(x){
  
  return(5*x -0.02*x^3 + 1.5*sin(4*x))
  
}

{
  
  {
    x <- seq(0,10,0.001)
    fx <- F(x)
    
    
    xstar <- x[1500] 
    fstar <- f(xstar)
    V_e <- 0.1 
    z <- fstar + rnorm(1 , 0 , sqrt(V_e))
    
    
    x11(20,14)
    p1 <- ggplot(data.frame(x , fx) ) +
      geom_line( aes(x , fx, color="F(x)") ) +
      xlab(TeX('x')) +
      ylab(TeX('z')) + 
      ggtitle(TeX('Toy Simulator and Observation'))  + geom_line(data = data.frame(x,z = (z+0*fx)) , aes(x,z,color = 'z') ) +
      geom_line(data = data.frame(x, z= z + 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) +
      geom_line(data = data.frame(x, z= z - 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) +
      scale_color_manual(values = c(
        "F(x)" = 'blue',
        'z' = 'black',
        'Credible Interval z' = 'red')) +
      labs(color = 'Legend')
    print(p1)
  }
  
  numberofpoints <- 10
  
  
  X <- seq(min(x)+0.5 , max(x)-0.5 , abs(min(x) - max(x))/(numberofpoints-1) )
  fX <- F(X)
  
  Validation_X <- sample(x , numberofpoints)
  Validation_fX <- F(Validation_X)
  
  { 
    PwaveEmulatorParameters <- BE_CreateDefaultEmulationClass()
    PwaveEmulatorParameters$X <- X
    PwaveEmulatorParameters$Y <- fX
    PwaveEmulatorParameters$w <- function(X){
      return(0.000001*var(PwaveEmulatorParameters$Y)*diag(dim(as.matrix(X))[1]))
    }
    PwaveEmulatorParameters$MeanFunction <- function(X){
      X <- as.matrix(X)
      H = cbind(as.matrix(1 + 0*X[,1]) , X , X^2 ) 
      return(H)}
    PwaveEmulatorParameters$CorrelationLength <- function(X , n){
      return(0.05*(apply(X , 2 , max) - apply(X , 2 , min)) )
    }
    { 
      Xstar <- Validation_X
      EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = PwaveEmulatorParameters  )
      
      V_std <- var((Validation_fX - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) ))
      R_sqM <- (1- var(Validation_fX - EmulatorOutput$E_D_MX)/var(Validation_X))
      R_sq <- (1- var(Validation_fX - EmulatorOutput$E_D_fX)/var(Validation_X))
      
      
    }
  }
  
  EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = x ,EmulatorSettings = PwaveEmulatorParameters  )
  
  x11(20,14)
  p2 <- p1 + geom_line(data = data.frame(x, z = EmulatorOutput$E_D_fX),aes(x,z,color = 'E_D[F(x)]')) +
    geom_line(data = data.frame(x, z = EmulatorOutput$E_D_MX),aes(x,z,color = 'mean function')) +
    geom_line(data = data.frame(x, z = EmulatorOutput$E_D_fX- + 3*sqrt(diag(EmulatorOutput$V_D_fX))),aes(x,z,color = 'Credible Interval E_D[F(x)]')) +
    geom_line(data = data.frame(x, z = EmulatorOutput$E_D_fX + 3*sqrt(diag(EmulatorOutput$V_D_fX))),aes(x,z,color = 'Credible Interval E_D[F(x)]')) +
    geom_point(data = data.frame(x=X, z = fX),aes(x , z , color = 'Training Points'))+
    scale_color_manual(values = c(
      "F(x)" = 'blue',
      'z' = 'black',
      'Credible Interval z' = 'red',
      'E_D[F(x)]' = 'darkblue',
      'Training Points' = 'purple',
      'Credible Interval E_D[F(x)]' = 'green',
      'mean function' = 'yellow')) +
    labs(color = 'Legend')
  print(p2)
  
  #Im <- abs(z - EmulatorOutput$E_D_fX)/sqrt(diag(EmulatorOutput$V_D_fX) + V_e)
  #lines(x , DP_RescaleZeroOneToab(Im , min(fx) , max(fx)) , col = 'purple')
}

R1 <- RTHM_EstimateNumberofPoints(x ,ztest =  fx[1000] , fx ,V_e , PwaveEmulatorParameters , doplots = 1,qq=0.01)


R2 <- RTHM_EstimateNumberofPoints(x ,ztest =  z , fx ,V_e , PwaveEmulatorParameters , doplots = 1,qq=0.01)

{
  nn <- 20
  nstruct <- matrix(0,length(seq(1,length(x),nn) ) , 1)
  zstruct <- matrix(0,length(seq(1,length(x),nn) ) , 1)
  muStruct <- matrix(0,length(seq(1,length(x),nn) ) , 1)
  sigmaStruct <- matrix(0,length(seq(1,length(x),nn) ) , 1)
  xstruct <- matrix(0,length(seq(1,length(x),nn) ) , 1)
  
  counter <- 1
  for( ii in seq( 1 , length(x) , nn ) ){
    
    output <- RTHM_EstimateNumberofPoints(x ,ztest =  fx[ii] , fx ,V_e , PwaveEmulatorParameters , doplots = 0,qq=0.01)
    
    nstruct[counter] <- output$nproposed
    zstruct[counter] <- fx[ii]
    muStruct[counter] <- output$muxstar
    sigmaStruct[counter] <- output$sigmaxstar
    xstruct[counter] <- x[ii]
    
    counter <- counter +1
    if(mod(ii , 11) ==0){
      DP_WaitBar(ii/length(x))
    }
  }
}

x11(20,14)
p3 <- ggplot(data.frame(x = xstruct , Nx = nstruct), aes(x , Nx))+
  geom_line(col = 'blue') +
  xlab('x') +
  ylab('N(x)')+
  ggtitle('Calculated N(x)')
p4 <- ggplot(data.frame(x = xstruct , Nx = muStruct), aes(x , Nx))+
  geom_line(col = 'blue') +
  xlab('x') +
  ylab('E[C[x]]')+
  ggtitle('Center Penultimate Non-Implausible Set') 
p5 <- ggplot(data.frame(x = xstruct , Nx = sigmaStruct), aes(x , Nx))+
  geom_line(col = 'blue') +
  xlab('x') +
  ylab('Var[C[x]]') +
  ggtitle('Variance Penultimate Non-Implausible Set') 
grid.arrange(p3,p4,p5,ncol = 1,nrow =3)

# DataBase design

{DataBaseStruct <- matrix(0,0,1)
  
  for(i in 1:length(nstruct)){
    {
      numberofpointsneeded <- nstruct[length(nstruct) +1 - i]  
      ErrorBar <- 3*sqrt(sigmaStruct[length(nstruct) +1 - i])
      LowerBound <- max(0,muStruct[length(nstruct) +1 - i] - ErrorBar)
      UpperBound <- min(10,muStruct[length(nstruct) +1 - i] + ErrorBar)
      
      if(nrow(DataBaseStruct) >0 ){
        numberindatabase <- sum( ( (DataBaseStruct>=LowerBound) & (DataBaseStruct<=UpperBound) ) )
        numberofpointsneeded <- max(0, numberofpointsneeded - numberindatabase)
      }
      
      if(numberofpointsneeded > 0 & (sigmaStruct[length(nstruct) +1 - i]>0) ){
        tmp <- seq(LowerBound,UpperBound , abs(UpperBound - LowerBound)/numberofpointsneeded)
        DataBaseStruct <- rbind(DataBaseStruct,as.matrix(tmp[1:numberofpointsneeded]) )
      }
    }
  }
  
  
  TestStruct <- matrix(0 , length(nstruct),1)
  for(i in 1:length(nstruct)){
    ErrorBar <- 3*sqrt(sigmaStruct[length(nstruct) +1 - i])
    LowerBound <- max(0,muStruct[length(nstruct) +1 - i] - ErrorBar)
    UpperBound <- min(10,muStruct[length(nstruct) +1 - i] + ErrorBar)
    
    TestStruct[length(nstruct) +1 - i] <- sum( ( (DataBaseStruct>=LowerBound) & (DataBaseStruct<=UpperBound) ) )
    
  }
  
  x11(20,14)
  p6 <- ggplot(data.frame(DataBaseStruct) , aes(DataBaseStruct)) + geom_histogram()+
    xlab('x')+
    ylab('Frequency')+
    ggtitle('Histogram of Database')
  
p7 <- ggplot() +
  geom_line(data=data.frame(x=xstruct,y=TestStruct),aes( x , y,color = 'Number of points in C(x) in database'))+
  geom_line(data=data.frame(x=xstruct,y=nstruct),aes( x , y,color = 'N(x)')) +
  xlab('x')+
  ylab('N(x)') +
  ggtitle('Database Validation')
  
grid.arrange(p6,p7,nrow=2)

  }



fxSample <- EmulatorOutput$E_D_fX + BE_SampleGP(DP_AddNugget(EmulatorOutput$V_D_fX , 0.005*diag(diag(EmulatorOutput$V_D_fX)) ))


output <- RTHM_EstimateNumberofPoints(x ,ztest =  z ,fx =  fxSample ,V_e , PwaveEmulatorParameters , doplots = 1,qq=0.01)

x11(20,14)

p8 <- ggplot(data.frame(x = zstruct , Nx = nstruct), aes(x , Nx))+
  geom_line(col = 'blue') +
  xlab('z') +
  ylab('N(x)')+
  ggtitle('Calculated N(x)')

p9 <- ggplot(data.frame(x = zstruct , Nx = muStruct), aes(x , Nx))+
  geom_line(col = 'blue') +
  xlab('z') +
  ylab('E[C[x]]')+
  ggtitle('Center Penultimate Non-Implausible Set') 


p10 <- ggplot(data.frame(x = zstruct , Nx = sigmaStruct), aes(x , Nx))+
  geom_line(col = 'blue') +
  xlab('z') +
  ylab('Var[C[x]]') +
  ggtitle('Variance Penultimate Non-Implausible Set') 

grid.arrange(p8,p9,p10,ncol = 1,nrow =3)


{zz <- z
  
  index <- which.min(abs(zz-zstruct))
  
  ErrorBar <- 6*sqrt( sigmaStruct[index])
  LowerBound <- max( 0 , muStruct[index]  - ErrorBar)
  UpperBound <- min( zz , muStruct[index] + ErrorBar)
  
  
  xx <- sort(DataBaseStruct[(DataBaseStruct>=LowerBound) & (DataBaseStruct<=UpperBound) ])
  fxx <- F(xx)
  
  PwaveEmulatorParameters2 <- BE_CreateDefaultEmulationClass()
  PwaveEmulatorParameters2$X <- xx
  PwaveEmulatorParameters2$Y <- fxx
  PwaveEmulatorParameters2$w <- function(X){
    return(0.000001*var(PwaveEmulatorParameters2$Y)*diag(dim(as.matrix(X))[1]))
  }
  PwaveEmulatorParameters2$MeanFunction <- function(X){
    X <- as.matrix(X)
    H = cbind(as.matrix(1 + 0*X[,1]) , X , X^2 ) 
    return(H)}
  PwaveEmulatorParameters2$CorrelationLength <- function(X , n){
    return(0.1*(apply(X , 2 , max) - apply(X , 2 , min)) )
  }
  Xstar <- seq(min(xx) , max(xx) , (max(xx)-min(xx))/1000)
  EmulatorOutput <- BE_BayesLinearEmulatorLSEstimates(xstar = Xstar ,EmulatorSettings = PwaveEmulatorParameters2  )
  
  
  p11 <- ggplot(data.frame(x =Xstar , fx=F(Xstar)) ) +
    geom_line( aes(x , fx, color="F(x)") ) +
    xlab(TeX('x')) +
    ylab(TeX('z')) + 
    ggtitle(TeX('J(z) Real Time History Match'))  + 
    geom_line(data = data.frame(x=Xstar,z = (z+0*Xstar)) , aes(x,z,color = 'z') ) +
    geom_line(data = data.frame(x=Xstar, z= z + 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) +
    geom_line(data = data.frame(x=Xstar, z= z - 3*sqrt(V_e)),aes(x,z,color = 'Credible Interval z')) + 
    geom_line(data = data.frame(x=Xstar, z = EmulatorOutput$E_D_fX),aes(x,z,color = 'E_D[F(x)]')) +
    geom_line(data = data.frame(x =Xstar , z = EmulatorOutput$E_D_fX- + 3*sqrt(diag(EmulatorOutput$V_D_fX))),aes(x,z,color = 'Credible Interval E_D[F(x)]')) +
    geom_line(data = data.frame(x =Xstar , z = EmulatorOutput$E_D_fX + 3*sqrt(diag(EmulatorOutput$V_D_fX))),aes(x,z,color = 'Credible Interval E_D[F(x)]')) +
    geom_point(data = data.frame(x = xx , z = F(xx)),aes(x , z , color = 'Training Points'))+
    scale_color_manual(values = c(
      "F(x)" = 'blue',
      'z' = 'black',
      'Credible Interval z' = 'red',
      'E_D[F(x)]' = 'darkblue',
      'Training Points' = 'purple',
      'Credible Interval E_D[F(x)]' = 'green')) +
    labs(color = 'Legend') + 
    ylim(c(z - 3*sqrt(V_e) , z + 3*sqrt(V_e)))
x11(20,14)
    print(p11)

}
