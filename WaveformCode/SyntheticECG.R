# Script to build a simultor of an ECG signal

ECGSimulator <- function(t, state, parameters)
{
  with(
    as.list(c(state, parameters)), {
    alpha <- 1 - sqrt(X^2 +Y^2) 
    dX <- alpha*X - omega*Y
    dY <- alpha*Y + omega*X
    theta <- atan2(Y , X)
    dZ <- -( ai[1]*( (theta - thetai[1]) %% (2*pi)  )*exp(-(( (theta - thetai[1]) %% (2*pi)  )^2)/(2*bi[1]^2)) )
    -( ai[2]*( (theta - thetai[2]) %% (2*pi)  )*exp(-(( (theta - thetai[2]) %% (2*pi)  )^2)/(2*bi[2]^2) ) )
    -( ai[3]*( (theta - thetai[3]) %% (2*pi)  )*exp(-(( (theta - thetai[3]) %% (2*pi)  )^2)/(2*bi[3]^2) ) )
    -( ai[4]*( (theta - thetai[4]) %% (2*pi)  )*exp(-(( (theta - thetai[4]) %% (2*pi)  )^2)/(2*bi[4]^2) ) )
    -( ai[5]*( (theta - thetai[5]) %% (2*pi)  )*exp(-(( (theta - thetai[5]) %% (2*pi)  )^2)/(2*bi[5]^2) ) ) 
    - Z  
    list(c(dX, dY, dZ))
  }
  )
}



state <-  c(X = 0 , Y=-1 , Z=0)
times  <-  seq(-1 , 1 , 0.005)
parameters <- c(thetai <- c(-(1/3)*pi, -(1/12)*pi , 0 , (1/12)*pi , 0.5*pi ) 
               , ai <- c(1.2 , -5 , 30 ,-7.5 , 0.75) 
               , bi <- c(0.25 , 0.1 , 0.1 , 0.1 , 0.4 )
               , omega <- 2*pi )

out <- ode(y = state, times = times, func = ECGSimulator, parms = parameters) 

par(mfrow = c( 1 , 1))
scatterplot3d(x =out[,2] , y=out[,3] , z = out[,4] )
plot(out[,1] , out[,4])
