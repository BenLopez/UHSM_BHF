pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))

source("LibrariesAndSettings.R" , print.eval  = TRUE )

t <- seq(-pi,pi,0.001)

Rcen <- 0
Rwidth <- 0.05
RA <- 80

Qcen = -0.25
Qwidth <- 0.05
QA <- - 10

Scen =   0.25
Swidth <- 0.05
SA <- -20

Pcen <- -2
Pwidth <- 0.1
PA <- 10
Pslant <- 3
  
Tcen <- 2
Twidth <- 0.2
TA <- 20
Tslant <- -3

Baseline <- 1

ECG <- Baseline + RA*ECGSim_Rpeak( t , Rcen , Rwidth ) +
       QA*ECGSim_Rpeak( t , Qcen , Qwidth  ) +
       SA*ECGSim_Rpeak( t , Scen , Swidth  ) +
       PA*ECGSim_SkewGaussian(t , Pcen , Pwidth , Pslant ) +
       TA*ECGSim_SkewGaussian(t , Tcen , Twidth , Tslant)

plot(t ,  ECG , xlab ='t' , ylab ='Hz')
title('Simulated ECG')


DP_ChooseDataReps()
DP_choosepatient(listAllPatients)
WaveData <- DP_LoadECGReduced(path , subList , numberrep , 1 )

t_observation <- as.numeric((WaveData$Date[100:250] -WaveData$Date[1]))
z <- WaveData$Value[100:250]           
plot(t_observation , z , type = 'l', xlab ='t' , ylab ='Hz')

# Define input names
{
InputNames <- c('Baseline',
                'Rcen',
                'Rwidth',
                'RA',
                'Qcen',
                'Qwidth',
                'QA',
                'Scen',
                'Swidth',
                'SA',
                'Pcen',
                'Pwidth',
                'PA',
                'Pslant',
                'Tcen',
                'Twidth',
                'TA',
                'Tslant')
}

# Manual history match
{ plot(t_observation , z , type = 'l', xlab ='t' , ylab ='Hz')

x = rep(0 , 18)

x[1] <- -8
x[2] <- 0.775000095367432
x[3] <- 0.01
x[4] <- 158
x[5] <- 0.75
x[6] <- 0.01
x[7] <- -30
x[8] <-   0.795
x[9] <- 0.012
x[10] <- -38
x[11] <- 0.67
x[12] <- 0.02
x[13] <- 14
x[14] <- 1
x[15] <- 1.005
x[16] <- 0.05
x[17] <- 20
x[18] <- -2.5

x <- setNames(x ,  InputNames)

ECG <- ECGSim_WrapperSingleBeat(x , t_observation )

lines(t_observation , ECG , col ='red')
title('Simulated ECG')

V_md <- 1*abs( c(0, diff(z) ) )
V_me <- 2

lines(t_observation , z + 3*sqrt(V_md + V_me) , col ='blue')
lines(t_observation , z - 3*sqrt(V_md + V_me) , col ='blue')

Im <- HM_MeanStdError(z , ECG , V_me , V_md)
}

# define bounds for prior nonimplausible set.
{
min_x <- c( -16      ,
            0       ,
            0.0001  ,
            100     ,
            0       ,
            0.0001  ,     
            -60     , 
            0       ,
            0.0001  ,
            -60     ,
            0       ,
            0.00001 ,
            0       , 
            -3      ,
             0      ,
             0.0001 ,
             0      , 
             -3   )

max_x <- c( 16    ,
           1.2    ,
           1      ,
           200    ,
           1.2    ,
           1.5    ,     
           0      , 
           1.2    ,
           1      ,
           0      ,
           1.2    ,
           0.1    ,
           30     , 
           3      ,
           1.2    ,
           0.2    ,
           40     , 
           3   )            

}

{ x = rep(0 , 18)
  
  x[1] <- -8
  x[2] <- 0.775000095367432
  x[3] <- 0.01
  x[4] <- 158
  x[5] <- 0.75
  x[6] <- 0.01
  x[7] <- -30
  x[8] <-   0.795
  x[9] <- 0.012
  x[10] <- -38
  x[11] <- 0.67
  x[12] <- 0.02
  x[13] <- 14
  x[14] <- 1
  x[15] <- 1.005
  x[16] <- 0.05
  x[17] <- 20
  x[18] <- -2.5
  
  min_x <- (x - 0.3*abs(x)) 
  max_x <- (x + 0.3*abs(x)) }

reductionfunction <- function(X ,  InputNames)
{
  X[(   ( X[ , which(InputNames == 'Pcen')] < X[ , which(InputNames == 'Qcen')] )*
        ( X[ , which(InputNames == 'Qcen')] < X[ , which(InputNames == 'Rcen')] )*
        ( X[ , which(InputNames == 'Rcen')] < X[ , which(InputNames == 'Scen')] )*
        ( X[ , which(InputNames == 'Scen')] < X[ , which(InputNames == 'Tcen')] )) ==1 , ]
}


Non_implausible <- t(as.matrix(c(x , Im)))
counter <- 1
Im_thresh <- 2
while( dim(Non_implausible)[1] < 300 )
{
  
numberofsamples <- 1000000
LHSample <- HM_LHD_Reduced(numberofsamples , min_x , max_x , reductionfunction = function(X){ reductionfunction(X , InputNames) } )
colnames(LHSample) <- InputNames

F_matrix <- t(apply(LHSample , 1 , function(X){ECGSim_WrapperSingleBeat( X , t_observation )} ))

Implausability <- apply( F_matrix , 1 , function(X){(HM_MeanStdError(z , X , V_me , V_md)) } )

if(sum(Implausability < Im_thresh) > 1){
  print('Non-implausible point found.')
  Non_implausible <- rbind(Non_implausible , cbind(as.matrix(LHSample[Implausability < Im_thresh,]) , as.matrix(Implausability[Implausability < Im_thresh]))  )}
counter <- counter + 1
print( counter )
}


plot(  t_observation , F_matrix[which.min(apply( F_matrix , 1 , function(X){(HM_MeanStdError(z , X , V_me , V_md))} )),] , type = 'l' , col = 'red')
lines( t_observation , z )
lines( t_observation , z + 3*sqrt( V_md + V_me ) , col ='blue')
lines( t_observation , z - 3*sqrt( V_md + V_me ) , col ='blue')
title( paste0('Minimum Implausability = ' , min(apply( F_matrix , 1 , function(X){HM_MeanStdError(z , X , V_me , V_md)} ))) )
