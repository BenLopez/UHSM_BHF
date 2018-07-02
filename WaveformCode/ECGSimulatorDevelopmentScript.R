{pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
source("LibrariesAndSettings.R" , print.eval  = TRUE )}

t <- seq(-pi,pi,0.04)

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

ECG <- Baseline + RA*ECGSim_Gaussian( t , Rcen , Rwidth ) +
       QA*ECGSim_Gaussian( t , Qcen , Qwidth  ) +
       SA*ECGSim_Gaussian( t , Scen , Swidth  ) +
       PA*ECGSim_SkewGaussian(t , Pcen , Pwidth , Pslant ) +
       TA*ECGSim_SkewGaussian(t , Tcen , Twidth , Tslant)

plot(t ,  ECG , xlab ='t' , ylab ='Hz' , type ='l' )
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
                'PcenL',
                'PwidthL',
                'PAL',
                'PcenR',
                'PwidthR',
                'PAR',
                'TcenL',
                'TwidthL',
                'TAL',
                'TcenR',
                'TwidthR',
                'TAR')
}

# Manual history match
{ plot(t_observation , z , type = 'l', xlab ='t' , ylab ='Hz')

x = rep(0 , 18)

x[1] <- -8 # Baseline
x[2] <- 0.775000095367432 #Rceb
x[3] <- 0.01 #Rwidth
x[4] <- 164 #RA
x[5] <- 0.755 #Qcen
x[6] <- 0.01 #Qwidth
x[7] <- -23 #QA
x[8] <-   0.79 #Scen
x[9] <- 0.012 #Swidth
x[10] <- -30 #SA
x[11] <- 0.62 # PcenL
x[12] <- 0.01 # PWidthL
x[13] <- 8 # PAL
x[14] <- 0.67 #PcenR
x[15] <- 0.02 #PenR
x[16] <- 12 #PAR
x[17] <- 1.01 #TcenL
x[18] <- 0.03 #TwidthL
x[19] <- 11 #TwidthL
x[20] <- 1.02 #TcenR
x[21] <- 0.03 #TwidthR
x[22] <- 8 #TAR

x <- setNames(x ,  InputNames)

ECG <- ECGSim_WrapperSingleBeat_m2(x , t_observation )

lines(t_observation , ECG , col ='red')
title('Simulated ECG')

V_md <- 1*abs( c(0, diff(z) ) )
V_me <- 2

#lines(t_observation , z + 3*sqrt(V_md + V_me) , col ='blue')
#lines(t_observation , z - 3*sqrt(V_md + V_me) , col ='blue')

legend(1, 95, legend=c("Non-implausible Simulated Signal", "Observed Patient ECG Signal"),
       col=c("red", "black"), lty=c(1 , 1), cex=0.8)

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
  
  x[1] <- -8 # Baseline
  x[2] <- 0.775000095367432 #Rceb
  x[3] <- 0.01 #Rwidth
  x[4] <- 164 #RA
  x[5] <- 0.755 #Qcen
  x[6] <- 0.01 #Qwidth
  x[7] <- -23 #QA
  x[8] <-   0.79 #Scen
  x[9] <- 0.012 #Swidth
  x[10] <- -30 #SA
  x[11] <- 0.62 # PcenL
  x[12] <- 0.01 # PWidthL
  x[13] <- 8 # PAL
  x[14] <- 0.67 #PcenR
  x[15] <- 0.02 #PenR
  x[16] <- 12 #PAR
  x[17] <- 1.01 #TcenL
  x[18] <- 0.03 #TwidthL
  x[19] <- 11 #TwidthL
  x[20] <- 1.02 #TcenR
  x[21] <- 0.03 #TwidthR
  x[22] <- 8 #TAR
  
  x <- setNames(x ,  InputNames)
  
  
  min_x <- (x - 0.3*abs(x)) 
  max_x <- (x + 0.3*abs(x)) }

reductionfunction <- function(X ,  InputNames)
{
  X[(   ( X[ , which(InputNames == 'PcenL')] < X[ , which(InputNames == 'PcenR')] )*   
        ( X[ , which(InputNames == 'PcenR')] < X[ , which(InputNames == 'Qcen')] )*
        ( X[ , which(InputNames == 'Qcen')] < X[ , which(InputNames == 'Rcen')] )*
        ( X[ , which(InputNames == 'Rcen')] < X[ , which(InputNames == 'Scen')] )*
        ( X[ , which(InputNames == 'Scen')] < X[ , which(InputNames == 'TcenL')] )*
        ( X[ , which(InputNames == 'TcenL')] < X[ , which(InputNames == 'TcenR')] )) ==1 , ]
}


Non_implausible <- t(as.matrix(c(x , Im)))
counter <- 1
Im_thresh <- 2
while( dim(Non_implausible)[1] < 500 )
{
  
numberofsamples <- 1000000
LHSample <- HM_LHD_Reduced(numberofsamples , min_x , max_x , reductionfunction = function(X){ reductionfunction(X , InputNames) } )
colnames(LHSample) <- InputNames

F_matrix <- t(apply(LHSample , 1 , function(X){ECGSim_WrapperSingleBeat_m2( X , t_observation )} ))

Implausability <- apply( F_matrix , 1 , function(X){(HM_MeanStdError(z , X , V_me , V_md)) } )

if(sum(Implausability < Im_thresh) > 1){
  print('Non-implausible point found.')
  Non_implausible <- rbind(Non_implausible , cbind(as.matrix(LHSample[Implausability < Im_thresh,]) , as.matrix(Implausability[Implausability < Im_thresh]))  )
  }
counter <- counter + 1
DP_WaitBar( dim(Non_implausible)[1]/500 )
}

MeanNonIm <- apply( Non_implausible[, 1:18] , 2, mean )
NonIm_F <- apply( Non_implausible[, 1:18] , 1 , function(X){ECGSim_WrapperSingleBeat(X , t_observation)} )


x11(25,20)
colorvector = matrix(0,  1128  , 1)
colorvector[1:128] = rgb(1,0,0 , alpha = 0.1)
colorvector[129:length(colorvector)] = rgb(1,0,0 , alpha = 0.0001)

pairs( rbind(Non_implausible[ , 1:5] , LHSample[1:1000 , 1:5]) , pch = 16 ,  col = colorvector )
title('Posterior Non-Implausible Sets')

pairs( rbind(LHSample[1:1000 , 1:5]) , pch = 16 ,  col = rgb(1,0,0, 0.1) )
title('Posterior Non-Implausible Sets')


plot(  t_observation ,  apply(NonIm_F , 1 , mean) , type = 'l' , col = 'red' , ylim = c(-40,170) , xlab = c('t') , ylab = 'Hz')
lines( t_observation , z )
lines( t_observation , z + 3*sqrt( V_md + V_me ) , col ='blue')
lines( t_observation , z - 3*sqrt( V_md + V_me ) , col ='blue')
title( paste0('Mean Non-Implausibe Output')  )

abline( v = MeanNonIm[2] , col ='yellow' )
abline( v =  MeanNonIm[5] , col ='yellow')
abline( v =  MeanNonIm[8] , col ='yellow')
abline( v =  MeanNonIm[11] , col ='yellow')
abline( v =  MeanNonIm[15] , col ='yellow')

abline( h = MeanNonIm[4] , col ='green' )
abline( h =  MeanNonIm[7] , col ='green')
abline( h =  MeanNonIm[8] , col ='green')
abline( h =  MeanNonIm[10] , col ='green')
abline( h =  MeanNonIm[13] , col ='green')
abline( h =  MeanNonIm[17] , col ='green')






