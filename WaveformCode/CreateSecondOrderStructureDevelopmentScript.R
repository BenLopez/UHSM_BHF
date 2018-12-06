# Script to assess correlation stcruture.


ProportioninHighRisk2  <-  matrix(0,11,1)
ProportioninLowRisk2   <-  matrix(0,11,1)
SensitivityHighRisk2   <-  matrix(0,11,1)
SpecificityHighRisk2   <-  matrix(0,11,1)
SensitivityLowRisk2    <-  matrix(0,11,1)
SpecificityLowRisk2    <-  matrix(0,11,1)
PPVHighRisk2 <-  matrix(0,11,1)
PPVLowRisk2  <-  matrix(0,11,1)

for(kkk in 1:11){
variableindex <- kkk
rangeoftimes <- 10000:10500
nugget = 0.001 

Data1 <-  BLU_VectoriseData( DataBaseMaster$AFPatientsDatabase[ ,rangeoftimes , variableindex] )
AFSOS <-  BLU_CreateSOS( as.matrix(apply(Data1 , 2 , mean))  ,  cov(Data1) )
AFSOS$WeightMatrix <- BLU_CreateWeightMatrix(AFSOS , nugget*diag(diag(AFSOS$Covariance )))

Data2 <- BLU_VectoriseData( DataBaseMaster$NAFPatientsDatabase[ ,rangeoftimes , variableindex ] )
NAFSOS <-  BLU_CreateSOS( as.matrix(apply(Data2 , 2 , mean)) , cov( Data2) )
NAFSOS$WeightMatrix <- BLU_CreateWeightMatrix(NAFSOS , nugget*diag(diag(NAFSOS$Covariance )))

AdjustedExpectations <-  array(0 , c(dim(Data1)[2] , dim(Data1)[1] , 2 ))
AdjustedVersions <-  array(0 , c(dim(Data1)[2] , dim(Data1)[1] , 2 ))
for(jj in 1:dim(Data1)[1]){
  D <- as.matrix( Data1[jj ,  ] )
  AFSOS <-  BLU_CreateSOS( as.matrix(apply(Data1[-jj,] , 2 , mean))  ,  cov(Data1[-jj,])  )
  print( 'Creating weight matrix.' )
  AFSOS$WeightMatrix <- BLU_CreateWeightMatrix(AFSOS,  nugget = nugget*diag(diag(AFSOS$Covariance )))
  print( 'Weight matrix created.' )
  AdjustedExpectations[ ,jj,1] = BLU_CalculateAdjustedExpectation( AFSOS , D )
  AdjustedExpectations[ ,jj,2] = BLU_CalculateAdjustedExpectation( NAFSOS , D )
  AdjustedVersions[,jj,1] = D - AdjustedExpectations[, jj , 1]
  AdjustedVersions[,jj,2] = D - AdjustedExpectations[, jj , 2] 
  DP_WaitBar(jj/dim(Data1)[1])
}

AFSOS <-  BLU_CreateSOS( as.matrix(apply(Data1 , 2 , mean))  ,  cov(Data1)  )
AFSOS$WeightMatrix <- BLU_CreateWeightMatrix(AFSOS, nugget = nugget*diag(diag(AFSOS$Covariance )))

AdjustedExpectations2 <-  array(0 , c(dim(Data2)[2] , dim(Data2)[1] , 2 ))
AdjustedVersions2 <-  array(0 , c(dim(Data2)[2] , dim(Data2)[1] , 2 ))
jj <- 1

for(jj in 1:dim(Data2)[1]){
  D <- as.matrix( Data2[jj ,  ] )
  NAFSOS <-  BLU_CreateSOS( as.matrix(apply(Data2[-jj,] , 2 , mean))  ,  cov(Data2[-jj,])  )
  print( 'Creating weight matrix.' )
  NAFSOS$WeightMatrix <- BLU_CreateWeightMatrix(NAFSOS,  nugget = nugget*diag(diag(NAFSOS$Covariance )))
  print( 'Weight matrix created.' )
  AdjustedExpectations2[ ,jj,1] = BLU_CalculateAdjustedExpectation( AFSOS , D )
  AdjustedExpectations2[ ,jj,2] = BLU_CalculateAdjustedExpectation( NAFSOS , D )
  AdjustedVersions2[,jj,1] = D - AdjustedExpectations2[,jj,1]
  AdjustedVersions2[,jj,2] = D - AdjustedExpectations2[,jj,2] 
  DP_WaitBar(jj/dim(Data2)[1])
}


AFAdjustedCovariance <- BLU_CalulateAdjustedCovariance(AdjustedVersions[ , , 1] , variableindex , rangeoftimes)
NAFAdjustedCovariance <- BLU_CalulateAdjustedCovariance(AdjustedVersions2[ , , 2] , variableindex , rangeoftimes)

DiscrepanciesAF_AF <- t(matrix(0 , dim(Data1)[1] ,  length(rangeoftimes) ))
DiscrepanciesAF_NAF <- t(matrix(0 , dim(Data1)[1] ,  length(rangeoftimes) ))
DiscrepanciesNAF_NAF <- t(matrix(0 , dim(Data2)[1] ,  length(rangeoftimes) ))
DiscrepanciesNAF_AF <- t(matrix(0 , dim(Data2)[1] ,  length(rangeoftimes) ))

for(jj in 1:max(dim(Data1)[1] , dim(Data2)[1])){
if(jj <= dim(Data1)[1]){
    DiscrepanciesAF_AF[ , jj] <-  BLU_CalulateMultivariateDiscrepancyTimeSeries(as.matrix( AdjustedVersions[  , jj , 1 ] ) ,rangeoftimes ,variableindex ,  AFAdjustedCovariance )
    DiscrepanciesAF_NAF[ , jj] <- BLU_CalulateMultivariateDiscrepancyTimeSeries(as.matrix( AdjustedVersions[  , jj , 2 ] ) ,rangeoftimes ,variableindex , NAFAdjustedCovariance )
}
if(jj <= dim(Data2)[1]){
  DiscrepanciesNAF_AF[ , jj] <-  BLU_CalulateMultivariateDiscrepancyTimeSeries(as.matrix( AdjustedVersions2[  , jj , 1 ] ) ,rangeoftimes ,variableindex ,  AFAdjustedCovariance )
  DiscrepanciesNAF_NAF[ , jj] <- BLU_CalulateMultivariateDiscrepancyTimeSeries(as.matrix( AdjustedVersions2[  , jj , 2 ] ) ,rangeoftimes ,variableindex , NAFAdjustedCovariance )
}
}


plot(apply(DiscrepanciesNAF_NAF ,1 , mean))
points(apply(DiscrepanciesAF_AF ,1 , mean) , col = 'red')

# Calculate culmulative means
CulmulativeStructAF <- apply(abs(apply(DiscrepanciesAF_AF[1:501 ,  ] , 2 , DP_cummean ) - length(variableindex) )  , 1 , mean)
CulmulativeStructNAF <- apply(abs(apply(DiscrepanciesNAF_NAF[1:501 ,  ] , 2 , DP_cummean) - length(variableindex) ) , 1 , mean)
CulmulativeStructAF_NAF <- apply(abs(apply(DiscrepanciesAF_NAF[1:501 ,  ] , 2 , DP_cummean ) - length(variableindex) )  , 1 , mean)
CulmulativeStructNAF_AF <- apply(abs(apply(DiscrepanciesNAF_AF[1:501 ,  ] , 2 , DP_cummean) - length(variableindex) ) , 1 , mean)


plot(CulmulativeStructAF , ylim = c(0,10))
points(CulmulativeStructNAF , col = 'red')
abline(0.175 , 0)

{
plot(c(1:501) ,0.175+0*c(1:501)  , type ='l' , col = rgb(0,0,0, alpha = 1), ylim = c(0,1) , ylab=c(TeX('Probabilities')) , xlab=c(TeX('Beat Number')) )
  
Probabilties <- matrix(0 , 501 , dim(Data2)[1] + dim(Data1)[1] )
Variances <- matrix(0 , 500 , dim(Data2)[1] + dim(Data1)[1] )
counter <- 1
for(j in 1:dim(Data1)[1]){

CulmulativeMean1 = DP_cummean( DiscrepanciesAF_AF[2:501,j] )
CulmulativeMean2 = DP_cummean( DiscrepanciesAF_NAF[2:501,j] )

f_i <- matrix(0 , dim(DiscrepanciesAF_AF)[1] , 2)

for(ii in 2:500){
  f_i[ii,1] <- CD_Tdensity( CulmulativeMean1[ii] , mean =  length(variableindex) , sd = sqrt( CulmulativeStructAF[ii-1] ) , df = ii )
  f_i[ii,2] <- CD_Tdensity( CulmulativeMean2[ii] , mean =  length(variableindex) , sd = sqrt( CulmulativeStructNAF[ii-1] ) , df = ii)
}  

Probabilties[ , counter] = 0.175*f_i[,1]/(0.175*f_i[,1] + (1-0.175)*f_i[,2])
lines(Probabilties[ , counter]  , type ='l' , col = rgb(0,0,1, alpha = 0.1))
counter = counter +1
}

j <- 1
for(j in 1:dim(Data2)[1]){
  #j = sample(1:dim(Data2)[1] , 1)
  CulmulativeMean1 = DP_cummean(DiscrepanciesNAF_AF[2:501,j] )
  CulmulativeMean2 = DP_cummean( DiscrepanciesNAF_NAF[2:501,j] )
  f_i <- matrix(0 , dim(DiscrepanciesNAF_AF)[1] , 2)
  
  for(ii in 2:501){
    f_i[ii,1] <- CD_Tdensity( CulmulativeMean1[ii] , mean =  length(variableindex) , sd = sqrt( CulmulativeStructAF[ii-1] ) , df = ii )
    f_i[ii,2] <- CD_Tdensity( CulmulativeMean2[ii] , mean =  length(variableindex) , sd = sqrt( CulmulativeStructNAF[ii-1] ) , df = ii)
  }  
  Probabilties[ , counter] = 0.175*f_i[,1]/(0.175*f_i[,1] + (1-0.175)*f_i[,2])
  lines(Probabilties[ , counter]  , type ='l' , col = rgb(1,0,0, alpha = 0.1))
  abline(0.175,0)
  counter = counter +1
}
}

plot(apply(Probabilties[,1:dim(Data1)[1]] , 1 , mean), ylim = c(0 , 1) , ylab=c(TeX('Mean Probabilities')), xlab=c(TeX('Beat Number')))
points(apply(Probabilties[,(dim(Data1)[1] + 1):(dim(Data1)[1]+dim(Data2)[1])] , 1 , mean) , col = 'red')
title(TeX('Analysis of Adjusted Probabilities'))
abline( 0.175 , 0)


ProportioninHighRisk2[kkk] <-   sum( apply( Probabilties[100:(length(rangeoftimes) -1),] > 0.95 , 2 , sum ) > 50 )/(dim(Data1)[1] + dim(Data2)[1])
ProportioninLowRisk2[kkk]  <-   sum( apply( Probabilties[100:(length(rangeoftimes) -1),] < 0.05 , 2 , sum ) > 50 )/(dim(Data1)[1] + dim(Data2)[1])
SensitivityHighRisk2[kkk]  <-   sum( apply( Probabilties[100:(length(rangeoftimes) -1),1:dim(Data1)[1]] > 0.95 , 2 , sum ) > 50 )/dim(Data1)[1]
SpecificityHighRisk2[kkk]  <-   (dim(Data2)[1] -sum( apply( Probabilties[100:(length(rangeoftimes) -1),(dim(Data1)[1] + 1):(dim(Data1)[1]+dim(Data2)[1])] > 0.95 , 2 , sum ) > 50  ))/dim(Data2)[1]
SensitivityLowRisk2[kkk]   <-   sum( apply( Probabilties[100:(length(rangeoftimes) -1),(dim(Data1)[1] + 1):(dim(Data1)[1]+dim(Data2)[1])] < 0.05 , 2 , sum ) > 50  ) /dim(Data2)[1] 
SpecificityLowRisk2[kkk]   <-   (dim(Data1)[1] -sum( apply( Probabilties[100:(length(rangeoftimes) -1),1:dim(Data1)[1]] < 0.05 , 2 , sum ) > 50  ))/dim(Data1)[1]
PPVHighRisk2[kkk] <- 0.175*SensitivityHighRisk / (0.175*SensitivityHighRisk + (1-0.175)*(1-SpecificityHighRisk) )
PPVLowRisk2[kkk]  <- (1-0.175)*SensitivityLowRisk / ((1-0.175)*SensitivityLowRisk + (0.175)*(1-SpecificityLowRisk) )
}

CorrectPatientsinHighRisk <- DataBaseMaster$AFPatinetsNames[apply( Probabilties[100:(length(rangeoftimes) -1),1:dim(Data1)[1]] > 0.95 , 2 , sum ) > 50 ]
InCorrectPatientsinHighRisk <- DataBaseMaster$NAFPatinetsNames[apply( Probabilties[100:(length(rangeoftimes) -1),(dim(Data1)[1] + 1):(dim(Data1)[1]+dim(Data2)[1])] > 0.95 , 2 , sum ) > 50 ]
CorrectPatientsinLowRisk <- DataBaseMaster$NAFPatinetsNames[apply( Probabilties[100:(length(rangeoftimes) -1),(dim(Data1)[1] + 1):(dim(Data1)[1]+dim(Data2)[1])] < 0.05 , 2 , sum ) > 50]
InCorrectPatientsinLowRisk <- DataBaseMaster$AFPatinetsNames[ apply( Probabilties[100:(length(rangeoftimes) -1),1:dim(Data1)[1]] < 0.05 , 2 , sum ) > 50 ]


