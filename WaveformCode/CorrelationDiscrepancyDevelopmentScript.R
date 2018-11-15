# Script to analyse the affect of discrepancy on the calulations

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

{mu1 <- 0.2
 v1 <-  0.4

mu2 <- mu1 + 0.8
v2 <-  1.1*v1

X <- matrix(c(1:500) ,500 , 1)
l <- 10
p <- 1.45
alpha = 0.1
df = 5
KXX <- CF_ExponentialFamily(X , X , l , p)

SampleGP1 <- mu1 + BE_SampleGP( v1*KXX )
SampleGP2 <- mu2 + BE_SampleGP( v2*KXX )

plot(SampleGP1 , type='l' , col = rgb(1,0,0,alpha = 0.5))
lines(SampleGP2 , type='l' , col = rgb(0,0,1,alpha = 0.5))

invD = CD_CalculateInverseVarDStack(KXX)

# Individual probabilities

f_i = matrix(0, dim(SampleGP1)[1] , 2)
f_i[ , 1] <- dnorm( SampleGP1 , mean = mu1  , sd = sqrt(v1))
f_i[ , 2] <- dnorm( SampleGP1 , mean = mu2  , sd = sqrt(v2))

IndividualProbabilities = (alpha*f_i[,1]) / ((1 - alpha)*f_i[,2] + alpha*f_i[,1])

# Create specification structures
Specification1 <- CD_CreateDefaultSpecification1D( l , p , mu1 , v1 ) 
Specification2 <- CD_CreateDefaultSpecification1D( l , p , mu2 , v2 ) 
}

# Calulate Actual Probailties
ActualProbabilities <- CD_CalculateActualUpdatedProbabilities( SampleGP1 , alpha , KXX , mu1 , mu2 , v1 , v2 )

p1 <- ggplot( data = data.frame(time = c(1:501) , Probabilities = ActualProbabilities[,1] ) , aes(time , Probabilities) ) +
      geom_line( color = rgb(0,0,1 , alpha = 0.5)) +
      ggtitle('Actual Probabilties')

x11(20 , 14)
print(p1)

numberofsamples <- 10
SetofPs <- seq( 1.25 , 1.65, 0.05 )
Setofls <- seq( 6 , 14, 1 )

Probabilities <- matrix(0,Specification1$n + 1 ,  length(SetofPs)*length(Setofls) )
counter <- 1
for(i in 1:length(SetofPs ) ){
for(j in 1:length( Setofls ) ){
Specification1$p  <-  SetofPs[ i ]
Specification1$l  <-  Setofls[ j ]
Specification2$p <- SetofPs[ i ]
Specification2$l <- SetofPs[ i ]
Probabilities[  , counter] <- CD_CalculateUpdatedProbabilitiesFromSpecifications(SampleGP1 , Specification1 , Specification2 , alpha = alpha)
counter <- counter +1
DP_WaitBar(counter/(length(SetofPs)*length(Setofls)))
}
}

p2 <- ggplot(  ) +
      geom_line(data = data.frame(time = c(1:501) , Probabilities = ActualProbabilities[,1] ) , aes(time , Probabilities) ,  color = rgb(0,0,1 , alpha = 0.5)) +
      ggtitle('Sensitivity of Conditional Probabilties to Discrepancy in Correlation Specification')
for(i in 1:dim(Probabilities)[2]){
p2 <- p2 +
  geom_line(data = data.frame(time = c(1:501) , Probabilities =Probabilities[,i] ) , aes(time , Probabilities) ,  color = rgb(1,0,0 , alpha = 0.1))
}
print(p2)

Specification1 <- CD_CreateDefaultSpecification1D( l , p , mu1 , v1 ) 
Specification2 <- CD_CreateDefaultSpecification1D( l , p , mu2 , v2 ) 
SOS <- CD_CalculateSecondOrderSpecification( Specification1 , Specification2 , alpha = alpha  , numberofsamples = 1000 , numberinupdate = 25)
AdjustedBeliefPrevision <- CD_CalulateAdjustedBeliefsForPrevision1D( SOS , as.matrix(SampleGP1) , Specification1 , Specification2 ,  alpha = alpha )

{x11(20,14)
p3 <- ggplot(  ) +
  geom_line(data = data.frame(time = c(1:501) , Probabilities = ActualProbabilities[,1] ) , aes(time , Probabilities) ,  color = rgb(0,0,1 , alpha = 0.5)) +
  ggtitle('Adjusted Beliefs for Actual Probabilities')
tmp1 <- AdjustedBeliefPrevision[,1] + 3*sqrt(AdjustedBeliefPrevision[,2]) 
tmp2 <- AdjustedBeliefPrevision[,1] - 3*sqrt(AdjustedBeliefPrevision[,2]) 
tmp1[tmp1 > 1] <- 1
tmp2[tmp2 < 0] <- 0
p3 <- p3 +
  geom_line(data = data.frame(time = c(1:501) , Probabilities = AdjustedBeliefPrevision[,1] ) , aes(time , Probabilities) ,  color = rgb(1,0,0 , alpha = 0.5))+
  geom_line(data = data.frame(time = c(1:501) , Probabilities =  tmp1), aes(time , Probabilities) ,  color = rgb(0,0,0 , alpha = 0.5))+
  geom_line(data = data.frame(time = c(1:501) , Probabilities =  tmp2 ) , aes(time , Probabilities) ,  color = rgb(0,0,0 , alpha = 0.5))
print(p3)}

# Fix specification vary update

counter = 1
SetofPs <- seq( 1.25 , 1.65, 0.05 )
Setofls <- seq( 6 , 14, 1 )
ActualProbabilities <- array(0,c(Specification1$n + 1 , numberofsamples,  length(SetofPs)*length(Setofls)) )
PredictedProbabilities <- array(0,c(Specification1$n + 1 , numberofsamples,  length(SetofPs)*length(Setofls)) )
AdjustedBeliefPrevision <- array(0,c(Specification1$n + 1 , numberofsamples,  length(SetofPs)*length(Setofls)) )
Specification1 <- CD_CreateDefaultSpecification1D( l , p , mu1 , v1 ) 
Specification2 <- CD_CreateDefaultSpecification1D( l , p , mu2 , v2 ) 
numberofsamples <- 10
counter <- 1
for(i in 1:length(SetofPs ) ){
  for(j in 1:length( Setofls ) ){
    Specification1$p  <-  SetofPs[ i ]
    Specification1$l  <-  Setofls[ j ]
    Specification2$p  <-  SetofPs[ i ]
    Specification2$l  <-  Setofls[ j ]
    SetofSamples <- cbind(CD_SampleDataFromSpecification1D(Specification = Specification1 , numberofsamples = round(alpha*numberofsamples)) , CD_SampleDataFromSpecification1D(Specification = Specification2 , numberofsamples = round((1-alpha)*numberofsamples)))
    ActualProbabilities[ , , counter] <- CD_CalculateUpdatedProbabilitiesFromSpecifications( as.matrix(SetofSamples) , Specification1 , Specification2 , alpha = alpha )
    Specification1 <- CD_CreateDefaultSpecification1D( l , p , mu1 , v1 ) 
    Specification2 <- CD_CreateDefaultSpecification1D( l , p , mu2 , v2 ) 
    PredictedProbabilities[ , , counter] <-  CD_CalculateUpdatedProbabilitiesFromSpecifications( as.matrix(SetofSamples) , Specification1 , Specification2 , alpha = alpha )
    AdjustedBeliefPrevision[, ,counter] <- CD_CalulateAdjustedBeliefsForPrevision1D( SOS , as.matrix(SetofSamples) , Specification1 , Specification2 ,  alpha = alpha )[,1]
    DP_WaitBar(counter /  (length(SetofPs)*length(Setofls)) )
    counter = counter +1 
}
} 
  
{x11(20,14)
p4 <- ggplot() +
ggtitle(TeX('Difference Between Estimate for Probability and Actual: Discrepancy in Correlation Function'))

for(j in 1:numberofsamples){
for(i in 1:length(SetofPs)*length(Setofls)){
  p4 <- p4 + geom_line(data = data.frame(time = c(1:501) , Probabilities = ActualProbabilities[,j,i] - PredictedProbabilities[,j,i]) ,  aes(time , Probabilities)  , color = rgb(0,0,1 ,alpha = 0.1))
}
}

for(j in 1:numberofsamples){  
for(i in 1:length(SetofPs)*length(Setofls)){
  p4 <- p4 + geom_line(data = data.frame(time = c(1:501) , Probabilities = ActualProbabilities[,j,i] - AdjustedBeliefPrevision[,j,i]) ,  aes(time , Probabilities)  , color = rgb(1,0,0 ,alpha = 0.1))
}
}
print(p4)

muhat1 <- apply(matrix(ActualProbabilities - PredictedProbabilities , dim(ActualProbabilities)[1] ,dim(ActualProbabilities)[2]*dim(ActualProbabilities)[3] ) , 1 , mean)
vhat1 <- apply(matrix(ActualProbabilities - PredictedProbabilities , dim(ActualProbabilities)[1] ,dim(ActualProbabilities)[2]*dim(ActualProbabilities)[3] ) , 1 , var)
muhat2 <- apply(matrix(ActualProbabilities - AdjustedBeliefPrevision , dim(ActualProbabilities)[1] ,dim(ActualProbabilities)[2]*dim(ActualProbabilities)[3] ) , 1 , var)
vhat2 <- apply(matrix(ActualProbabilities - AdjustedBeliefPrevision , dim(ActualProbabilities)[1] ,dim(ActualProbabilities)[2]*dim(ActualProbabilities)[3] ) , 1 , var)

x11(20,14)
p5 <- ggplot() +
  ggtitle(TeX('Sample Mean of Difference Between Estimated and Actual Probabilities')) +
  geom_line(data = data.frame(time = c(1:501) , Probabilities = muhat1) ,  aes(time , Probabilities)  , color = rgb(0,0,1 ,alpha = 0.5)) +
  geom_line(data = data.frame(time = c(1:501) , Probabilities = muhat2) ,  aes(time , Probabilities)  , color = rgb(1,0,0 ,alpha = 0.5)) 
print(p5)
x11(20,14)
p6 <- ggplot() +
  ggtitle(TeX('Sample Variance of Difference Between Estimated and Actual Probabilities')) +
  geom_line(data = data.frame(time = c(1:501) , Probabilities = vhat1) ,  aes(time , Probabilities)  , color = rgb(0,0,1 ,alpha = 0.5)) +
  geom_line(data = data.frame(time = c(1:501) , Probabilities = vhat2) ,  aes(time , Probabilities)  , color = rgb(1,0,0 ,alpha = 0.5)) 
print(p6)

}

# T process. 
numberofsamples <- 100
SetofSamples <- cbind(CD_SampleDataFromSpecificationTP(Specification = Specification1 , numberofsamples = round(alpha*numberofsamples)) , CD_SampleDataFromSpecificationTP(Specification = Specification2 , numberofsamples = round((1-alpha)*numberofsamples)))
ActualProbabilities <- CD_CalculateUpdatedProbabilitiesFromSpecificationsTP(SetofSamples , Specification1 , Specification2 , alpha = alpha , df = df )
EstimatedProbabilities <- CD_CalculateUpdatedProbabilitiesFromSpecifications(SetofSamples , Specification1 , Specification2 , alpha = alpha )
AdjustedProbabilties <- CD_CalulateAdjustedBeliefsforPrevisonFromSpecifications(SOS,SetofSamples , Specification1 , Specification2 , alpha = alpha )

{
x11(20,14)
p7 <- ggplot() +
  ggtitle( TeX('Difference Between Estimate for Probability and Actual: Discrepancy in Likelihood') )

  for(i in 1:numberofsamples){
    p7 <- p7 + geom_line(data = data.frame(time = c(1:501) , Probabilities = ActualProbabilities[,i] - EstimatedProbabilities[,i]) ,  aes(time , Probabilities)  , color = rgb(0,0,1 ,alpha = 0.1))
  }


for(i in 1:numberofsamples){
    p7 <- p7 + geom_line(data = data.frame(time = c(1:501) , Probabilities = ActualProbabilities[,i] - AdjustedProbabilties[,i]) ,  aes(time , Probabilities)  , color = rgb(1,0,0 ,alpha = 0.1))
  }
print(p7)
}






