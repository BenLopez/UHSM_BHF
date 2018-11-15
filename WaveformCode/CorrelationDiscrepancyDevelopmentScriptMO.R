

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

# Multiple output    
X <- matrix(c(1:500) ,500 , 1)
l <- 100
p <- 1
KXX <- CF_ExponentialFamily(X , X , l , p)

mu1 <-  colMeans( DP_RemoveNaRows(DataBase[[3]]) ) 
mu2 <-  colMeans( DP_RemoveNaRows(DataBase[[3]]) ) 
Sigma1 <-  cov( DP_RemoveNaRows(DataBase[[3]]) ) 
Sigma2 <-  1.5*cov( DP_RemoveNaRows(DataBase[[3]]) ) 

SampleGP <- t(apply(BE_SampleSeparableMVGP(KXX , Sigma1)   , 1 , function(X){X + mu1}))
SampleGP <- DataBaseMaster$AFPatientsDatabase[2 , 10001:10500 , 1:11]

test1 <- CD_CalculateActualUpdatedProbabilitiesMOGP(SampleGP , alpha , KXX , mu1 , mu2 , Sigma1 , Sigma2 )
test2 <- CD_CalculateIndividualDensitiesMOPG( SampleGP , mu1 , mu2 , Sigma1 , Sigma2 )

plot(alpha*test2[,1] / (alpha*test2[,1] + (1-alpha)*test2[,2]) , type = 'l' , ylim = c(0,1))
lines(test1 , type = 'l', col = 'red')
