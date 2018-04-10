BayesLinearDynamicUpdateAutomatedPrior <- function(RWaveExtractedData , PriorTimePeriod , n)
{
  # Function to produce an automated prior specification from R-R wave time and R wave amplitude data.
  # See tipping points chapter
  # Inputs: RWaveExtractedData: dataframe containing time , R-R wave time and R-wve amplitude.
  #         PriorTimePeriod: time from start used to produce prior specification. 
  #         n: number of points used to calulate adjusted expectation for next time period.
  # Ouputs: List containing second order specification
  Datatocreateprior <-  RWaveExtractedData[RWaveExtractedData$t < (RWaveExtractedData$t[1] + 60*PriorTimePeriod ) ,]
  
  N <- length(Datatocreateprior$t)
  # producing the priors
  
  # pre-allocating memory for structures
  twonmatrix <- matrix(0 , N[1] - (n) ,  2*n)
  
  # Extract data from from automated prior training set
  for(i in 1:(N[1] - n))
  {
    
    twonmatrix[i, 1:n ] <- Datatocreateprior$RA[(1 + (i-1)) :(n + (i-1))]
    twonmatrix[i, (n+1):(2*n) ] <- Datatocreateprior$RR[(1 + (i-1)) :(n + (i-1))]
    
  }
  
  # Calculate second order specification 
  E_D <- as.matrix(apply(twonmatrix , 2 , mean))
  V_D <- cov(twonmatrix)
  E_z <- cbind(E_D[1:(n-1),] ,E_D[ (n+1) : (2*n-1)] )
  V_Z <- V_D[c(-n , -2*n) , ]
  V_Z <- V_Z[ , c(-n , -2*n) ]
  E_x <- E_D[c(n , 2*n)]
  V_x <- V_D[c(n , 2*n) , c(n , 2*n)]
  C_xz <- V_D[c(1:(n-1), (n+1):(2*n -1)) , c(n , 2*n)]
  
  L = t(chol(V_Z))
  W = t(C_xz)%*%solve(t(L) , solve(L))
  V_D_z <- V_x - W%*%(C_xz)
  inv_V_D_z <- solve(V_D_z  )  
  
  outputnames <- list('E_z' , 'V_Z' , 'E_x' , 'V_x' , 'C_xz' , 'V_D_z' , 'inv_V_D_z' , 'W')
  outputs <- list(E_z , V_Z , E_x , V_x , C_xz , V_D_z , inv_V_D_z , W)
  outputs <- setNames(outputs , outputnames)
  return(outputs)
  
}

BayesLinearDynamicUpdateCalulateAdjustedBeliefs <- function(RWaveExtractedData , SecondOrderSpecifiction , n)
{
# Function to calulated adjusted beliefs and diagnstic measures sequentially on an ECG
# Inputs: RWaveExtractedData: dataframe containing time , R-R wave time and R-wve amplitude.
#         Second order specification: list containing a second order specification (output of BayesLinearDynamicUpdateAutomatedPrior)
#         n: number of points used to calulate adjusted expectation for next time period.
# Ouputs: List containing adjusted expectations, standardised residuals and discrepancies.
  
  E_D_z <- matrix(0 , length(RWaveExtractedData$RA) , 2)
  stdresid <-matrix(0 , length(RWaveExtractedData$RA) , 2)
  discrepancyadjustedversion <-matrix(0 , length(RWaveExtractedData$RA) , 1)
  
  
  for(i in 1:length(RWaveExtractedData$RA))
   {
    if(i < (length(RWaveExtractedData$RA) - (n-1)))
    {  
      
      diff <- as.vector(cbind( RWaveExtractedData$RA[(1 + (i-1)):((n-1) + (i-1))] , RWaveExtractedData$RR[(1 + (i-1)):((n-1) + (i-1))])
                        - SecondOrderSpecifiction[["E_z"]] )
      
      E_D_z[i , 1:2] <- SecondOrderSpecifiction[["E_x"]] +  SecondOrderSpecifiction[["W"]]%*%diff
      
      diff <- (E_D_z[i , 1:2] - cbind(RWaveExtractedData$RA[(n-1) + i] , RWaveExtractedData$RR[(n-1) + i]))
      stdresid[i,1:2] <- diff/(sqrt(cbind(SecondOrderSpecifiction[["V_D_z"]][1 , 1] , SecondOrderSpecifiction[["V_D_z"]][2 , 2])))
      discrepancyadjustedversion[i,1] <-diff%*%SecondOrderSpecifiction[["inv_V_D_z"]]%*%t(diff)
    
    }
  }
  
  outputs <- list(E_D_z , stdresid , discrepancyadjustedversion)
  outputnames <- list('E_D_z' , 'std_D_z' ,'dis_D_z')
  outputs <- setNames(outputs , outputnames)
  return(outputs)
}
