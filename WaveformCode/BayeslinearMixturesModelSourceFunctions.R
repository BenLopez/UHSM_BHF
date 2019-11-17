{
  BLMM_SampleX <- function(E_pi,E_mu,E_Sigma,V_pi,V_mu,V_Sigma){
    # Sample parameters
    #pi <- BLMM_SampleBeta(E_pi , V_pi)
    #pi <- pi/sum( pi )
    pi <- rdirichlet(1 , E_pi)
    mu <-  rmvnorm(1 ,E_mu ,  V_mu) 
    #sigma <-  abs(rmvnorm(1 ,E_Sigma ,  V_Sigma))
    sigma <- BLMM_SampleGamma(E_Sigma ,  V_Sigma)
    X <- c(pi,mu,sigma,0)  
    return(X)
  }
  BLMM_CaluculateGammaParameters <- function(E_sigma , V_sigma){
    output <- matrix(0,2,length(E_sigma))
    output[1,] <- (E_sigma^2)/diag(V_sigma)
    output[2,] <- diag(V_sigma)/(E_sigma)
    return(output)
  }
  BLMM_SampleGamma <- function(E_Sigma , V_Sigma){
    ncomps <- length(E_Sigma)
    output <- matrix(0 , ncomps,1)
    GammaParameters <- BLMM_CaluculateGammaParameters(E_Sigma , V_Sigma)
    for(i in 1:ncomps){
      output[i] <- rgamma(1,shape =GammaParameters[1,i] , scale = GammaParameters[2,i])
    }
    return(output)
  }
  BLMM_CaluculateBetaParameters <- function(E_pi , V_pi){
    output <- matrix(0,2,length(E_pi))
    output[1,] <- E_pi*(((E_pi*(1-E_pi))/diag(V_pi)) -1)
    output[2,] <- (1 - E_pi)*(((E_pi*(1-E_pi))/diag(V_pi)) -1)
    return(output)
  }
  BLMM_SampleBeta <- function(E_pi , V_pi){
    ncomponents <- length(E_pi)
    BetaParameters <- BLMM_CaluculateBetaParameters(E_pi , V_pi)
    output <- matrix(0 , ncomponents , 1)
    for( i in 1:ncomponents ){
      output[i] <- rbeta(1 , BetaParameters[1 , i] , BetaParameters[2 , i])
    }
    return(output)
  }
  
  BLMM_CalulateSecondorderspecification <- function(E_pi,E_mu,E_Sigma,V_pi,V_mu,V_Sigma , N , n){
    ncomponents <- length(E_pi)
    output <- matrix(0 ,0, 3*ncomponents )
    for(i in 1:n){
      X <- BLMM_SampleX(E_pi , E_mu , E_Sigma , V_pi , V_mu , V_Sigma)
      Z = BLMM_SampleGMM(X , N)
      
      output <- rbind(output , 
                      cbind(Z[,1:ncomponents] ,(matrix(Z$z , N , ncomponents ) - E_mu)/sqrt(E_Sigma + diag(V_mu) ), abs(matrix(Z$z , N , ncomponents ) - E_mu)/sqrt(E_Sigma + diag(V_mu) ))  
                      ) 
    }
    
    output <- list( mu = apply(output , 2 , mean) , CV = cov(output))
    return(  output)
  }
}  
BLMM_SampleGMM <- function( X , N = 250 ){
  
  # sample multinomial distribution for mixtures model
  SampleMultinomial <- t(rmultinom(N  , size = 1 , X[1:3]))
  
  
  output <- matrix(0 , N , 1)
  for(i in 1:dim(SampleMultinomial)[2]){
    if(sum(SampleMultinomial[,i] == 1) ==0){next} # If not sample for componet i: skip
    if( X[10] == 1 ){
      # If GMM 
      output[SampleMultinomial[,i] == 1 , ] <- rnorm(n = sum(SampleMultinomial[,i]) , mean  = X[3+i] , sd = X[6 + i])
    }else{
      # If mixture of GMM and laplace 
      output[SampleMultinomial[,i] == 1 , ] <- FM_rlaplace(n = sum(SampleMultinomial[,i] ) , mu  = X[3+i] , sigma = X[6 + i] , alpha = X[10])
    }
    
  }
  return(data.frame(pi = SampleMultinomial,z = output))
}
