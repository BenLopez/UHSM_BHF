DF_mvnpdf <- function(x , mu , Sigma){
  #Sigma <- Sigma + 0.0000000001*diag(dim(Sigma)[1])
  detSigma = det(Sigma) 
  
  diff <- mahalanobis(x , center = mu , cov = Sigma)
  return( -0.5*( log(detSigma) + diff + size(x)[2]*log(2*pi) ) )
}

mvnpdf <- function(x , mu , Sigma){
  return(DF_mvnpdf(x , mu , Sigma))
}


