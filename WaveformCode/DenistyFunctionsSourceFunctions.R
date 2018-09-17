mvnpdf <- function(x , mu , Sigma){
  Sigma <- Sigma + 0.0000000001*diag(dim(Sigma)[1])
  A = solve(Sigma )
  detSigma = det(Sigma) 
  
  diff <- apply( x  , 1 , function(X){ ((X - mu)%*%A%*%(X-mu)) } )
  return( -0.5*( log(detSigma) + diff + size(x)[2]*log(2*pi) ) )
}
  


