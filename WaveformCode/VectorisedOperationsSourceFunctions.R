VO_CalulateDistanceMatrix <- function(A , B){
  
 return( matrix(rep(A , length(B)), ncol = length(B)) - t( matrix(rep(B , length(A)), ncol = length(A))) )
 
}
VO_CalulateMinDistance <- function(A , B){
  
  return(apply( abs(AO_CalulateDistanceMatrix(A , B))  , 1 , min))
  
}
VO_rollvar <- function(X , k = 1000){
  output <- rollmean(X^2 , k = k , na.pad = TRUE , align = 'center') - rollmean(X , k = k , na.pad = TRUE , align = 'center')^2
  output[output <0] <- 0
  return(output)
}