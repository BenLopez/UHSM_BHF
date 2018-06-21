VO_CalulateDistanceMatrix <- function(A , B){
  
 return( matrix(rep(A , length(B)), ncol = length(B)) - t( matrix(rep(B , length(A)), ncol = length(A))) )
 
}
VO_CalulateMinDistance <- function(A , B){
  
  return(apply( abs(AO_CalulateDistanceMatrix(A , B))  , 1 , min))
  
}
