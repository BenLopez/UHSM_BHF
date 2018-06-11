HM_LHD <- function(n , min_x , max_x ){
# Function to create lating hyper cude sample.
# Inputs: n = number of samples min_x:lower bound max_x:upper bound
# output : n xp matrix

  
  if(length(min_x) != length(max_x)){ stop('Dimesnions of bounds do notmatch.') }
  if(sum(min_x > max_x) != 0 ){ stop('Max bound is larger than min bound.') }
  
  return( apply(randomLHS(n , length(min_x) ) , 1 , function(X){ X*(max_x - min_x) + min_x  }) )
  
}

HM_LHD_Reduced<- function( n , min_x , max_x , reductionfunction = function(X){X} ){
 
  if(length(min_x) != length(max_x)){ stop('Dimesnions of bounds do notmatch.') }
  if(sum(min_x > max_x) != 0 ){ stop('Max bound is larger than min bound.') }
  
  Sample <- reductionfunction( t(apply(randomLHS(n , length(min_x) ) , 1 , function(X){ X*(max_x - min_x) + min_x  })) )
  
  while( dim(Sample)[1] < n )
  {
    Sample <- rbind(Sample , reductionfunction( t(apply(randomLHS(n , length(min_x) ) , 1 , function(X){ X*(max_x - min_x) + min_x  })) ))
  } 
  return(Sample)
}

HM_StdError <- function(z , f_x , V_me = 0 , V_md = 0 , V_f = 0){
  if(V_me !=0 || V_md!=0 || V_f !=0 ){ return((z - f_x)/sqrt(V_me + V_md + V_f))}
  if(V_me == 0 & V_md== 0 & V_f == 0 ){ return( (z - f_x) )}
}

HM_MeanStdError <- function(z , f_x , V_me = 0 , V_md = 0 , V_f = 0){
  return(mean(abs(HM_StdError(z , f_x , V_me  , V_md  , V_f )) , na.rm = T))
}

HM_monoIm <- function(z , f_x , V_me = 0 , V_md = 0 , V_f = 0){
  return(abs(HM_StdError(z , f_x , V_me  , V_md  , V_f ) ))
}

HM_maxIm <- function(z , f_x , V_me = 0 , V_md = 0 , V_f = 0){
 return( max(HM_monoIm(z , f_x , V_me , V_md , V_f ) , na.rm = TRUE) )
} 
