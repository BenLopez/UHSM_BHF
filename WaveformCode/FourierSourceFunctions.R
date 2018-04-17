fftshift <- function(x)
{
  # Function to shift zero frequency component to centre.
  n <- length(x)
  x <- c( x[(n/2+1):n] , x[1:(n/2)] )
  return(x)
}

SetElementsOfListoToZero <- function(A , a)
{
  for (i in 1:length(a))
  {
    A[[a[i]]] <- 0*A[[a[i]]]
  }
  return(A)
}

imfiter1D <- function(f_tt , Filter)
{
 # Function to convolve a signal with a filter.
  output <- abs(convolve(f_tt , Filter , type = "open"))
  output <- output[(length(Filter)/2):(length(output)-(length(Filter)/2))]
  return(output)
} 