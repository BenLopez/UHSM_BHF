




plot(x, t(L) %*% rnorm(dim(x)[1])  , type ='l' , col = rgb(1 , 0 , 0 , alpha = 0.5) , ylim = c(-2,2))
lines(x, t(L) %*% rnorm(dim(x)[1])  , type ='l' , col = rgb(0 , 0 ,1 , alpha = 0.5))

