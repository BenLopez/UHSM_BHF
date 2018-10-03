

x <- seq(-1,1 , 0.01)
KXX <- CF_NeuralNetwork(x = x , xstar = x  , Sigma = diag(c(0.001,100))) 

plot(x , BE_SampleGP(KXX )  , type ='l' , col = rgb(0 , 0 , 1 , alpha = 0.25) , ylim = c(-3,3) , xlab = 'x' , ylab = 'f(x)')
for(i in 1:100){
lines(x , BE_SampleGP(KXX ) , type ='l' , col = rgb(0 , 0 ,1 , alpha = 0.1))
}
