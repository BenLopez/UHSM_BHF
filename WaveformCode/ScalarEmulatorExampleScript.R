source('EmulationSourceCode.R')

x <- as.matrix(c(1:5))/5
y <- as.matrix( 1 + 0.5*x + 0.1*sinpi(2*x) )
xstar <- as.matrix(c(1:1000))/1000

par(mfrow = c(1 , 1))
plot(x , y , xlab = 't' , ylab = 'f(t)')
lines(xstar ,  1 + 0.5*xstar + 0.1*sinpi(2*xstar)  , type = 'l')
h <- function(x){cbind(1 +0*x , x)}
l <- 0.2
p <- 2
w <- 0.00000001

AdjustedBeliefs <- BayesLinearEmulatorGLSEstimates(y , x , xstar , w , l , p , h )

lines( xstar , AdjustedBeliefs[[1]] , col = 'blue')
lines( xstar , AdjustedBeliefs[[1]] + 2*sqrt(diag(AdjustedBeliefs[[2]])) , col = 'red')
lines( xstar , AdjustedBeliefs[[1]] - 2*sqrt(diag(AdjustedBeliefs[[2]])), col = 'red')
title('Adjusted Beliefs for Toy Simulator')

legend(0.2, 1.5, legend=c("Credible Interval", "Adjusted Expectation" , "Known Truth" , "Training Points"),
       col=c("red", "blue" , "black" , 'black'),lwd=1, pch = c(NA , NA , NA , 1), lty=c(1 , 1 , 1 , 0 ), cex=0.8)

