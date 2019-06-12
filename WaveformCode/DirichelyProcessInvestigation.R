
{
  if(file.exists('CheckforDefaultsScript.R')){
    source('CheckforDefaultsScript.R')
  }else{
    pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
    source("LibrariesAndSettings.R" , print.eval  = TRUE )
    DP_LoadPatientIndex()
    DP_ChooseDataReps()
    FilestoProcess <- DP_ChooseECGstoProcess() 
    HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
  }
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}
FilestoProcess = DP_ChooseECGstoProcess()

source('FM_CreateRhythumPriors.R')

PNIS <- PriorNonImplausibleSetRegularyIreRegular
F_x <- F_x_ReIre
f_x <- f_x_ReIre
  
{

{
  c <- 1000
  l <- 0.005
  n <- 500 - 19
  G_0 <- function(n) { FM_SampleGMM( X = PNIS[index,] , n) } 
  index <- 1
}
  
 
  {P1 <- ggplot(data.frame(x , y = F_x[index,]) , aes(x , y) ) + geom_line(col = 'blue') +
        geom_line(data = data.frame(x , y = F_x[index,]+ 3*sqrt( (F_x[index,]*(1-F_x[index,]))/(c + 1) + (F_x[index,]*(1-F_x[index,]))/(n) + 0.005^2)) , aes(x , y) , col ='red' ) +
        geom_line(data = data.frame(x , y = F_x[index,]- 3*sqrt( (F_x[index,]*(1-F_x[index,]))/(c + 1) + (F_x[index,]*(1-F_x[index,]))/(n) + 0.005^2)) , aes(x , y) , col ='red' ) +
        ggtitle('Samples of Observations for Heart-rhythm') +
        xlab('RR') +
        ylab('Culmulative Probability')
  for(i in 1:10){
    y <- G_0(n)
    theta <- FM_SampleDP(c , l , n , G_0)
    
    A <- sqrt(var(theta)/(var(theta) + l^2))
    theta <- A*(theta + rnorm(n = n  , mean = 0 , l )) + (1 - A)*mean(theta)
    P1 <- P1 + geom_line(data = data.frame(x , y = FM_CalculateCDFS(theta , xx = x)) , aes(x , y) , col =rgb(0,0,0,alpha = 0.5))
  }
    x11(10,10)
    print(P1)
    } 
  
}

x11()
par(mfrow = c(2 , 2))

{
 
{
  c <- 10000
  l <- 0.005
  G_0 <- function(n) { FM_SampleGMM( X = PNIS[index,] , n)}
  n <- 5000 - 19
  
}
  index <- 1
  plot(x   , f_x[index,] , type ='l' , ylim = c(0,max(f_x[index,])+0.5 ) , xlab = 'RR' , ylab = 'Density', xlim = c(0.25,2))
  title(TeX(paste0('Samples of DP: $\\alpha = ' , c  , '$')))
{
  p2 <- ggplot(data.frame(x = x, y = f_x[index,]) , aes(x , y) ) + geom_line(col = 'black' , size = 1) + xlab('t') + ylab('Density')  
  p3 <- p2
  
  c <- 100
for(i in 1:10){
  y <- G_0(n)
  theta <- FM_SampleDP(c , l , n , G_0)
  
  kde1 <- kde( theta )
  kde2 <- kde( y )
  
  p3 <- p3 + geom_line(data = data.frame(x = x, y =  predict( object = kde1 , x = x)) , aes(x , y) , col =rgb(0,0,1,alpha = 0.5))
  #lines(x - 0.005 , predict( object = kde2 , x = x- 0.005 ) , type = 'l' , col =rgb(1,0,0,alpha = 0.1))
}
  p4 <- p2
  
  c <- 500
  for(i in 1:10){
    y <- G_0(n)
    theta <- FM_SampleDP(c , l , n , G_0)
    
    kde1 <- kde( theta )
    kde2 <- kde( y )
    
    p4 <- p4 + geom_line(data = data.frame(x = x, y =  predict( object = kde1 , x = x )) , aes(x , y) , col =rgb(0,0,1,alpha = 0.5))
    #lines(x - 0.005 , predict( object = kde2 , x = x- 0.005 ) , type = 'l' , col =rgb(1,0,0,alpha = 0.1))
  }
  c <- 1000
  p5 <- p2
  
for(i in 1:10){
    y <- G_0(n)
    theta <- FM_SampleDP(c , l , n , G_0)
    
    kde1 <- kde( theta )
    kde2 <- kde( y )
    
    p5 <- p5 + geom_line(data = data.frame(x = x, y =  predict( object = kde1 , x = x )) , aes(x , y) , col =rgb(0,0,1,alpha = 0.5))
    #lines(x - 0.005 , predict( object = kde2 , x = x- 0.005 ) , type = 'l' , col =rgb(1,0,0,alpha = 0.1))
}
  c <- 10000
  p6 <- p2
  for(i in 1:10){
    y <- G_0(n)
    theta <- FM_SampleDP(c , l , n , G_0)
    
    kde1 <- kde( theta )
    kde2 <- kde( y )
    
    p6 <- p6 + geom_line(data = data.frame(x = x, y =  predict( object = kde1 , x = x- 0.005 )) , aes(x , y) , col =rgb(0,0,1,alpha = 0.5))
    #lines(x - 0.005 , predict( object = kde2 , x = x- 0.005 ) , type = 'l' , col =rgb(1,0,0,alpha = 0.1))
  }
  }
}

x11()
grid.arrange( p3 + ggtitle(TeX('$\\alpha = 100$')),p4+ ggtitle(TeX('$\\alpha = 500$')),p5+ ggtitle(TeX('$\\alpha = 1000$')),p6+ ggtitle(TeX('$\\alpha = 10000$')) , nrow = 2, ncol = 2 ,top = textGrob("Samples from DPs",gp=gpar(fontsize=20,font=3)) )

