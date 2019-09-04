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

source('FM_CreatePWavePriors.R')


{
  sigma_f = 0.5
x <- seq(0,1,0.001)
f_x <- 20*PsimulatorFunction(PriorNonImplausibleSet[2,] , x)
nugget <-  0.00001
l <- 0.01
KXX <- sigma_f*DP_AddNugget(CF_ExponentialFamily(x , x , l , 2), nugget)
y <- f_x + BE_SampleGP(KXX)

p1 <- ggplot(data.frame(x , y) , aes(x , y)) + geom_line(col = rgb(0,0,1,alpha = 0.5)) + ggtitle('l = 0.01') + xlab('t') + ylab(TeX('$\\epsilon_{md}(t)$'))
for(i in 1:9){
  y <- f_x + BE_SampleGP(KXX)
  p1 <- p1 + geom_line(data = data.frame(x , y) , aes(x , y) , col = rgb(0,0,1,alpha = 0.5))  
}


l <- 0.02
KXX <- sigma_f*DP_AddNugget(CF_ExponentialFamily(x , x , l , 2), nugget)
y <- f_x+BE_SampleGP(KXX)

p2 <- ggplot(data.frame(x , y) , aes(x , y)) + geom_line(col = rgb(0,0,1,alpha = 0.5))  + ggtitle('l = 0.02') + xlab('t') + ylab(TeX('$\\epsilon_{md}(t)$'))
for(i in 1:9){
  y <- f_x+BE_SampleGP(KXX)
  p2 <- p2+ geom_line(data = data.frame(x , y) , aes(x , y) , col = rgb(0,0,1,alpha = 0.5))  
}

l <- 0.04
KXX <- sigma_f*DP_AddNugget(CF_ExponentialFamily(x , x , l , 2), nugget)
y <- f_x+BE_SampleGP(KXX )

p3 <- ggplot(data.frame(x , y) , aes(x , y)) + geom_line(col = rgb(0,0,1,alpha = 0.5))  + ggtitle('l = 0.04') + xlab('t') + ylab(TeX('$\\epsilon_{md}(t)$'))
for(i in 1:9){
  y <- f_x+BE_SampleGP(KXX)
  p3 <- p3 + geom_line(data = data.frame(x , y) , aes(x , y) , col = rgb(0,0,1,alpha = 0.5))  
}

l <- 0.1
KXX <- sigma_f*DP_AddNugget(CF_ExponentialFamily(x , x , l , 2), nugget)
y <- f_x+BE_SampleGP(KXX)

p4 <- ggplot(data.frame(x , y) , aes(x , y)) + geom_line(col = rgb(0,0,1,alpha = 0.5))  + ggtitle('l = 0.1') + xlab('t') + ylab(TeX('$\\epsilon_{md}(t)$'))
for(i in 1:9){
  y <- f_x+BE_SampleGP(KXX)
  p4 <- p4 + geom_line(data = data.frame(x , y) , aes(x , y) , col = rgb(0,0,1,alpha = 0.5))  
}


x11()
grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2,top = textGrob("Samples from GPs",gp=gpar(fontsize=20,font=3)))
}
