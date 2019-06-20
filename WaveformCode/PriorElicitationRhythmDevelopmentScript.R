{
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
    set.seed( 1 )
  }
}
# SOme functions
FM_SampleGMMBigeminy <- function(X , N = 250 ){
  output <- matrix(0 , N , 1)
  for(i in 1:N){
    if(mod(i , 2) == 1){
      output[i , ] <- rnorm(1 , mean = X[4] , sd = (X[7]) )
    }  
    if(mod(i , 2) == 0){
      output[i , ] <- rnorm(1 , mean = X[6] , sd = (X[9]) )
    }  
  }
  return(output)
}



# Prior elicitation.


{X <- rep(0 , 10)

X[1] <- 0.5
X[2] <- 0
X[3] <- 0.5
X[4] <- 0.4
X[5] <- 0.0000001
X[6] <- 0.7
X[7] <- 0.01
X[8] <- 0.0000001
X[9] <- 0.01
X[10] <- 1}
  
x <- seq(0.5,1,0.001)  
f_x <- FM_EvaluateDenistyEstimate(x , X)

x11(20,20)
p1 <- ggplot(data.frame(RR = x , f_x = f_x) , aes(RR , f_x)) + geom_line(col = 'blue') +ggtitle('Distribution of RRTimes')
p1

RRTimes <- FM_SampleGMMBigeminy(X , 25000)

x11(20 , 4)
p2 <- ggplot(data.frame(t = cumsum(RRTimes) , RRTimes = RRTimes) , aes(t , RRTimes)) + geom_point(col =rgb(0,0,1,0.01)) + ylim(c(0,2)) +ggtitle('RRTimes')
p2 

RRTimes <- FM_SampleGMMBigeminy(X , 20)
t = cumsum(RRTimes)
t_observation = seq(0.25  , 10 , 0.005)

x11(20*0.984252 , 5*0.393701)
ECG <- PER_CreateECG( t , t_observation )
p3 <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
p3 <- p3 + theme(
  panel.background = element_rect(fill = rgb(1,0,0,alpha = 0.08), colour = "pink",
                                  size = 2, linetype = "solid"),
  panel.grid.major = element_line(size = 1, linetype = 'solid',
                                  colour = rgb(1,0,0,0.25)), 
  panel.grid.minor = element_line(size = 0.1, linetype = 'solid',
                                  colour = rgb(1,0,0,0.25))
) 

p3 <-p3 + scale_x_continuous(limits = c(0.4,9.6) , minor_breaks = seq(0, 10, 0.04)[-seq(1,251,5)] , breaks  = seq(0, 10, 0.2) ) + scale_y_continuous(minor_breaks = seq(-50, 200, 10) , breaks = seq(-50, 200, 50))
p3

x11(20,4)
ECG <- PER_CreateECG( t , t_observation )
p3 <- ggplot(data.frame(t = t_observation , V = ECG ) , aes(t , V)) + geom_line(col =rgb(0,0,0,0.9) , size = 0.7)
p3
