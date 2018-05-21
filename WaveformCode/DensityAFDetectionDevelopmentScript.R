pathFiles <- choose.dir(caption="Select folder with source code.")
pathFiles <- paste0(pathFiles, "\\")
setwd(pathFiles)

source("LibrariesAndSettings.R" , print.eval  = TRUE )

load(choose.files(caption = "Select patientindexmaster.RData"))
DP_ChooseDataReps()

DP_choosepatient(listAllPatients)
sub_pat = subset( PatIndex2017, PseudoId %in% subList )

print('Loading ECG.')
WaveData <- DP_LoadECG(path , subList , numberrep , ECGNum = 1)
print('ECG Loaded.')

interestingtimepoint <- DP_SelectInterestingTimePoint(WaveData[ seq(from = 1 , to = length(WaveData[ , 1]), by = 1000) , ] , sub_pat)
interestingindex <- which.min( abs(WaveData[ , 1] - as.POSIXct(interestingtimepoint)))
WaveData <- DP_CropWaveData(WaveData , interestingindex , DP_SelectHoursBeforeandAfter() )

print('Calulating Rpeaks')
RWaveExtractedData <- CleanRpeaks(RPeakExtractionWavelet( WaveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8) , 2)
print('Rpeaks Calulated')


AFScore <- ExtractIHVAFScore(RWaveExtractedData ,  binlims <- c(0, seq(from = 0.25  , to = 1.8  , 0.05  ) , 3))
VBM <- AFScore$IHAVFScore

p <- ggplot(RWaveExtractedData , aes(t , RR)) +
  geom_point(colour="blue", alpha=0.01) +
  xlab("t") +
  ylab("RR") + coord_cartesian(ylim = c(0, 1.2)) +
  geom_vline( xintercept = as.numeric( as.POSIXct( sub_pat$FirstNewAF) ) , linetype="dashed" , color = "purple" )+
  ggtitle( paste0(subList , " RR-times" ) ) 
p2 <- ggplot(data.frame(x=t , y=VBM), aes(x,y)) + geom_line(colour = "blue") + ggtitle('AF Score') + ylim(c(0,600))

tmp <- ASWF_GetStartEndAF(RWaveExtractedData$t , logicaltimeseries = (VBM > 150) )

for( i in ( 1:length(tmp$Start) ) )
{
p <- p + annotate("rect" , xmin = tmp$Start[i], xmax = tmp$End[i], ymin = -1000, ymax= 1000 , fill = 'pink' , alpha = 0.5)
}

print(grid.arrange( p  + 
                      geom_vline( xintercept = as.numeric( as.POSIXct(sub_pat$FirstNewAF[1] ) ) ,linetype="dashed" , color = "black" ) ,
                    p2 + geom_hline( yintercept = 150 ,linetype="dashed" , color = "black" ) , nrow = 2 , ncol = 1))


par(mforw = c(1 , 1))
plot(binmatrix[ , 1] , type ='l' , ylim = c(0,1))
for(i in 1:dim(binmatrix)[2])
{
  lines(binmatrix[ , i] , col = i)
}


timeintervalforhist <- 5 # time interval in minutes
starttime <- RWaveExtractedData[1 , 1]


timelogical <- ( (RWaveExtractedData[ , 1] > starttime)*(RWaveExtractedData[ , 1] < (starttime + (timeintervalforhist*60)) )*(RWaveExtractedData$RR < 2) ) == 1
p2 <- ggplot(RWaveExtractedData[ timelogical , ] , aes( RR)) +
    geom_histogram(colour="blue")

#dev.off()
grid.arrange( p  + 
                geom_vline( xintercept = as.numeric( as.POSIXct(starttime ) ) ,linetype="dashed" , color = "black" ) + 
                geom_vline( xintercept = as.numeric( as.POSIXct(starttime + (timeintervalforhist*60) )) ,linetype="dashed" , color = "black" ) ,
               p2 , nrow = 2 , ncol = 1)

starttime <- starttime - (timeintervalforhist*60)

Beta <- cumsum(as.numeric(RWaveExtractedData$RR))/cumsum(rep(1,length(as.numeric(RWaveExtractedData$RR))))
W_t <- smth(RWaveExtractedData$RR - Beta , method = 'sma' , n = 500 )

tmp <- as.numeric(RWaveExtractedData$RR - Beta - W_t)
t <- RWaveExtractedData$t[!is.na(tmp)]
tmp <- tmp[!is.na(tmp)]
#tmp <- tmp[tmp < quantile(tmp , probs = 0.999 , na.rm = 'TRUE') ]
#tmp <- tmp[tmp > quantile(tmp , probs = 0.001 , na.rm = 'TRUE') ]

mu <- smth((tmp) , method = 'sma' , n = 100 )
V_t <- smth(as.numeric(tmp)^2 , method = 'sma' , n = 100 ) - mu^2
S_t <- (smth(as.numeric(tmp)^3 , method = 'sma' , n = 100 ) - 3*mu*V_t)/(V_t^(1.5))
S_t <- S_t[!is.na(S_t)]
K_t <- smth( (tmp)^4 , method = 'sma' , n = 100 )/(V_t^2)
K_t <- K_t[!is.na(K_t)]
mu <- mu[!is.na(K_t)]

t <- t[!is.na(V_t)]
V_t <- V_t[!is.na(V_t)]
Beta_V <- cumsum(V_t)/cumsum(rep(1,length(as.numeric(V_t))))

#H <- as.matrix(cbind(rep(1 , length(t)) , 1:length(t)))
#Beta_V <- solve(t(H)%*%H)%*%(t(H)%*%as.matrix(V_t))


#Sigma_W <- cumsum(V_t[!is.na(V_t)] - Beta_V)/
W_V_t <- smth(V_t[!is.na(V_t)] - Beta_V , method = 'sma' , n = 100 )
t <- t[!is.na(W_V_t)]
W_V_t <- W_V_t[!is.na(W_V_t)]

Sigma_W <- var(W_V_t[ difftime(t , t[1] , units = 'secs') < 30*60])
stdresid <- abs(W_V_t/sqrt(Sigma_W))
  
p3 <- ggplot(data.frame( x = t , y = stdresid ) , aes(x,y) ) + geom_line( colour = "blue")
grid.arrange( p  + 
                geom_vline( xintercept = as.numeric( as.POSIXct(sub_pat$FirstNewAF[1] ) ) ,linetype="dashed" , color = "black" ) ,
              p3 + geom_hline(yintercept =3 ), nrow = 2 , ncol = 1)




