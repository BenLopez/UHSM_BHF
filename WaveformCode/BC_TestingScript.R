
{Patinettotest <- DP_choosepatient(listAllPatients = listAllPatients)
BCParameters$TS_Likelihood_clique <- 200
source('BC_LoadDataandTestSinglePatient.R')
if(nrow(AFLocations)>0){
  timetoview <- AFLocations$Start[1]
  }else{
  timetoview <- DP_SelectTimetoview(ECGs$ECGI$Date)
  }
lengthtoview <- 10
source('BC_CreatePlots.R')

x11(15,12)
print(grid.arrange( ECGIPlot,
                    ECGIIPlot,
                    ECGIIIPlot,
                    RAPlot,
                    RRPlot,
                    Inferenceplot,
                    nrow = 6 ,
                    ncol = 1  , 
                    top = paste0('Sensitivity= ' , Performance$Sensitvity , ' Specifictity= ' , Performance$Specifictity ,' PPV= ' ,Performance$NPV , ' NPV= ' , Performance$PPV ) ) )
}

#AFVev <- AdjustedBeliefs$W[ RPeaksStruct$RRCombined$t > DP_StripTime(MetaData$ConfirmedFirstNewAF) , ]
AFVev <- AdjustedBeliefs$W[ (RPeaksStruct$RRCombined$t >  AFLocations$Start[1])*((RPeaksStruct$RRCombined$t <  AFLocations$End[1])) == 1 , ]

x11(30,20)
SAMPLENUM <- dim(AFVev)[1]
variable = c(1:11)
tmp <- rbind(DataBase[[1]][sample(1:size(DataBase[[1]])[1] , SAMPLENUM ),variable] , 
             DataBase[[2]][sample(1:size(DataBase[[2]])[1] , SAMPLENUM ),variable],
             AFVev[sample(1:size(AFVev)[1] , SAMPLENUM ),variable])
al = 0.01
colvector <- c(rep(rgb(0 ,0 , 1, alpha = al)   , SAMPLENUM) , 
               rep(rgb(1 , 0 , 0, alpha = al)   , SAMPLENUM),
               rep(rgb(0 , 1 , 0, alpha = al)   , SAMPLENUM))
pairs(tmp , pch = 16 , col = colvector , labels = AFD_CreateDefaultSettings()$BinlimsScore[variable])

