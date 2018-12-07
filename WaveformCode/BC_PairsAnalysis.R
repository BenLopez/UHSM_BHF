AFVev <- AdjustedBeliefs$W[ RPeaksStruct$RRCombined$t > DP_StripTime(MetaData$ConfirmedFirstNewAF) , ]
AFVev <- AdjustedBeliefs$W[ (RPeaksStruct$RRCombined$t >  AFLocations$Start[1])*((RPeaksStruct$RRCombined$t <  AFLocations$End[1])) == 1 , ]

x11(30,20)
SAMPLENUM <- dim(AFVev)[1]
variable = c(1:11)
tmp <- rbind(DataBase[[1]][sample(1:size(DataBase[[1]])[1] , SAMPLENUM ),variable] , 
             DataBase[[2]][sample(1:size(DataBase[[2]])[1] , SAMPLENUM ),variable],
             AFVev[sample(1:size(AFVev)[1] , SAMPLENUM ),variable])
al = 0.01
colvector <- c(rep(rgb(1 ,0 , 0, alpha = al)   , SAMPLENUM) , 
          rep(rgb(0 , 0 , 1, alpha = al)   , SAMPLENUM),
         rep(rgb(0 , 1 , 0, alpha = al)   , SAMPLENUM))
pairs(tmp , pch = 16 , col = colvector , labels = AFD_CreateDefaultSettings()$BinlimsScore[variable])
