
for(Patinettotest in unlist(DataBaseMaster$AFPatinetsNames)){ 
source('BC_LoadDataandTestSinglePatient.R')
if(nrow(AFLocations) == 0){
timetoview <- RPeaksStruct$RRCombined$t[1000]}else{
timetoview <- AFLocations$Start[1]}
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

#x11()
#plot(RPeaksStruct$RRCombined$t, rollapply(RPeaksStruct$RRCombined$RR , width = 250 , na.pad = TRUE , FUN = function(X){IQR(X, na.rm= TRUE)} ) , type ='l'  )
#title(Patinettotest)
}

for(Patinettotest in unlist(DataBaseMaster$AFPatinetsNames)){ 
  source('BC_LoadDataandTestSinglePatient.R')
  if(nrow(AFLocations) == 0){
    print(Patinettotest)
    timetoview <- RPeaksStruct$RRCombined$t[10000] }else{
 timetoview <- AFLocations$Start[1]}
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
