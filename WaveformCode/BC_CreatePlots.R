
RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeaksStruct , MetaData = MetaData )
RAPlot <- BC_PlotCreateRATimesPlots(RPeaksStruct = RPeaksStruct , ECGs = ECGs ,  MetaData = MetaData )
ECGIPlot <- BC_PlotCreateECGPlots(RPeaksStruct ,  ECGs$ECGI  , timestart = timetoview , ECGindex = 1, timeindex = lengthtoview)
ECGIPlot <- BC_PlotECGAddTitle(ECGIPlot , MetaData , timetoview)
ECGIIPlot <- BC_PlotCreateECGPlots(RPeaksStruct ,  ECGs$ECGII  , timestart = timetoview , ECGindex = 2 , color = 'red' , timeindex = lengthtoview)
ECGIIPlot <- BC_PlotECGAddTitle(ECGIIPlot , ECGindex = 2)
ECGIIPlot <- BC_PlotECGAlignAxis(ECGIIPlot , c(timetoview, timetoview + lengthtoview))
ECGIIIPlot <- BC_PlotCreateECGPlots(RPeaksStruct ,  ECGs$ECGIII  , timestart = timetoview , ECGindex = 3 , color = 'green', timeindex = lengthtoview)
ECGIIIPlot <- BC_PlotECGAddTitle(ECGIIIPlot , ECGindex = 3)
ECGIIIPlot <- BC_PlotECGAlignAxis(ECGIIIPlot , c(timetoview , timetoview + lengthtoview))
RAPlot <- BC_PlotAddViewingRegionLines(RAPlot , c(timetoview, timetoview + lengthtoview))
RRPlot <- BC_PlotAddViewingRegionLines(RRPlot , c(timetoview, timetoview + lengthtoview))
Inferenceplot <- BC_PlotCreatePosteriorBeliefPlot(t = RPeaksStruct$RRCombined$t , Implausability , PosteriorProbabilities)
RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = AFLocations , fillcolor = 'pink')
RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = BadDataLocations , fillcolor = 'Black')


x11(15,12)
print(grid.arrange( ECGIPlot,
                    ECGIIPlot,
                    ECGIIIPlot,
                    RRPlot,
                    Inferenceplot,
                    nrow = 5 ,
                    ncol = 1  , 
                    top = paste0('Sensitivity= ' , Performance$Sensitvity , ' Specifictity= ' , Performance$Specifictity ,' PPV= ' ,Performance$NPV , ' NPV= ' , Performance$PPV ) ) )

UserResponse <- winDialog(type = c('yesnocancel') , message = 'Would you like to view another time period?')

if(UserResponse == 'CANCEL')
{
}
if(UserResponse == 'YES')
{
  timetoview <- BC_TimetoViewChange(timetoview)
  dev.off()
  source('BC_CreatePlots.R')
}
