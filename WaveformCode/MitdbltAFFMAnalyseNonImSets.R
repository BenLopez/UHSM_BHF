

DP_ChooseDataReps()


i <- 3
AFNonImplausibleSets <- matrix(0 , 0 , 10)
#length(listAllPatients)
for(i in 1:length(listAllPatients)){

load(paste0("D:\\ltafdb_outputs\\HMOutput" , listAllPatients[[i]] , '.RData'))
load(paste0(path , '\\' , listAllPatients[[i]] , '\\Zip_out\\' ,  listAllPatients[[i]] , '_', "AFTimes" , '.RData'))
RPeakData <- DP_LoadRpeaksfile( path ,  listAllPatients[[i]]  )

AFTimes <-  FM_ltafdbExtractStartandEndAF(RPeakData , AFlocations , minutethreshold = 1)
AFLogical <- BC_CreateAnnotationFromInference(RPeakData$RRCombined$t ,AFTimes)

MetaData = DP_CreateDummyMetaData(PatIndex2017 = PatIndex2017 , Name = listAllPatients[[i]] )  
RRPlot <- BC_PlotCreateRRTimesPlots(RPeaksStruct = RPeakData , MetaData = MetaData) + ggtitle( paste0( listAllPatients[[i]]  , ' RRTimes' ) )
RRPlot <- BC_plotAddColouredRegions(plotstruct = RRPlot , Locations = AFTimes , fillcolor = 'red')
x11()
print(RRPlot)

#for(j in 1:length(outputstruct)){
#  if(FM_ltafCheckSegmentisAF(outputstruct[[j]][[4]] ,  AFLogical , RPeakData , n = 500)){
#    if(!is.null(outputstruct[[j]][[2]]$NonImplausibleSets)){
#    AFNonImplausibleSets <- unique( rbind(AFNonImplausibleSets,outputstruct[[j]][[2]]$NonImplausibleSets ))    
#    }
#    }
#}
}

BC_PlotPairs(AFNonImplausibleSets[sample(1:dim(AFNonImplausibleSets)[1] , 1000) , ] , alpha = 0.05)
BC_PlotPairsFromTwoVariables(AFNonImplausibleSets[sample(1:dim(AFNonImplausibleSets)[1] , 1000) , ] , PriorNonImplausibleSetRegular[sample(100000 :dim(PriorNonImplausibleSetRegular)[1] , 1000) ,]  , alpha = 0.05)
BC_PlotPairsFromThreeVariables(AFNonImplausibleSets[sample(1:dim(AFNonImplausibleSets)[1] , 1000) , ] , PriorNonImplausibleSetRegularyIreRegular[1:1000,] , rbind(PriorNonImplausibleSetRegular[sample(1:dim(PriorNonImplausibleSetRegular)[1] , 500) ,] ,PriorNonImplausibleSetRegular[sample(100000 :dim(PriorNonImplausibleSetRegular)[1] , 500) ,] )  , alpha = 0.05)

