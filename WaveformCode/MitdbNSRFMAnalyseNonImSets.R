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

SetOfNonImplausibleSets <- FM_GetNonImplausibleSetsFromlabelledDataset()

BC_PlotCompareTwoHists(SetOfNonImplausibleSets , SetOfNonImplausibleSets)
BC_PlotPairs(SetOfNonImplausibleSets , alpha = 0.01)

BC_PlotCompareTwoHists(PriorNonImplausibleSetRegularyIreRegular[1:1000,]  , unique(SetOfNonImplausibleSets)[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),] )
BC_PlotCompareTwoHists(PriorNonImplausibleSetRegular[1:1000,] , unique(SetOfNonImplausibleSets)[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),])

BC_PlotPairsFromTwoVariables(PriorNonImplausibleSetRegularyIreRegular[1:1000,]  , unique(SetOfNonImplausibleSets)[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),] , alpha = 0.05)
BC_PlotPairsFromTwoVariables(PriorNonImplausibleSetRegular[1:1000,] , unique(SetOfNonImplausibleSets)[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),] , alpha = 0.05)

