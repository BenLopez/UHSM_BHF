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
DP_ChooseDataReps()


SetOfNonImplausibleSets <- matrix(0 , 0 , 10)

for(j in 1:length(listAllPatients)){
  load(paste0("D:\\nsrdb_outputs\\HMOutput" , listAllPatients[[j]] , '.RData'))
for( i in 1:length(outputstruct)){
  SetOfNonImplausibleSets <-  rbind(SetOfNonImplausibleSets , outputstruct[[i]][[2]]$NonImplausibleSets )
}
}

BC_PlotPairs(unique(SetOfNonImplausibleSets[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),]))
BC_PlotPairsFromThreeVariables(PriorNonImplausibleSetRegularyIreRegular[1:1000,] , PriorNonImplausibleSetRegular[1:1000,] , unique(SetOfNonImplausibleSets)[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),] , alpha = 0.1)  

#BC_PlotCompareTwoHists(PriorNonImplausibleSetRegularyIreRegular[1:1000,]  , unique(SetOfNonImplausibleSets)[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),] )
#BC_PlotCompareTwoHists(PriorNonImplausibleSetRegular[1:1000,] , unique(SetOfNonImplausibleSets)[sample(1:dim(unique(SetOfNonImplausibleSets))[1] , 1000),])


             