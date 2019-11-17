
pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))

load("MasterData.RData")

source("PreOpModelSourceFunctions.R")
source('StepWiseModelSelectionSourceFunctions.R')

outputDischargeStatus <- POM_CategoricalUnivariateAnalysis(CategoricalData = as.factor(MasterData[ , which(names(MasterData) == 'DischargeStatus')]),
                                  AFLogical = MasterData$AFLogical , MasterData , 'DischargeStatus' )  

outputICUStay <- POM_ContinuousUnivariateAnalysis(ContinuousData = as.factor(MasterData[ , which(names(MasterData) == 'CICUStay')]),
                                            AFLogical = MasterData$AFLogical , MasterData , 'CICUStay' )  
