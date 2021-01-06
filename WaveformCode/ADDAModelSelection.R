{
  NameofVariable <- NamesofPeropCategoricalVariables[1]
  PreOpCategoricalUnivariateResults <- POM_CategoricalUnivariateAnalysis(CategoricalData = as.factor(MasterData[ , which(names(MasterData) == NameofVariable)]) ,
                                                                         AFLogical = MasterData$AFLogical , MasterData , NameofVariable )
  for(i in 2:length(NamesofPeropCategoricalVariables)){
    NameofVariable <- NamesofPeropCategoricalVariables[i]
    PreOpCategoricalUnivariateResults <- rbind(PreOpCategoricalUnivariateResults,
                                               POM_CategoricalUnivariateAnalysis(CategoricalData = MasterData[ , which(names(MasterData) == NameofVariable)] ,
                                                                                 AFLogical = MasterData$AFLogical , MasterData , NameofVariable ))
  }
  
  #PreOpCategoricalUnivariateResults[is.na(PreOpCategoricalUnivariateResults)] <- 'NA'
  #PreOpCategoricalUnivariateResults[ which(is.na(PreOpCategoricalUnivariateResults[,9])) , 9] <- 'NA'
  #PreOpCategoricalUnivariateResults[ -which(PreOpCategoricalUnivariateResults[,9] < (1e-06)) , 9] <- 'NA' # no idea about -which() here
  #PreOpCategoricalUnivariateResults[ which(PreOpCategoricalUnivariateResults[,10] > (1e02)) , 10] <- 'NA'
  #PreOpCategoricalUnivariateResults[ which(PreOpCategoricalUnivariateResults[,8] > (1e02)) , 8] <- 'NA'
  
  PreOpCategoricalUnivariateResults <- FB_CorrectUnivariateTableNames(PreOpCategoricalUnivariateResults)
  PreOpCategoricalUnivariateResultsLatex <- FB_CorrectUnivariateTableValues(UnivariateTableLatex = xtable(PreOpCategoricalUnivariateResults) )
  
  # Pre_op Continuous Data Structure
  
  NameofVariable <- NamesofPreopContinuousVariables[1]
  ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
  PreOpContinuousUnivariateResults <- POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable)
  for(i in 2:length(NamesofPreopContinuousVariables)){
    NameofVariable <- NamesofPreopContinuousVariables[i]
    ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
    PreOpContinuousUnivariateResults <- rbind(PreOpContinuousUnivariateResults ,POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable) )
  }
  
  PreOpContinuousUnivariateResults <- FB_CorrectUnivariateTableNames(PreOpContinuousUnivariateResults)
  PreOpContinuousUnivariateResultsLatex <- FB_CorrectUnivariateTableValues(xtable(PreOpContinuousUnivariateResults))
  
  # Post_op Catagorical Data Structure
  NamesofPostopCategoricalVariables <- c("AKIUODayOne",
                                         'IABP2')

  NameofVariable <- NamesofPostopCategoricalVariables[1]
  PostOpCategoricalUnivariateResults <- POM_CategoricalUnivariateAnalysis(CategoricalData = MasterData[ , which(names(MasterData) == NameofVariable)] ,
                                                                          AFLogical = MasterData$AFLogical , MasterData , NameofVariable )
  
  for(i in 2:length(NamesofPostopCategoricalVariables)){
    NameofVariable <- NamesofPostopCategoricalVariables[i]
    PostOpCategoricalUnivariateResults <- rbind(PostOpCategoricalUnivariateResults,
                                               POM_CategoricalUnivariateAnalysis(CategoricalData = MasterData[ , which(names(MasterData) == NameofVariable)] ,
                                                                                 AFLogical = MasterData$AFLogical , MasterData , NameofVariable ))
  }
  
  
  #PostOpCategoricalUnivariateResults[is.na(PostOpCategoricalUnivariateResults)] <- 'NA'
  #PostOpCategoricalUnivariateResults[ which(is.na(PostOpCategoricalUnivariateResults[,9])) , 9] <- 'NA'
  #PostOpCategoricalUnivariateResults[ -which(PostOpCategoricalUnivariateResults[,9] < (1e-06)) , 9] <- 'NA' # no idea about -which() here
  #PostOpCategoricalUnivariateResults[ which(PostOpCategoricalUnivariateResults[,10] > (1e02)) , 10] <- 'NA'
  #PostOpCategoricalUnivariateResults[ which(PostOpCategoricalUnivariateResults[,8] > (1e02)) , 8] <- 'NA'
  
  PostOpCategoricalUnivariateResults <- FB_CorrectUnivariateTableNames(PostOpCategoricalUnivariateResults)
  PostOpCategoricalUnivariateResultsLatex <- FB_CorrectUnivariateTableNames(xtable(PostOpCategoricalUnivariateResults))
  
  
  # Post_op Continuous Data Structure
  
  
  NameofVariable <- NamesofPostopContinuousVariables[1]
  ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
  PostOpContinuousUnivariateResults <- POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable)
  for(i in 2:length(NamesofPostopContinuousVariables)){
    NameofVariable <- NamesofPostopContinuousVariables[i]
    ContinuousData <- MasterData[ ,which(names(MasterData) == NameofVariable)]
    PostOpContinuousUnivariateResults <- rbind(PostOpContinuousUnivariateResults ,POM_ContinuousUnivariateAnalysis(ContinuousData , AFLogical , MasterData , NameofVariable) )
  }
  
  #PostOpContinuousUnivariateResults[is.na(PostOpContinuousUnivariateResults)] <- 'NA'
  #PostOpContinuousUnivariateResults[ which(is.na(PostOpContinuousUnivariateResults[,9])) , 9] <- 'NA'
  #PostOpContinuousUnivariateResults[ -which(PostOpContinuousUnivariateResults[,9] < (1e-06)) , 9] <- 'NA' # no idea about -which() here
  #PostOpContinuousUnivariateResults[ which(PostOpContinuousUnivariateResults[,10] > (1e02)) , 10] <- 'NA'
  #PostOpContinuousUnivariateResults[ which(PostOpContinuousUnivariateResults[,8] > (1e02)) , 8] <- 'NA'
  
  PostOpContinuousUnivariateResults <- FB_CorrectUnivariateTableNames(PostOpContinuousUnivariateResults)
  PostOpContinuousUnivariateResultsLatex <- FB_CorrectUnivariateTableValues(xtable(PostOpContinuousUnivariateResults))
}


{
  print( PreOpCategoricalUnivariateResultsLatex )
  print( PreOpContinuousUnivariateResultsLatex )
  print( PostOpCategoricalUnivariateResultsLatex ) 
  print( PostOpContinuousUnivariateResultsLatex ) 
}

#### Group Models #####


#### Comparison of risk models with models built using underlying parameters ####

# LES

  {
    
    NamesofLogisticEuroScoreVariables <- POM_AddSquarestonames(NamesofLogisticEuroScoreVariables , ListNumericVariables)
    EuroComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofLogisticEuroScoreVariables ,MasterData = MasterData , ControlModel  = 'LogisticEUROScore')

    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = EuroComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    summary(model)  
    LogisticEuroscoreModelTableLatex <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
  }
  
  # SOFA
  
  {
    
    NamesofSofaScoreVariables <- POM_AddSquarestonames(NamesofSofaScoreVariables , ListNumericVariables)
    SOFAComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofSofaScoreVariables ,MasterData = MasterData , ControlModel  = 'SOFA')
    
    model <-  glm( formula = SOFAComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    summary(model)  
    SOFAEuroscoreModelTableLatex <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
  }
  
  # Log Casus
  
  {
    NamesofLogCasusScoreVariables <- POM_AddSquarestonames(NamesofLogCasusScoreVariables , ListNumericVariables)
    LogCasusComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesofLogCasusScoreVariables ,MasterData = MasterData , ControlModel  = 'logCASUS')
    
    model <-  glm( formula = LogCasusComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    summary(model)  
    CASUSEuroscoreModelTableLatex <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))

  }
  
  #### Group Analysis #####
  
  {
    
    NamesPreOpBleedingVariables <- POM_AddSquarestonames(NamesPreOpBleedingVariables , ListNumericVariables)
    PreOpBleedingComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesPreOpBleedingVariables ,MasterData = MasterData )
    
    model <-  glm( formula = PreOpBleedingComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpBleedingModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopBleedingModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
  
    NamesPreOpCardioVascularVariables <- POM_AddSquarestonames(NamesPreOpCardioVascularVariables , ListNumericVariables)
    PreOpCardioComparisonOutputs <- POM_GroupComparison(NamesofVariables = NamesPreOpCardioVascularVariables ,MasterData = MasterData )
    
    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = PreOpCardioComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpCardioVascularModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopCardioVascularModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
    
    NamesPreOpDemographicVariables <- POM_AddSquarestonames(NamesPreOpDemographicVariables , ListNumericVariables)
    PreOpDemographicComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpDemographicVariables ,MasterData = MasterData )
    
    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = PreOpDemographicComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpDemographicModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopDemographicModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
    
    NamesPreOpElectrolytesVariables <- POM_AddSquarestonames(NamesPreOpElectrolytesVariables , ListNumericVariables)
    PreOpElectrolytesComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpElectrolytesVariables ,MasterData = MasterData )
    
    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = PreOpElectrolytesComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpElectrolytesModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopElectrolytesModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
 
    NamesPreOpInflamatoryVariables <- POM_AddSquarestonames(NamesPreOpInflamatoryVariables , ListNumericVariables)
    PreOpInflamatoryComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpInflamatoryVariables ,MasterData = MasterData )
    
    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = PreOpInflamatoryComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpInflamatoryModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopInflamatoryModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
    NamesPreOpLiverVariables <- POM_AddSquarestonames(NamesPreOpLiverVariables , ListNumericVariables)
    PreOpLiverComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpLiverVariables ,MasterData = MasterData )
    
    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = PreOpLiverComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpLiverModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopLiverModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
    
    NamesPreOpOperativeVariables <- POM_AddSquarestonames(NamesPreOpOperativeVariables , ListNumericVariables)
    PreOpOperativeComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpOperativeVariables ,MasterData = MasterData )
    
    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = PreOpOperativeComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpOperativeModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopOperativeModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
    NamesPreOpRenalVariables <- POM_AddSquarestonames(NamesPreOpRenalVariables , ListNumericVariables)
    PreOpRenalComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpRenalVariables ,MasterData = MasterData )
    
    DataForLogistic <- POM_SampledImputation(MasterData)
    model <-  glm( formula = PreOpRenalComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
    MasterData$PreOpRenalModel <- predict(model , DataForLogistic ,  type = c('response'))
    summary(model)  
    PreopRenalModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
    
    
    NamesPreOpRespiratoryVariables <- c('FiO26num'
    )
    #PreOpRespiratoryComparisonOutputs  <- POM_GroupComparison(NamesofVariables = NamesPreOpRespiratoryVariables ,MasterData = MasterData )
    #}
    
    #### Post Op Groups
    {
      {
       
        NamesPostOpBleedingVariables <- POM_AddSquarestonames(NamesPostOpBleedingVariables , ListNumericVariables)
        PostOpBleedingComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpBleedingVariables,NamesPreOpBleedingVariables) ,MasterData = MasterData )
        
        DataForLogistic <- POM_SampledImputation(MasterData)
        model <-  glm( formula = PostOpBleedingComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
        MasterData$PostOpBleedingModel <- predict(model , DataForLogistic ,  type = c('response'))
        summary(model)  
        PostopBleedingModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
        
      }
      {
       
        NamesPostOpCardioVascularVariables <- POM_AddSquarestonames(NamesPostOpCardioVascularVariables , ListNumericVariables)
        PostOpCardioVascularComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpCardioVascularVariables,NamesPreOpCardioVascularVariables) ,MasterData = MasterData )
        
        DataForLogistic <- POM_SampledImputation(MasterData)
        model <-  glm( formula = PostOpCardioVascularComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
        MasterData$PostOpCardioVascularModel <- predict(model , DataForLogistic ,  type = c('response'))
        summary(model)  
        PostopCardioVascularModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
        
        {
         
          NamesPostOpElectrolytesVariables <- POM_AddSquarestonames(NamesPostOpElectrolytesVariables , ListNumericVariables)
          PostOpElectrolytesComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpElectrolytesVariables,NamesPreOpElectrolytesVariables) ,MasterData = MasterData )
          
          DataForLogistic <- POM_SampledImputation(MasterData)
          model <-  glm( formula = PostOpElectrolytesComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
          MasterData$PostOpElectrolytesModel <- predict(model , DataForLogistic ,  type = c('response'))
          summary(model)  
          PostopElectrolytesModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
          
          NamesPostOpInflamatoryVariables <- POM_AddSquarestonames(NamesPostOpInflamatoryVariables , ListNumericVariables)
          PostOpInflamatoryComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpInflamatoryVariables,NamesPreOpInflamatoryVariables) ,MasterData = MasterData )
          
          DataForLogistic <- POM_SampledImputation(MasterData)
          model <-  glm( formula = PostOpInflamatoryComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
          MasterData$PostOpInflamatoryModel <- predict(model , DataForLogistic ,  type = c('response'))
          summary(model)  
          PostopInflamatoryModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
          
          
          NamesPostOpLiverVariables <- POM_AddSquarestonames(NamesPostOpLiverVariables , ListNumericVariables)
          PostOpLiverComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpLiverVariables,NamesPreOpLiverVariables) ,MasterData = MasterData )
          DataForLogistic <- POM_SampledImputation(MasterData)
          model <-  glm( formula = PostOpLiverComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
          MasterData$PostOpLiverModel <- predict(model , DataForLogistic ,  type = c('response'))
          summary(model)  
          PostopLiverModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
          
      
          #NamesPostOpNeuroVariables <- c('GCSnum'
          #                               )
          #PostOpNeuroComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpNeuroVariables) ,MasterData = MasterData )
          
          
         
          NamesPostOpOperativeVariables <- POM_AddSquarestonames(NamesPostOpOperativeVariables , ListNumericVariables)
          PostOpOperativeComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpOperativeVariables,NamesPreOpOperativeVariables) ,MasterData = MasterData )
          
          DataForLogistic <- POM_SampledImputation(MasterData)
          model <-  glm( formula = PostOpOperativeComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
          MasterData$PostOpOperativeModel <- predict(model , DataForLogistic ,  type = c('response'))
          summary(model)  
          PostopOperativeModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
      
          
          NamesPostOpRenalVariables <- POM_AddSquarestonames(NamesPostOpRenalVariables , ListNumericVariables)
          PostOpRenalComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpRenalVariables,NamesPreOpRenalVariables) ,MasterData = MasterData )
          
          DataForLogistic <- POM_SampledImputation(MasterData)
          model <-  glm( formula = PostOpRenalComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
          MasterData$PostOpRenalModel <- predict(model , DataForLogistic ,  type = c('response'))
          summary(model)  
          PostopRenalModelTable <- FB_CorrectModelTableNames(modeltable = FC_AlterModelSummary(model,MasterData ))
          
          
          NamesPostOpRespiratoryVariables <- POM_AddSquarestonames(NamesPostOpRespiratoryVariables , ListNumericVariables)
          PostOpRespiratoryComparisonOutputs <- POM_GroupComparison(NamesofVariables = c(NamesPostOpRespiratoryVariables,NamesPreOpRespiratoryVariables) ,MasterData = MasterData )
          
          DataForLogistic <- POM_SampledImputation(MasterData)
          model <-  glm( formula = PostOpRespiratoryComparisonOutputs$ModelAUC , family=binomial(link='logit') , data=DataForLogistic)
          MasterData$PostOpRespiratoryModel <- predict(model , DataForLogistic ,  type = c('response'))
          summary(model)  
          PostopRespiratoryModelTable <- FB_CorrectModelTableNames( modeltable = FC_AlterModelSummary(model,MasterData ) )
          
          
        }
      }
    }
  }   
  #### Pre_op Predictive Model #####
  
  {
  NamesofModels <- names(MasterData)[which(names(MasterData) == 'PreOpBleedingModel'):(which(names(MasterData) == 'PreOpBleedingModel') + 7)]
  NamesofPeropVariables <- c(NamesofPeropCategoricalVariables , POM_AddSquarestonames(c(NamesofModels, NamesofPreopContinuousVariables[c(-2,-4)]) , ListNumericVariables) )
  IndiciesofPeropVariables <- which(names(MasterData) %in% NamesofPeropVariables)
  NamesofPeropVariables <- names(MasterData)[IndiciesofPeropVariables]
  
  { matrixforstep <- diag(length(IndiciesofPeropVariables))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAUC))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAUC))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAUC))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAUC))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAUC))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAUC))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAUC))
    matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = EuroComparisonOutputs$ModelAUC))
    
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = EuroComparisonOutputs$ModelAIC))
  }
  
  PreOpStepOutput <- FC_StepWiseForwardBackwardAUC(PreoperativeIndices = IndiciesofPeropVariables , MasterData , matrixforstep)
  
  
  model <-  glm( formula = PreOpStepOutput[[2]] , family=binomial(link='logit') , data=DataForLogistic)
  MasterData$PreOpModel <- predict(model , DataForLogistic ,  type = c('response'))
  
  AUC2 <- FC_CalculateCrossValidatedROC(formulaformodel = PreOpStepOutput[[2]] , PreoperativeIndices = IndiciesofPeropVariables , MasterData)
  
  #### Post_op Predictive Model #####
  
  NamesofModels <- names(MasterData)[which(names(MasterData) == 'PreOpBleedingModel'):length(names(MasterData))]
  NamesofPostopVariables <- c(NamesofPostopCategoricalVariables ,POM_AddSquarestonames(c(NamesofModels,NamesofPostopContinuousVariables) , ListNumericVariables) )
  NamesofPostopVariables <- c(NamesofPostopVariables , NamesofPeropVariables)
  IndiciesofPostopVariables <- which(names(MasterData) %in% NamesofPostopVariables)
  NamesofPostopVariables <- names(MasterData)[IndiciesofPostopVariables]
  
  { matrixforstep <- diag(length(IndiciesofPostopVariables))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = EuroComparisonOutputs$ModelAUC))

    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpBleedingComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpDemographicComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpCardioComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpElectrolytesComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpInflamatoryComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpLiverComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = PreOpOperativeComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind(matrixforstep ,  FC_ExtractIndiciesFromFormula(NamesofVariables = NamesofPeropVariables , ModelFormula = EuroComparisonOutputs$ModelAIC))
    
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = SOFAComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = LogCasusComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpCardioVascularComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpElectrolytesComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpBleedingComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpInflamatoryComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpLiverComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpOperativeComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRenalComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRespiratoryComparisonOutputs$ModelAUC))
    matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PreOpStepOutput[[2]]))
    
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = SOFAComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = LogCasusComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpCardioVascularComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpElectrolytesComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpBleedingComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpInflamatoryComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpLiverComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpOperativeComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRenalComparisonOutputs$ModelAIC))
    #matrixforstep <- rbind( matrixforstep ,  FC_ExtractIndiciesFromFormula( NamesofVariables = NamesofPostopVariables , ModelFormula = PostOpRespiratoryComparisonOutputs$ModelAIC))
    
  }
  #PostOpStepOutput <- FC_StepWiseForwardAUC(PreoperativeIndices = IndiciesofPostopVariables , MasterData , matrixforstep)
  
  PostOpStepOutput <- FC_StepWiseForwardBackwardAUC(PreoperativeIndices = IndiciesofPostopVariables , MasterData , matrixforstep)
  }
  
  ##### Model validation
  {
    {
      DataForLogistic <- POM_SampledImputation(MasterData)
      model <- glm(formula = PreOpStepOutput[[2]]
                   , family = binomial(link = "logit"),  data = DataForLogistic)
      summary(model)

      LogisticProbility <- predict(model , DataForLogistic , type = c('response'))
      auc( MasterData$AFLogical , c(LogisticProbility) )
      ci.auc( MasterData$AFLogical , c(LogisticProbility) )
      
      PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(as.matrix(LogisticProbility) ,MasterData$AFLogical == 1))
      ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1) + xlab('1 - Specificity')
      NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-sum(AFLogical)/length(AFLogical))) + geom_hline(yintercept = sum(AFLogical)/length(AFLogical))+  ggtitle('PPV vs NPV Curves')
      PPVSenPlot <- BC_PlotPPVSen(PerformanceSweep3) +ggtitle('Precision Recall Curves')
      
      LogisticProbility <- FB_CalulateCrossValidatedProbilities(PreOpStepOutput[[2]],PreOpStepOutput[[1]], MasterData)
      auc( MasterData$AFLogical , c(LogisticProbility) )
      ci.auc( MasterData$AFLogical , c(LogisticProbility) )
      
      print(BC_CalculateAreaUnderCurve2(PerformanceSweep3))
      
      phist <- ggplot(data.frame( p =LogisticProbility[MasterData$AFLogical ==0]) , aes(p)) + geom_density(col = rgb(0,0,1,alpha = 0.5)) +
               geom_density(data = data.frame( p =LogisticProbility[MasterData$AFLogical ==1]) , aes(p), col = rgb(1,0,0,alpha = 0.5)) + 
               ggtitle('Histogram of Model Output')
      
      EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$AFLogical == 1), BinWidth = 0.1))
      EmpericalProbabilityPlot <- BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities')
      
      x11(20,20)
      grid.arrange(phist , EmpericalProbabilityPlot , ROCplot , NPVPPVPlot  , ncol =2 , nrow = 2)
      x11(20,20)
      PPVSenPlot
      }
    Preopfullmodeltable <- FB_CorrectModelTableNames( modeltable = FC_AlterModelSummary( model , MasterData ) )
    
    {
      model <- glm(formula = PostOpStepOutput[[2]]
                   , family = binomial(link = "logit"),  data = DataForLogistic)
      summary(model)

      LogisticProbility <- predict(model , DataForLogistic , type = c('response'))
      #LogisticProbility <- FB_CalulateCrossValidatedProbilities(PostOpStepOutput[[2]],PostOpStepOutput[[1]], MasterData)
      
      PerformanceSweep3 <- BC_PerformanceSweep(GlobalProbCalibrationStruct = cbind(as.matrix(LogisticProbility) ,MasterData$AFLogical == 1))
      ROCplot <- BC_PlotsCreateROC(PerformanceSweep3) + ggtitle('ROC Curves') + geom_abline(intercept = 0 , slope = 1) + xlab('1 - Specificity')
      NPVPPVPlot <- BC_PlotsCreateNPVPPV(PerformanceSweep3) + geom_vline(xintercept =( 1-(sum(AFLogical)/length(AFLogical)) )) + geom_hline(yintercept = (sum(AFLogical)/length(AFLogical)))+  ggtitle('PPV vs NPV Curves')
      PPVSenPlot <- BC_PlotPPVSen(PerformanceSweep3) + ggtitle('Precision Recall Curves (%)') +xlab('Precision (%)') +ylab('Recall/Sensitivity')
      
      auc( MasterData$AFLogical , c(LogisticProbility) )
      ci.auc( MasterData$AFLogical , c(LogisticProbility) )
      
      LogisticProbility <- FB_CalulateCrossValidatedProbilities(PostOpStepOutput[[2]],PostOpStepOutput[[1]], MasterData)
      auc( MasterData$AFLogical , c(LogisticProbility) )
      ci.auc( MasterData$AFLogical , c(LogisticProbility) )
         
      print(BC_CalculateAreaUnderCurve2(PerformanceSweep3))
      
      phist <- ggplot(data.frame( p =LogisticProbility[MasterData$AFLogical ==0]) , aes(p)) + geom_density(col = rgb(0,0,1,alpha = 0.5)) +
               geom_density(data = data.frame( p =LogisticProbility[MasterData$AFLogical ==1]) , aes(p), col = rgb(1,0,0,alpha = 0.5)) + 
               ggtitle('Histogram of Model Output')
      
      EmpericalProbabilityStructure <- BC_CleanProbCalibrationOutput(BC_CreateCalibrationStructure(cbind(LogisticProbility , MasterData$AFLogical == 1), BinWidth = 0.1))
      EmpericalProbabilityPlot <- BC_PlotCreateProbabilityCalibrationPlot(EmpericalProbabilityStructure) + ggtitle('Emperical Probabilities')
      
      x11(20,20)
      grid.arrange(phist , EmpericalProbabilityPlot , ROCplot , NPVPPVPlot  , ncol =2 , nrow = 2)
      }
    Postopfullmodeltable <-  FB_CorrectModelTableNames( modeltable = FC_AlterModelSummary(model,MasterData ) )
    
  }
  

index20 <- which.min(abs(PerformanceSweep3[,3] - 0.2))
index33 <- which.min(abs(PerformanceSweep3[,3] - 0.33))
index50 <- which.min(abs(PerformanceSweep3[,3] - 0.5))

x11(20,20)
grid.arrange(ROCplot ,PPVSenPlot   , ncol=2 , nrow =1)
