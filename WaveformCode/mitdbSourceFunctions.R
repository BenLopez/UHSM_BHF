mitdb_ComputePeakresults <- function(waveData, RWaveExtractedDataI , Annotation ){

  # Grab subset of annotation indicies
  AnnIn <- Annotation$Index[  ((Annotation$Code == 'N') 
                             | (Annotation$Code == 'A') 
                             | (Annotation$Code == 'Q') 
                             | (Annotation$Code == 'V') 
                             | (Annotation$Code == 'L')
                             | (Annotation$Code == 'R')
                             | (Annotation$Code == 'a')
                             | (Annotation$Code == 'J')
                             | (Annotation$Code == 'S')
                             | (Annotation$Code == 'r')
                             | (Annotation$Code == 'F')
                             | (Annotation$Code == 'e')
                             | (Annotation$Code == 'j')
                             | (Annotation$Code == 'n')
                             | (Annotation$Code == 'E')
                             | (Annotation$Code == '/')
                             | (Annotation$Code == 'f')
                             | (Annotation$Code == '?')
                             | (Annotation$Code == 'B')
                             | (Annotation$Code == 'I')) ]  
  
  AnnCode <- Annotation$Code[   ((Annotation$Code == 'N') 
                                 | (Annotation$Code == 'A') 
                                 | (Annotation$Code == 'Q') 
                                 | (Annotation$Code == 'V') 
                                 | (Annotation$Code == 'L')
                                 | (Annotation$Code == 'R')
                                 | (Annotation$Code == 'a')
                                 | (Annotation$Code == 'J')
                                 | (Annotation$Code == 'S')
                                 | (Annotation$Code == 'r')
                                 | (Annotation$Code == 'F')
                                 | (Annotation$Code == 'e')
                                 | (Annotation$Code == 'j')
                                 | (Annotation$Code == 'n')
                                 | (Annotation$Code == 'E')
                                 | (Annotation$Code == '/')
                                 | (Annotation$Code == 'f')
                                 | (Annotation$Code == '?')
                                 | (Annotation$Code == 'B')
                                 | (Annotation$Code == 'I')) ]  
  
  
  # Calulate positive values.
  P <- length( AnnIn )
  N <- length(waveData$Date) - P  
  
  distmatrix <- abs( as.matrix( pdist( as.matrix(as.numeric(RWaveExtractedDataI$t)) , as.matrix(as.numeric(waveData$Date[AnnIn])) ) ) )

  minvalue <- unlist(lapply(waveData$Date[AnnIn]  , function(X){min(abs(X - RWaveExtractedDataI$t ))}))
  
  #minvalue <- apply(distmatrix , 2 ,  min)
  #minindex <- apply(distmatrix , 2 ,  which.min)
  MissedPeaks <- AnnIn[minvalue > (0.1)] 
  MissedPeakCodes <- AnnCode[minvalue > (0.1)]
  
  minvalue <- unlist(lapply(RWaveExtractedDataI$t , function(X){min(abs(X - waveData$Date[AnnIn]  ))}))
  minindex <- unlist(lapply(RWaveExtractedDataI$t  , function(X){which.min(abs(X - waveData$Date[AnnIn]  ))}))

  ExtraPeaks <- which(minvalue > 0.1)
  
  FN <- length(MissedPeaks)
  FP <- length(ExtraPeaks)
  TP <- P - length(MissedPeaks)
  TN <- N - FP
  
    
  Sen  <-  TP / P 
  Spec <-  TN / N
  Acc  <-  (TP + TN)/(TP + FP +FN + TN)   
  PPV  <-  TP / (TP + FP)
  NPV  <-  TN / ( TN + FN)  
  
  
  BeatTypes <- DP_FindNumberUniques(AnnCode)
  SufficientStatistics <- setNames(list(P , N , TP , TN , FP , FN ) , c('P' , 'N' , 'TP' , 'TN' , 'FP' , 'FN' ))
  
  
  output <- setNames( list(Sen , Spec , Acc , PPV , NPV , MissedPeaks , ExtraPeaks , MissedPeakCodes , SufficientStatistics , BeatTypes)  ,
                      c('Sensitivity' , 'Specificity' , 'Accuracy' , 'PPV' , 'NPV' ,  'Missedpeaks' , 'Extrapeaks' , 'Missedpeakcodes' , 'SufficientStatistics' ,  'BeatTypes') )
  
  
  
  return(output)
}
mitdb_createlistofannotations <- function(){
  return(c('N'  ,
                                'L'  , 
                                'R'  , 
                                'B'  , 
                                'A'  , 
                                'a'  , 
                                'J'  , 
                                'S'  , 
                                'V'  , 
                                'r'  , 
                                'F'  , 
                                'e'  , 
                                'j'  , 
                                'r'  , 
                                'E'  , 
                                '/' , 
                                'f'  , 
                                'Q'  , 
                                '?' ))
}  
