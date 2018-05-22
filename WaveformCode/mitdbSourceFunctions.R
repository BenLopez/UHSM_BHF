mitdb_ComputePeakresults <- function( WaveData , Annotation )
{
  
  RWaveExtractedDataI <- RPeakExtractionWavelet( waveData , wt.filter(filter = "d6" , modwt=TRUE, level=1) , nlevels = 12 , ComponetsToKeep = c(3,4) , stdthresh = 2.8)
  
  # Grab subset of annotation indicies
  AnnIn <- Annotation$Index[  ((Annotation$Code == 'N') 
                             + (Annotation$Code == 'A') 
                             + (Annotation$Code == 'Q') 
                             + (Annotation$Code == 'V') 
                             + (Annotation$Code == 'L')
                             + (Annotation$Code == 'R')
                             + (Annotation$Code == 'a')
                             + (Annotation$Code == 'J')
                             + (Annotation$Code == 'S')
                             + (Annotation$Code == 'r')
                             + (Annotation$Code == 'F')
                             + (Annotation$Code == 'e')
                             + (Annotation$Code == 'j')
                             + (Annotation$Code == 'n')
                             + (Annotation$Code == 'E')
                             + (Annotation$Code == '/')
                             + (Annotation$Code == 'f')
                             + (Annotation$Code == '?')
                             + (Annotation$Code == 'B')
                             + (Annotation$Code == 'I')) == 1 ]  
  
  AnnCode <- Annotation$Code[  ((Annotation$Code == 'N') 
                                           + (Annotation$Code == 'A') 
                                           + (Annotation$Code == 'Q') 
                                           + (Annotation$Code == 'V') 
                                           + (Annotation$Code == 'L')
                                           + (Annotation$Code == 'R')
                                           + (Annotation$Code == 'a')
                                           + (Annotation$Code == 'J')
                                           + (Annotation$Code == 'S')
                                           + (Annotation$Code == 'r')
                                           + (Annotation$Code == 'F')
                                           + (Annotation$Code == 'e')
                                           + (Annotation$Code == 'j')
                                           + (Annotation$Code == 'n')
                                           + (Annotation$Code == 'E')
                                           + (Annotation$Code == '/')
                                           + (Annotation$Code == 'f')
                                           + (Annotation$Code == '?')
                                           + (Annotation$Code == 'B')
                                           + (Annotation$Code == 'I')) == 1 ]  
  
  
  # Calulate positive values.
  P <- length( AnnIn )
  N <- length(waveData$Date) - P  
  
  distmatrix <- abs(as.matrix(pdist(as.matrix(RWaveExtractedDataI$t) , as.matrix(waveData$Date[AnnIn]))))

  minvalue <- apply(distmatrix , 2 ,  min)
  minindex <- apply(distmatrix , 2 ,  which.min)
  MissedPeaks <- minindex[minvalue> (4*0.005)] 
  FN <- length(MissedPeaks)
  TN <- N - FN
  MissedPeakCodes <- AnnCode[MissedPeaks]
  
  minvalue <- apply(distmatrix , 1 ,  min)
  minindex <- apply(distmatrix , 1 ,  which.min)
  
  ExtraPeaks <- minindex[minvalue > (4*0.005)] 
  FP <- length(ExtraPeaks)
  TP <- P - length(ExtraPeaks)
  
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