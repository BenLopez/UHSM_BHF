
interactivemode <- 1
while(interactivemode == 1){
  if(nrow(PatientRecord) ==0){PatientRecord <- outputdata$Meta_Data}
  if(nrow(PatientRecord) ==0){PatientRecord <- DP_CreateDummyMetaData(PatIndex2017 , Name = subList , FirstNewAF = NA)}
{  p3 <- ggplot(ECGI[regionofinterest , ] , aes(Date , Value)) +
    geom_line(colour="blue") + ylab('Hz') + xlab('t') +
    geom_point(data = outputdata$ECGI[ ((outputdata$ECGI$t > ECGI$Date[regionofinterest[1]])*(outputdata$ECGI$t < ECGI$Date[regionofinterest[length(regionofinterest)]]) == 1) , ] , aes(t , RA) ) +
    xlim(ECGI[regionofinterest[1] , 1] , ECGI[regionofinterest[length(regionofinterest)] , 1] )  
  p5 <- ggplot(ECGII[regionofinterest2 , ] , aes(Date , Value)) +
    geom_point(data = outputdata$ECGII[ ((outputdata$ECGII$t > ECGII$Date[regionofinterest2[1]])*(outputdata$ECGII$t < ECGII$Date[regionofinterest2[length(regionofinterest2)]]) == 1) , ] , aes(t , RA) ) +
    geom_line(colour="red")+ ylab('Hz') + xlab('t') +
    xlim(ECGI[regionofinterest[1] , 1] , ECGI[regionofinterest[length(regionofinterest)] , 1] )
  p6 <- ggplot(ECGIII[regionofinterest3 , ] , aes(Date , Value)) +
    geom_point(data = outputdata$ECGIII[ ((outputdata$ECGIII$t > ECGIII$Date[regionofinterest3[1]])*(outputdata$ECGIII$t < ECGIII$Date[regionofinterest3[length(regionofinterest)]]) == 1) , ] , aes(t , RA) ) +
    geom_line(colour="black")+ ylab('Hz') + xlab('t') +
    xlim(ECGI[regionofinterest[1] , 1] , ECGI[regionofinterest[length(regionofinterest)] , 1] ) + 
    ggtitle( 'ECGIII' )
  
  plot(1)
  dev.off()
  x11(15,10)
  print(grid.arrange(  p3 + ggtitle(paste0('ECGI ' , PatientRecord$PseudoId , '  ' ,  ECGI$Date[regionofinterest[1]])) ,
                       p5 + ggtitle('ECGII') , 
                       p6 ,
                       p1 + geom_vline( xintercept = as.numeric(ECGI[regionofinterest[1] , 1]) , linetype="dashed" , color = "black" ) + 
                         geom_vline( xintercept = as.numeric(ECGI[regionofinterest[length(regionofinterest)] , 1]) , linetype="dashed" , color = "black" ) +
                         geom_vline( xintercept = as.numeric( as.POSIXct( DP_StripTime(PatientRecord$FirstNewAF)) ) , linetype="dashed" , color = "purple" ) +
                         geom_vline( xintercept = as.numeric( as.POSIXct( DP_StripTime(PatientRecord$ConfirmedFirstNewAF)) )  , color = "purple" ),
                       p2 + geom_vline( xintercept = as.numeric(ECGI[regionofinterest[1] , 1]) , linetype="dashed" , color = "black" ) + 
                         geom_vline( xintercept = as.numeric(ECGI[regionofinterest[length(regionofinterest)] , 1]) , linetype="dashed" , color = "black" ) +
                         geom_vline( xintercept = as.numeric( as.POSIXct( DP_StripTime(PatientRecord$FirstNewAF)) ) , linetype="dashed" , color = "purple" ) +
                         geom_vline( xintercept = as.numeric( as.POSIXct( DP_StripTime(PatientRecord$ConfirmedFirstNewAF)) )  , color = "purple" ) +
                         ggtitle('R-R times ' ) ,
                       
                       nrow = 5 ,
                       ncol = 1) )
  }
  Sys.sleep(0.1)
  UserResponse <- winDialog(type = c('yesnocancel') , message = 'Would you like to view another time period?')
  
  if(UserResponse == 'CANCEL')
  {
    interactivemode = 0
    break
  }
  
if(UserResponse == 'YES')
{
    
    jumpschoices <- ASWF_CreateJumpChoices()
    jump <- select.list(  jumpschoices
                          , preselect = jumpschoices[1]
                          , multiple = FALSE
                          , title = 'Select number of hours before.'
                          , graphics = TRUE )
    
    regionofinterest <- ASWF_SegmentChange(regionofinterest , jump)
    regionofinterest <- ASWF_Truncatetoregionwithdata(regionofinterest , ECGI)
    regionofinterest2 <- ASWF_Truncatetoregionwithdata(regionofinterest2 , ECGII)
    regionofinterest3 <- ASWF_Truncatetoregionwithdata(regionofinterest3 , ECGIII)
    
    
    if( round(difftime(ECGI[regionofinterest[1] , 1] , ECGII[regionofinterest2[1] , 1] , units = 'secs')) != 0  )
    {
      print('Missing data realigning region of interest')
      regionofinterest2  <- ASWF_AlignRegionofInterests(ECGI , ECGII , regionofinterest)
      regionofinterest3  <- ASWF_AlignRegionofInterests(ECGI , ECGIII , regionofinterest)
          }
    
    
    next 
}
  
if(UserResponse == 'NO'){
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to save the reduced waveforms?')
  if(UserResponse == 'YES'){
    WaveData <- ECGI
    save( WaveData , file = paste0(path , '\\' , subList , '\\Zip_out\\' ,   subList  , '_' , 'ECGI'  , '_reduced.RData' ) )
    WaveData <- ECGII
    save( WaveData , file = paste0(path , '\\' , subList , '\\Zip_out\\' ,   subList  , '_' , 'ECGII'  , '_reduced.RData' ) )
    rm(WaveData)
    WaveData <- ECGIII
    save( WaveData , file = paste0(path , '\\' , subList , '\\Zip_out\\' ,   subList  , '_' , 'ECGIII'  , '_reduced.RData' ) )
    rm(WaveData)
}
  
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to save the RpeaksFile?')
  if(UserResponse == 'YES'){
  print( 'Saving output.' )
  save( outputdata , file = paste0(path , '\\' , subList , '\\Zip_out\\' ,  subList  , '_RPeaks.RData' ) )
  print( 'Output saved.' )
  rm(outputdata)
  } 
  
  UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to view another patient?')
  if(UserResponse == 'YES'){source('ASWF_ChooseLoadandProcessPatient.R')}
  if(UserResponse == 'NO'){interactivemode <- 0 }
}
  
}