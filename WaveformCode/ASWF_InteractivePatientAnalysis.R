
interactivemode <- 1
while(interactivemode == 1){
  
  p3 <- ggplot(ECGI[regionofinterest , ] , aes(Date , Value)) +
    geom_line(colour="blue") +
    geom_point(data = RWaveExtractedDataI[ ((RWaveExtractedDataI$t > ECGI$Date[regionofinterest[1]])*(RWaveExtractedDataI$t < ECGI$Date[regionofinterest[length(regionofinterest)]]) == 1) , ] , aes(t , RA) ) +
    xlim(ECGI[regionofinterest[1] , 1] , ECGI[regionofinterest[length(regionofinterest)] , 1] )  
  p5 <- ggplot(ECGII[regionofinterest2 , ] , aes(Date , Value)) +
    geom_line(colour="blue")+
    xlim(ECGI[regionofinterest[1] , 1] , ECGI[regionofinterest[length(regionofinterest)] , 1] )  
  
  dev.off()
  x11(15,10)
  print(grid.arrange(  p3 + ggtitle(paste0('ECGI ' , DataSet$MetaData$PseudoId , '  ' ,  ECGI$Date[regionofinterest[1]])) ,
                       p5 + ggtitle('ECGII') , 
                       p1 + geom_vline( xintercept = as.numeric(ECGI[regionofinterest[1] , 1]) , linetype="dashed" , color = "black" ) + 
                         geom_vline( xintercept = as.numeric(ECGI[regionofinterest[length(regionofinterest)] , 1]) , linetype="dashed" , color = "black" ) +
                         geom_vline( xintercept = as.numeric( as.POSIXct( DataSet$MetaData$FirstNewAF) ) , linetype="dashed" , color = "purple" ) ,
                       p2 + geom_vline( xintercept = as.numeric(ECGI[regionofinterest[1] , 1]) , linetype="dashed" , color = "black" ) + 
                         geom_vline( xintercept = as.numeric(ECGI[regionofinterest[length(regionofinterest)] , 1]) , linetype="dashed" , color = "black" ) +
                         geom_vline( xintercept = as.numeric( as.POSIXct( DataSet$MetaData$FirstNewAF) ) , linetype="dashed" , color = "purple" ) +
                         ggtitle('R-R times ' ) ,
                       p4,
                       nrow = 5 ,
                       ncol = 1) )
  Sys.sleep(0.1)
  UserResponse <- winDialog(type = c('yesnocancel') , message = 'Would you like to view another time period?')
  
  if(UserResponse == 'CANCEL')
  {
    interactivemode = 0
    break
  }
  
if(UserResponse == 'YES')
{
    
    jumpschoices <- c('next Segment' , 'next 10' , 'next 100' , 'next 1000' , 'previous Segment' , 'previous 10' , 'previous 100' , 'previous 1000'  )
    jump <- select.list(  jumpschoices
                          , preselect = jumpschoices[1]
                          , multiple = FALSE
                          , title = 'Select number of hours before.'
                          , graphics = TRUE )
    
    
    regionofinterest <- ASWF_SegmentChange(regionofinterest , jump)
    regionofinterest2 <- ASWF_SegmentChange(regionofinterest2 , jump)
    regionofinterest <- ASWF_Truncatetoregionwithdata(regionofinterest , ECGI)
    regionofinterest2 <- ASWF_Truncatetoregionwithdata(regionofinterest2 , ECGII)
    
    if( round(difftime(ECGI[regionofinterest[1] , 1] , ECGII[regionofinterest2[1] , 1] , units = 'secs')) != 0  )
    {
      print('Missing data realigning region of interest')
      regionofinterest2  <- ASWF_AlignRegionofInterests(ECGI , ECGII , regionofinterest)
    }
    
    
    next 
}
  
if(UserResponse == 'NO')
{
    UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to view another patient?')
    if(UserResponse == 'YES'){source('ASWF_ChooseLoadandProcessPatient.R')}
    if(UserResponse == 'NO'){ interactivemode <- 0 }
}
  
}