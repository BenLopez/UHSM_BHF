
interactivemode <- 1
while(interactivemode == 1){
  
  p3 <- ggplot(WaveData[regionofinterest , ] , aes(Date , Value)) +
    geom_line(colour="blue") +
    geom_point(data = RWaveExtractedData[ ((RWaveExtractedData$t > WaveData$Date[regionofinterest[1]])*(RWaveExtractedData$t < WaveData$Date[regionofinterest[length(regionofinterest)]]) == 1) , ] , aes(t , RA) )
  dev.off(4)
  x11(12,7)
  print(grid.arrange(  p3 ,
                       p1 + geom_vline( xintercept = as.numeric(WaveData[regionofinterest[1] , 1]) , linetype="dashed" , color = "black" ) + 
                         geom_vline( xintercept = as.numeric(WaveData[regionofinterest[length(regionofinterest)] , 1]) , linetype="dashed" , color = "black" ) ,
                       p2 + geom_vline( xintercept = as.numeric(WaveData[regionofinterest[1] , 1]) , linetype="dashed" , color = "black" ) + 
                         geom_vline( xintercept = as.numeric(WaveData[regionofinterest[length(regionofinterest)] , 1]) , linetype="dashed" , color = "black" )+
                         ggtitle(paste0('R-R times ' , WaveData$Date[regionofinterest[1]] )) ,
                       p4,
                       nrow = 4 ,
                       ncol = 1))
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
    
    if( jump == jumpschoices[1] ){ regionofinterest <- regionofinterest + length(regionofinterest) }
    if( jump == jumpschoices[2] ){ regionofinterest <- regionofinterest + 10*length(regionofinterest) }
    if( jump == jumpschoices[3] ){ regionofinterest <- regionofinterest + 100*length(regionofinterest) }
    if( jump == jumpschoices[4] ){ regionofinterest <- regionofinterest + 1000*length(regionofinterest) }
    if( jump == jumpschoices[5] ){ regionofinterest <- regionofinterest - length(regionofinterest) }
    if( jump == jumpschoices[6] ){ regionofinterest <- regionofinterest - 10*length(regionofinterest) }
    if( jump == jumpschoices[7] ){ regionofinterest <- regionofinterest - 100*length(regionofinterest) }
    if( jump == jumpschoices[8] ){ regionofinterest <- regionofinterest - 1000*length(regionofinterest) }
    if( regionofinterest[1] <= 0 ){ regionofinterest <-  c(1:length(regionofinterest)) }
    if( regionofinterest[length(regionofinterest)]  >=  length(WaveData[ , 1]) ){regionofinterest  <- c(length(WaveData[ , 1]) - length(regionofinterest)):length(WaveData[ , 1])}
    next 
}
  
if(UserResponse == 'NO')
{
    UserResponse <- winDialog(type = c('yesno') , message = 'Would you like to view another patient?')
    if(UserResponse == 'YES'){source('ASWF_ChooseLoadandProcessPatient.R')}
    if(UserResponse == 'NO'){ interactivemode <- 0 }
}
  
}