
if( BCOptions[[2]] == 'AFClassifier' & BCOptions$DataType == 'DistributionSummaries' ){
{
x11(20,14)
  par(mfrow = c(2 , 5))
  
  for(variabletoview in c(1:10)){
    tmp <- hist(DataBase[[1]][ , variabletoview], col=rgb(0,0,1,alpha = 0.5) ,
                main= paste0(AFD_CreateDistributionSummaryNames()[variabletoview] , ' Histogram') , xlab = AFD_CreateDistributionSummaryNames()[variabletoview] , freq = FALSE)
    hist(DataBase[[2]][ , variabletoview], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(min(DataBase[[2]][!is.na(DataBase[[2]][ , variabletoview]) , variabletoview] ) , tmp$breaks, max(DataBase[[2]][!is.na(DataBase[[2]][ , variabletoview]) , variabletoview]) ))
  }
  title(' Local Histograms' , outer=TRUE)
  
}

{
x11(30,20)
SAMPLENUM <- 2000
variable = c(1:11)
tmp <- rbind(DataBase[[1]][sample(1:size(DataBase[[1]])[1] , SAMPLENUM ),variable] , DataBase[[2]][sample(1:size(DataBase[[2]])[1] , SAMPLENUM ),variable])
al = 0.05
colvector <- c(rep(rgb(0 ,0 , 1, alpha = al)   , SAMPLENUM) ,  rep(rgb(1 , 0 , 0, alpha = al)   , SAMPLENUM))
pairs(tmp , pch = 16 , col = colvector , labels = AFD_CreateDistributionSummaryNames()[variable])
}
}

if(BCOptions[[2]] == 'AFClassifier' &  BCOptions$GlobalUpdate == 'Yes' & BCOptions$DataType == 'DistributionSummaries'){
{
x11(20,14)
par(mfrow = c(2 , 5))
for(variabletoview in c(1:10)){
  
  tmp <- hist(DataBase[[3]][ , variabletoview], col=rgb(0,0,1,alpha = 0.5) ,
              main= paste0(AFD_CreateDistributionSummaryNames()[variabletoview] , ' Histogram') , xlab = AFD_CreateDistributionSummaryNames()[variabletoview] , freq = FALSE)
  hist(DataBase[[4]][ , variabletoview], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(min(DataBase[[4]][!is.na(DataBase[[4]][ , variabletoview]) , variabletoview] ) , tmp$breaks, max(DataBase[[4]][!is.na(DataBase[[4]][ , variabletoview]) , variabletoview]) ))
}
title(' Global Histograms' , outer=TRUE)
x11(30,20)
variable = c(1:10)
tmp <- rbind(DataBase[[3]][,variable] , DataBase[[4]][sample(1:dim(DataBase[[4]])[1] , dim(DataBase[[3]])[1]),variable])
al = 0.1
colvector <- c(rep(rgb(0 ,0 , 1, alpha = al)   , length(DataBase[[3]][,1])) ,  rep(rgb(1 , 0 , 0, alpha = al)   , length(DataBase[[3]][,1])))
pairs(tmp , pch = 16 , col = colvector , labels = AFD_CreateDistributionSummaryNames()[variable])
title(' Global Pairs' , outer=TRUE)
}
}

if(BCOptions[[2]] == 'AFClassifier' & BCOptions$DataType == 'CDFs'){
  {
    x11(20,14)
    par(mfrow = c(5 , 6))
    for(variabletoview in c(1:30)){
      
      tmp <- hist(DataBase[[3]][ , variabletoview], col=rgb(0,0,1,alpha = 0.5) ,
                  main= paste0(AFD_CreateDefaultSettings()$BinlimsScore[variabletoview] , ' Histogram') , xlab = AFD_CreateDistributionSummaryNames()[variabletoview] , freq = FALSE)
      hist(DataBase[[4]][ , variabletoview], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(min(DataBase[[4]][!is.na(DataBase[[4]][ , variabletoview]) , variabletoview] ) , tmp$breaks, max(DataBase[[4]][!is.na(DataBase[[4]][ , variabletoview]) , variabletoview]) ))
    }
    title(' Global Histograms' , outer=TRUE)
    x11(30,20)
    variable = c(6:16)
    tmp <- rbind(DataBase[[3]][,variable] , DataBase[[4]][sample(1:dim(DataBase[[4]])[1] , dim(DataBase[[3]])[1]),variable])
    al = 0.1
    colvector <- c(rep(rgb(0 ,0 , 1, alpha = al)   , length(DataBase[[3]][,1])) ,  rep(rgb(1 , 0 , 0, alpha = al)   , length(DataBase[[3]][,1])))
    pairs(tmp , pch = 16 , col = colvector , labels = AFD_CreateDefaultSettings()$BinlimsScore[variable])
    title(' Global Pairs' , outer=TRUE)
  }
}

if( BCOptions[[2]] == 'AFClassifier' & BCOptions$DataType == 'CDFs' ){
    x11(20,14)
    par(mfrow = c(5 , 6))
    
    for(variabletoview in c(1:30)){
      tmp <- hist(DataBase[[1]][ , variabletoview], col=rgb(0,0,1,alpha = 0.5) ,
                  main= paste0(AFD_CreateDefaultSettings()$BinlimsScore[variabletoview] , ' Histogram') ,
                  xlab = AFD_CreateDefaultSettings()$BinlimsScore[variabletoview] , freq = FALSE)
      hist(DataBase[[2]][ , variabletoview], col=rgb(1,0,0,alpha =0.5), add=T , freq = FALSE ,  breaks = c(min(DataBase[[2]][!is.na(DataBase[[2]][ , variabletoview]) , variabletoview] ) , tmp$breaks, max(DataBase[[2]][!is.na(DataBase[[2]][ , variabletoview]) , variabletoview]) ))
    }
    title(' Local Histograms' , outer=TRUE)
  
    {
      x11(30,20)
      SAMPLENUM <- 10000
      variable = c(6:16)
      tmp <- rbind(DataBase[[1]][sample(1:size(DataBase[[1]])[1] , SAMPLENUM ),variable] , DataBase[[2]][sample(1:size(DataBase[[2]])[1] , SAMPLENUM ),variable])
      al = 0.01
      colvector <- c(rep(rgb(0 ,0 , 1, alpha = al)   , SAMPLENUM) ,  rep(rgb(1 , 0 , 0, alpha = al)   , SAMPLENUM))
      pairs(tmp , pch = 16 , col = colvector , labels = AFD_CreateDefaultSettings()$BinlimsScore[variable])
    }
    title(' Local Pairs' , outer=TRUE)      
    
}
