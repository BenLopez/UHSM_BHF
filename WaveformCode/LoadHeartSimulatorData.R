HSim_ExtractMetaData <- function(filestoprocess){
 
  MetaData = c('a','b','c')
  
  if(grepl('SR' , filestoprocess)){
    
    MetaData[1] <- 'SR'
    tmp = regexpr('SR_' , filestoprocess )[[1]]
    MetaData[2] <- substr(filestoprocess,tmp+3 ,tmp[[1]] + 3 + 2)
    tmp = regexpr(paste0(MetaData[2],'_') , filestoprocess )[[1]]
    tmp2 = regexpr('.txt' , filestoprocess )[[1]]
    MetaData[3] <- substr(filestoprocess,tmp+nchar(paste0(MetaData[2],'_') ) ,tmp2 -1)
  }
  
  if(grepl('AF' , filestoprocess)){
    
    MetaData[1] <- 'AF'
    tmp = regexpr('AF_' , filestoprocess )[[1]]
    MetaData[2] <- substr(filestoprocess,tmp+3 ,tmp[[1]] + 3 + 2)
    tmp = regexpr(paste0(MetaData[2],'_') , filestoprocess )[[1]]
    tmp2 = regexpr('.txt' , filestoprocess )[[1]]
    MetaData[3] <- substr(filestoprocess,tmp+nchar(paste0(MetaData[2],'_') ) ,tmp2 -1)
  }
  
  
  data = names( read.delim( file = paste0(filelocation ,'\\', filestoprocess),sep = ","))
  data = apply( as.matrix( data ) , 1 , function( X ){gsub( 'X','',X)} )
  data = apply( as.matrix( data ) , 1 , function( X ){gsub( '.0.','-0.', X , fixed = T)} )
  data = as.numeric( data )

  
  return(list(data = data , MetaData = MetaData ))
    
}

filelocation = choose.dir()
setwd(filelocation)
FILES <- list.files( pattern = ".txt")

AFData = list()
for(i in 1:length(FILES)){
  
  AFData[[i]] <- HSim_ExtractMetaData(FILES[[i]])
  
}

filelocation = choose.dir()
setwd(filelocation)
FILES <- list.files( pattern = ".txt")

SRData = list()
for(i in 1:length(FILES)){
  
  SRData[[i]] <- HSim_ExtractMetaData(FILES[[i]])
  
}


# Create datastructure
ListofMeshs <- unique(unlist(lapply(SRData , function(X){X$MetaData[2]})))
ListofECGs <- unique(unlist(lapply(SRData , function(X){X$MetaData[3]})))

tmp <- rep( list(1),length(ListofECGs) )
tmp <- setNames( tmp ,ListofECGs )

SimulatedPwaveDataSet <- list()
for(i in 1:length(ListofMeshs)){
  SimulatedPwaveDataSet[[i]] <- setNames(list(tmp,tmp) , c('SR','AF'))
}

SimulatedPwaveDataSet <- setNames(SimulatedPwaveDataSet , ListofMeshs)


for(i in 1:length(SRData)){
  
  SimulatedPwaveDataSet[[which(names(SimulatedPwaveDataSet) == SRData[[i]]$MetaData[2] )]]$SR[[which(names( SimulatedPwaveDataSet[[which(names(SimulatedPwaveDataSet) == SRData[[i]]$MetaData[2] )]]$SR) == SRData[[i]]$MetaData[3]) ]] <- SRData[[i]]$data
  SimulatedPwaveDataSet[[which(names(SimulatedPwaveDataSet) == AFData[[i]]$MetaData[2] )]]$AF[[which(names( SimulatedPwaveDataSet[[which(names(SimulatedPwaveDataSet) == SRData[[i]]$MetaData[2] )]]$AF) == AFData[[i]]$MetaData[3]) ]] <- AFData[[i]]$data
  
}


plot(SimulatedPwaveDataSet[[1]]$SR$I , type ='l' , ylim = c(-0.1,0.1))
for(i in 2:length(SimulatedPwaveDataSet)){
  
  lines(SimulatedPwaveDataSet[[i]]$SR$I)
}



plot(SimulatedPwaveDataSet[[1]]$AF$I , type ='l' , ylim = c(-0.1,0.1))
for(i in 2:length(SimulatedPwaveDataSet)){
  
  lines(SimulatedPwaveDataSet[[i]]$AF$I)
}