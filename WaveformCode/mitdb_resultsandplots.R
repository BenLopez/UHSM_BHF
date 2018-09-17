{if(file.exists('CheckforDefaultsScript.R')){
  source('CheckforDefaultsScript.R')
}else{
  pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
  source("LibrariesAndSettings.R" , print.eval  = TRUE )
  DP_LoadPatientIndex()
  DP_ChooseDataReps()
}
}

Files <- choose.files()
patientlist <- as.matrix(lapply(Files , function(X){ substr(X , 30 , nchar(X) - 27)  }))
listofannotations <- mitdb_createlistofannotations()


datatable <- matrix(0, 6 , length(listofannotations) + 1)
rownames(datatable) <- c('Number Positives' , 'Number negatives' , 'True positive' , 'True negative' , 'False negative' , 'False postive')
colnames(datatable) <- c(listofannotations , 'Total')

Sensitivity <- matrix(0 , length(Files) , 1)
Specificity <- matrix(0 , length(Files) , 1)
Accuracy <- matrix(0 , length(Files) , 1)
PPV <- matrix(0 , length(Files) , 1)
NPV <- matrix(0 , length(Files) , 1)


for (i in 1:length(Files)){
load(Files[[i]])
  Sensitivity[i] <- Results$Sensitivity
  Specificity[i] <- Results$Specificity
  Accuracy[i] <- Results$Accuracy
  PPV[i] <- Results$PPV
  NPV[i] <- Results$NPV
  
datatable[ 1 , 20] <- datatable[ 1 , 20] + Results$SufficientStatistics$P
datatable[ 2 , 20] <- datatable[ 2 , 20] + Results$SufficientStatistics$N
datatable[ 3 , 20] <- datatable[ 3 , 20] + Results$SufficientStatistics$TP
datatable[ 4 , 20] <- datatable[ 4 , 20] + Results$SufficientStatistics$TN
datatable[ 5 , 20] <- datatable[ 5 , 20] + Results$SufficientStatistics$FP
datatable[ 6 , 20] <- datatable[ 6 , 20] + Results$SufficientStatistics$FN


for(j in 1:length(listofannotations) )
{
  
if(sum(Results$BeatTypes$values == listofannotations[j]) >0 ){
datatable[ 1 , j] <- datatable[ 1 , j]  +  Results$BeatTypes$n[as.vector(Results$BeatTypes$values) == listofannotations[j]]
datatable[ 3 , j] <- datatable[ 3 , j]  +  Results$BeatTypes$n[as.vector(Results$BeatTypes$values) == listofannotations[j]] - sum(as.vector(Results$Missedpeakcodes) == listofannotations[j])
datatable[ 5 , j] <- datatable[ 5 , j]  +  sum(Results$Missedpeakcodes == listofannotations[j])
}
}
}
listofannotations <- listofannotations[apply(datatable , 2 , sum) > 0]
datatable <- datatable[ , apply(datatable , 2 , sum) > 0]
datatable[datatable == 0] = NA


p1 <- ggplot( data.frame(x =  c(1:length(Sensitivity)) ,Sensitivity , Specificity , Accuracy , PPV , NPV ) , aes(x = x)) +
    geom_point(aes(y = Sensitivity, colour ='black')) +
    geom_point(aes(y = Specificity, colour ='red')) +
    geom_point(aes(y = Accuracy, colour ='blue')) +
    geom_point(aes(y = PPV , colour ='green')) +
    geom_point(aes(y = NPV, colour ='yellow')) +
  scale_colour_manual(labels = c('Sensitivity' , 'Specificity' , 'Accuracy' , 'PPV' , 'NPV'), values = c("black","red", "green", "blue" , 'yellow')) +
  ggtitle('Analysis by Patient') +
  xlab('Patient Index') +
  ylab('Score')  

listofannotations[length(listofannotations)] = 'Total'
p2 <- ggplot( data.frame( x = as.factor(c(1:dim(datatable)[2])),y = datatable[3 , ]/datatable[1 , ]  ) , aes(x ,  y)) +
      geom_point( ) +
      ggtitle('Analysis by Beat Type ') +
      ylab('Sensitivity') +
      scale_x_discrete( labels=listofannotations , name = 'Beat Type'  )

x11(14 ,7)
grid.arrange(p1 , p2 , nrow = 1 , ncol = 2)
