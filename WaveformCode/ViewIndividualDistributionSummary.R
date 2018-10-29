{
  if(file.exists('CheckforDefaultsScript.R')){
    source('CheckforDefaultsScript.R')
  }else{
    pathFiles <- setwd(paste0(choose.dir(caption="Select folder with source code."), "\\"))
    source("LibrariesAndSettings.R" , print.eval  = TRUE )
    DP_LoadPatientIndex()
    DP_ChooseDataReps()
    FilestoProcess <- DP_ChooseECGstoProcess() 
    HoursBeforeandAfter <- DP_SelectHoursBeforeandAfter()
  }
  listAllPatients <- DP_FilterPatients(listAllPatients , PatIndex2017 , HowtoFilterops , path , FilestoProcess)
  set.seed(1)
}

patient <- DP_choosepatient(listAllPatients = listAllPatients)
DistributionSummaries = DP_LoadDistributionSummaries(path , patient)

p1 <- ggplot( ) 
for(i in 1:(length(DistributionSummaries) -1)){
  p1 <- p1 +   geom_line(data = data.frame(x = DistributionSummaries$time , y = DP_NormaliseData(DistributionSummaries[[i]]) ) , aes(x , y)  , color = rgb(0 , 0 , i/11 , alpha = 0.5) )
}  

p1 <- p1 + ggtitle( TeX('Observation z for a patient') ) + xlab(TeX('time')) + ylab(TeX('Normalised $z_{(i)}$.'))
p1 <- p1 + geom_vline(xintercept = as.numeric(DP_StripTime(patientrecord$ConfirmedFirstNewAF)))
x11()
print(p1)