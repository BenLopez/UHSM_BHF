

GMMmodelbyPatient <- list()

for(i in 1:dim( DataBaseMaster$NAFPatientsDatabase)[1]){
TempData1 <- DataBaseMaster$NAFPatientsDatabase[i , DataBaseMaster$NAFPatientsDatabase[i ,  ,1] !=0  ,1:11]
if(length(TempData1) < 1000){next}
GMMmodelbyPatient[[i]] <- densityMclust( DP_RemoveNaRows( TempData1 ) , G = 100 )
DP_WaitBar(i/dim(DataBaseMaster$NAFPatientsDatabase)[1])
}




