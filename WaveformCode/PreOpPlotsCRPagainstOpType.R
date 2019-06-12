Proc <- PatIndex2017$ProcDetails
Proc2 <- as.matrix(Proc)
LevelNames <- c('Aortic' , 'Complex' , 'CABG' , 'Transplant' , 'Valve' , 'MCS' , 'SP' , 'PP' , 'Other')

Proc2[ Proc %in% c('Aortic','Aortic and CABG','Aortic and Valve')] <- LevelNames[1]
Proc2[ Proc %in% c(' Aortic and CABG and Other', 'Aortic and CABG and Other', ' Aortic and Valve and CABG and Other','Aortic and Valve and CABG and Other','Aortic and Valve and Other','CABG and Valve and Other' , "Aortic and Valve and CABG" )] <-  LevelNames[2]
Proc2[  Proc %in% c('CABG','CABG and Other','CABG and Valve') ] <-  LevelNames[3]
Proc2[ Proc %in% c('Cardiac transplant','Cardiac Transplant')] <-  LevelNames[4]
Proc2[ Proc %in% c("Valve and LVAD insertion",'valve and LVAD insertion', 'Valve and Other')] <-  LevelNames[5]
Proc2[ Proc %in% c('ecmo decannulation','ECMO decannulation','explantation of BiVAD',"BIVAD Insertion",'LVAD Insertion','LVAD Removal','LVAD Replacement','Perc RVAD Insertion',"LV aneurysmectomy;",'Peripheral VA ECMO','Primary VAD','Primary VAD;Other procedure not listed above;','PrimaryVAD;ASD Closure','LVAD Replacement','RVAD insertion','VA ECMO','VA ECMO (off on unit)','VA ECMO insertion','VA ECMO Insertion','VAD Explantation')] <-  LevelNames[6]
Proc2[ Proc %in% c('ASD closure','Atrial myxoma','Atrial myxoma and Other','LV aneurysmectomy')] <-  LevelNames[7]
Proc2[ Proc %in% c('drain removal','Major Sternal Debridement/VAC','exploration secondary to flow issues','Oxygenator removal','pericardial collection drainage','Pericardial Window','Pericardiectomy','removal of chest drain','Removal Of Packs/Secondary Chest Closure','Removal Sternal Wire/Debdridement','Sternal Flap','sternal rewire','Sternal wound flap','Tamponade/Re-Operation For Bleeding',"Surgical Tracheostomy")] <-  LevelNames[8]
Proc2[ Proc %in% c('Stenting carotid pseudo aneurysm','Leg fasciectomy','embolisation of gi bleed','Laparotomy' , 'Epicardial pacemaker;Other procedure not listed above;','BK amuptation')] <-LevelNames[9]

rownames(Proc2) <- PatIndex2017$PseudoId

Proc2 <- as.factor(Proc2)

# some summaries to check the data.
levels(Proc2)
summary(Proc2)

CRP <-  matrix(0 , dim(as.matrix(Proc2) )[1] , 1)
for(i in 1:dim(as.matrix(Proc2) )[1]){
  if(sum(names(BioChemIndex2017) == PatIndex2017$NewPseudoId[i]) ==0){
    CPB[i] <- NA
    next}
  CRP[i] <- BioChemIndex2017[[which(names(BioChemIndex2017) == PatIndex2017$NewPseudoId[i])]]$TimeSeriesData[1 , 5]
}

x11(width = 100,height = 40)
par(mfrow = c(1 ,1))
plot(Proc2[!is.na(PatIndex2017$FirstNewAF)] , CPB[!is.na(PatIndex2017$FirstNewAF)] , xlab = 'Operation Type' , ylab = 'CRP'  , col = rgb(1,0,0,alpha = 0.1))
title('CRP by Operation Type AF (Red) NAF (Blue)')


# plottwo
plot(Proc2[is.na(PatIndex2017$FirstNewAF)]  , CPB[is.na(PatIndex2017$FirstNewAF)]  , xlab = 'Operation Type' , ylab = 'CRP'  , add = T , col = rgb(0,0,1,alpha = 0.1))

# plotone
x11(width = 100,height = 40)
AFLogical  = !is.na(PatIndex2017$FirstNewAF)
AFLogical[AFLogical == T] = 'AF'
AFLogical[AFLogical == F] = 'NAF'
plot(Proc2 , as.factor(AFLogical) , xlab = 'Operation Type'  , ylab = 'Post Operative Atrial Fibrillation',xaxt="n" , cex=0.05)
title('Operation Type and AF Prevelance')
