## Clean Discrete Files
DiscData = lapply(Disc_files,function(x){readAnovFiles(pathCSV,x)})
DiscData = do.call("rbind",DiscData)
DiscData$Date = as.POSIXct(as.numeric(DiscData$Date),origin = "1970-01-01")
DiscData$YMD = format(DiscData$Date,"%Y-%m-%d")
DiscData$Time = format(DiscData$Date, "%H:%M:%S")
DiscData = subset(DiscData, !is.na(Value))

print("Discrete Combined")

#Use Patient Index to trim data
if(length(path_PatIndex)>0){
# Don't know why this is loaded in again.  
  if(filetype == 'csv'){PatIndex2017 = read.csv(file=path_PatIndex, stringsAsFactors = FALSE)}
  if(filetype == 'RData'){load(path_PatIndex)}
  
  sub_pat = subset(PatIndex2017, PseudoId %in% PatientCode)
  if(nrow(sub_pat) >0){
  FirstITUEntry = as.POSIXct(sub_pat$FirstITUEntry)
  LastITUEntry = as.POSIXct(sub_pat$LastITUEntry)
  } else{
    sub_pat = list()
    DiscData$FirstITUEntry = min(DiscData$Date)
    DiscData$LastITUEntry = max(DiscData$Date)
  }
} else {
  warning("No Patient Info provided")
  DiscData$FirstITUEntry = min(DiscData$Date)
  DiscData$LastITUEntry = max(DiscData$Date)
}
start_time = min(DiscData$Date)
stop_time = max(DiscData$Date)

pathOutFiles = paste0(path,"\\",PatientCode,"\\Disc_clean\\")
dir.create(pathOutFiles, showWarnings = FALSE)
if (choose_outputs[2] == 1){
  pathOutFilesMAT = paste0(path,"\\",PatientCode,"\\Disc_cleanMAT\\")
  dir.create(pathOutFilesMAT, showWarnings = FALSE)
}

#Remove excess data end
DiscData$InInterval = DiscData$Date > LastITUEntry[1]

DiscData = aggregate(Value ~ ., data = DiscData, mean, na.rm = TRUE)

ExtraDisc = subset(DiscData, DiscData$InInterval == TRUE)

if (dim(ExtraDisc)[1]>0){
  ExtraDisc$LastITUEntry = NULL
  ExtraDisc$FirstITUEntry = NULL
  ExtraDisc$InInterval = NULL
  if (choose_outputs[1] == 1){
    write.csv(ExtraDisc,paste0(pathOutFilesExtra,"DiscreteExtraPost.csv"), row.names=FALSE)
  }
  if (choose_outputs[2] == 1){
    writeMat(paste0(pathOutFilesExtra, "DiscreteExtraPost.mat"), Disc = ExtraDisc)
  }
  if (choose_outputs[3] == 1){
    save(file=paste0(pathOutFilesExtra, "DiscreteExtraPost.RData"), list = "ExtraDisc")
  }
}

#Remove excess data beginning
DiscData$InInterval = DiscData$Date < FirstITUEntry[1]
ExtraDisc = subset(DiscData, DiscData$InInterval == TRUE)
if (dim(ExtraDisc)[1]>0){
  ExtraDisc$LastITUEntry = NULL
  ExtraDisc$FirstITUEntry = NULL
  ExtraDisc$InInterval = NULL
  if (choose_outputs[1] == 1){
    write.csv(ExtraDisc,paste0(pathOutFilesExtra,"DiscreteExtraPre.csv"), row.names=FALSE)
  }
  if (choose_outputs[2] == 1){
    writeMat(paste0(pathOutFilesExtra, "DiscreteExtraPre.mat"), Disc = ExtraDisc)
  }
  if (choose_outputs[3] == 1){
    save(file=paste0(pathOutFilesExtra, "DiscreteExtraPre.RData"), list = "ExtraDisc")
  }
}

rm(ExtraDisc)

DiscData$InInterval = (DiscData$Date > LastITUEntry[1]) | (DiscData$Date < FirstITUEntry)

DiscDataTrim = subset(DiscData, DiscData$InInterval == FALSE)

DiscDataTrim = subset(DiscDataTrim, !is.na(Value))

DiscDataTrim$InInterval = NULL

seq_times = seq(start_time,stop_time,3*60*60)

char_name = paste0(PatientCode,"_disc_file_")

if (choose_outputs[1] == 1 | choose_outputs[2] == 1){
  for (lenX in 2:length(seq_times)){
    sub_Disc = subset(DiscDataTrim, Date<seq_times[lenX] & Date>=seq_times[lenX-1])
    if (dim(sub_Disc)[1]>0){
      subTest = aggregate(.~Date+VarName+Time+YMD, data = sub_Disc, mean)
      filename = gsub(" ","_",paste0(c(char_name,sprintf('%0.5d', lenX-1),"_",as.character(seq_times[lenX-1]),".csv"),sep = "",collapse = ""))
      filenameMAT = gsub(" ","_",paste0(c(char_name,sprintf('%0.5d', lenX-1),"_",as.character(seq_times[lenX-1]),".mat"),sep = "",collapse = ""))
      filename = gsub(":","_",filename)
      filenameMAT = gsub(":","_",filenameMAT)
      if (choose_outputs[1] == 1){
        write.csv(sub_Disc,file = paste0(pathOutFiles,filename, sep = "", collapse = ""), row.names=FALSE)
      }
      if (choose_outputs[2] == 1){
        writeMat(paste0(pathOutFilesMAT,filenameMAT, sep = "", collapse = ""), Disc = sub_Disc)
      }
    }
  }
}

if (choose_outputs[1] == 1){
  if(UseZip == 1){
  zipfiles("test.vbs", paste0(PatientCode,"_Disc_clean.zip"),pathZIPs,
           paste0(path, "\\", PatientCode, "\\Disc_clean"),
           gsub("\\","/",paste0(pathFiles,"concat_vb.txt"), fixed = TRUE))
  }
  if(Use7z == 1){
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
        "C:\\Program Files\\7-Zip", 
        sep = ";"))
    systCom = paste0('7z a ',paste0(pathZIPs, "\\",PatientCode,"_Disc_cleanCSV.7z -o"),paste0(path, "\\", PatientCode, "\\Disc_clean"), " ", paste0(path, "\\", PatientCode, "\\Disc_clean\\*"))
    system(systCom)
  }
}
if (choose_outputs[2] == 1){
  if(UseZip == 1){
  zipfiles("testMAT.vbs", paste0(PatientCode,"_Disc_cleanMAT.zip"),pathZIPs,
           gsub("\\","/",paste0(pathFiles,"concat_vb.txt"), fixed = TRUE))
  }
  if(Use7z == 1){
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
                            "C:\\Program Files\\7-Zip", 
                            sep = ";"))
    systCom = paste0('7z a ',paste0(pathZIPs, "\\",PatientCode,"_Disc_cleanMAT.7z -o"),paste0(path, "\\", PatientCode, "\\Disc_cleanMAT"), " ", paste0(path, "\\", PatientCode, "\\Disc_cleanMAT\\*"))
    system(systCom)
  }
}
if (choose_outputs[3] == 1){
  save(file = paste0(pathZIPs, "\\Discrete_",PatientCode,".RData"), list = "DiscDataTrim")
  if(Use7z == 1){
  Sys.setenv(
    PATH = paste(
      Sys.getenv("PATH"), 
      "C:\\Program Files\\7-Zip", 
      sep = ";"
    )
  )
  systCom = paste0('7z a ',paste0(pathZIPs, "\\Discrete_",PatientCode,"_RData.7z -o"),pathZIPs, " ", paste0(pathZIPs, "\\Discrete_",PatientCode,".RData"))
  system(systCom)
  }
}

rm(DiscData)
rm(DiscDataTrim)
