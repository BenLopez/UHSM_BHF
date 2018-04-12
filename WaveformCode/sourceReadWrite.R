readWriteWave = function(Wave_files, cleanWave, sub_pat, choose_outputs, wavename,path,pathOutFilesExtra,use7z,useZip){

  
  WaveData = lapply(Wave_files, cleanWave)
  print(paste0("Finished reading csv files - ",wavename))
  
  WaveDataT = do.call("rbind",WaveData)
  print("Rbind complete")
  print(dim(WaveDataT))
  
  if (nrow(WaveDataT)){
    WaveData = WaveDataT
    rm(WaveDataT)
    WaveData$Value = as.numeric(WaveData$Value)
    if (grepl("ECG", wavename)){
      qcut = c(-200,200) # unusual calibration threshold - ask Sam
      WaveData$Value[WaveData$Value >= qcut[2] | WaveData$Value <= qcut[1]] = NA
      WaveData = subset(WaveData, !is.na(WaveData$Value))
    }
    print(dim(WaveData))
    if(length(sub_pat)>0){
      FirstITUEntry = as.POSIXct(sub_pat$FirstITUEntry)
      LastITUEntry = as.POSIXct(sub_pat$LastITUEntry)
    } else {
      warning("No Patient Info provided")
      FirstITUEntry = min(WaveData$Date)
      LastITUEntry = max(WaveData$Date)
    }
    
    #Post Dates
    WaveData$InInterval = WaveData$Date > LastITUEntry
    ExtraWave = subset(WaveData, WaveData$InInterval == TRUE)
    print(dim(ExtraWave))
    
    if (dim(ExtraWave)[1]>0){
      ExtraWave$InInterval = NULL
      if (choose_outputs[1] == 1){
        write.csv(ExtraWave,paste0(pathOutFilesExtra,wavename,"ExtraPost.csv"), row.names=FALSE)
      }
      if (choose_outputs[2] == 1){
        writeMat(paste0(pathOutFilesExtra,wavename,"ExtraPost.mat"), Wave = ExtraWave)
      }
      if (choose_outputs[3] == 1){
        save(file=paste0(pathOutFilesExtra,wavename,"ExtraPost.RData"), list = "ExtraWave")
      }}
    
    
    #Pre Dates
    WaveData$InInterval = WaveData$Date < FirstITUEntry
    ExtraWave = subset(WaveData, WaveData$InInterval == TRUE)
    print(dim(ExtraWave))
    
    if (dim(ExtraWave)[1]>0){
      ExtraWave$InInterval = NULL
      if (choose_outputs[1] == 1){
        write.csv(ExtraWave,paste0(pathOutFilesExtra,wavename,"ExtraPre.csv"), row.names=FALSE)
      }
      if (choose_outputs[2] == 1){
        writeMat(paste0(pathOutFilesExtra,wavename,"ExtraPre.mat"), Wave = ExtraWave)
      }
      if (choose_outputs[3] == 1){
        save(file=paste0(pathOutFilesExtra,wavename,"ExtraPre.RData"), list = "ExtraWave")
      }
    }
    
    WaveData$InInterval = (WaveData$Date > LastITUEntry) | (WaveData$Date < FirstITUEntry)
    
    WaveData = subset(WaveData, WaveData$InInterval == FALSE)
    
    WaveData = subset(WaveData, !is.na(Value))
    
    WaveData$InInterval = NULL
    
    start_time = min(WaveData$Date, na.rm = TRUE)
    stop_time = max(WaveData$Date, na.rm = TRUE)
    
    seq_times = unique(c(seq(start_time,stop_time,3*60*60),stop_time))
    
    pathOutFiles = paste0(path,wavename,"_clean/")
    dir.create(pathOutFiles, showWarnings = FALSE)
    pathOutFilesMAT = paste0(path,wavename,"_cleanMAT/")
    dir.create(pathOutFilesMAT, showWarnings = FALSE)
    
    char_name = paste0(PatientCode,"_",wavename,"_file_")
    if (choose_outputs[1] == 1 | choose_outputs[2] == 1){
      for (lenX in 2:length(seq_times)){
        sub_Wave = subset(WaveData, Date<seq_times[lenX] & Date>=seq_times[lenX-1])
        if (dim(sub_Wave)[1]>0){
          filename = gsub(" ","_",paste0(c(char_name,sprintf('%0.5d', lenX-1),"_",as.character(seq_times[lenX-1]),".csv"),sep = "",collapse = ""))
          filenameMAT = gsub(" ","_",paste0(c(char_name,sprintf('%0.5d', lenX-1),"_",as.character(seq_times[lenX-1]),".mat"),sep = "",collapse = ""))
          filename = gsub(":","_",filename)
          filenameMAT = gsub(":","_",filenameMAT)
          if (choose_outputs[1] == 1){
            write.csv(sub_Wave,file = paste0(pathOutFiles,filename, sep = "", collapse = ""), row.names=FALSE)
          }
          if (choose_outputs[2] == 1){
            writeMat(paste0(pathOutFilesMAT,filenameMAT, sep = "", collapse = ""), Wave = sub_Wave)
          }
        }
      }
    }
    
    if (choose_outputs[3] == 1){
      save(file = paste0(pathZIPs, wavename, "_",PatientCode,".RData"), list = "WaveData", compress = "bzip2", compression_level = )
      if(Use7z == 1){
        Sys.setenv(
          PATH = paste(
            Sys.getenv("PATH"), 
            "C:\\Program Files\\7-Zip", 
            sep = ";"
          )
        )
        systCom = paste0('7z a ',paste0(pathZIPs, wavename,"_",PatientCode,"_RData.7z -o"),pathZIPs, " ", paste0(pathZIPs, wavename, "_",PatientCode,".RData"))
        system(systCom)
      }
    }
    
    if (choose_outputs[1] == 1){
      if(UseZip == 1){
        zipfiles("test.vbs", paste0(PatientCode,"_",wavename,"_clean.zip"),pathZIPs, paste0(path, wavename,"_clean"))
      }
      if(Use7z == 1){
        Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
                                "C:\\Program Files\\7-Zip", 
                                sep = ";"))
        systCom = paste0('7z a ',paste0(pathZIPs, "\\",PatientCode,"_",wavename,"_clean.7z -o"),paste0(path, wavename, "_clean"), " ", paste0(path, wavename,"_clean\\*"))
        system(systCom)
      }
    }
    if (choose_outputs[2] == 1){
      if(UseZip == 1){
        zipfiles("testMAT.vbs", paste0(PatientCode,wavename,"_cleanMAT.zip"),pathZIPs, paste0(path,wavename,"_cleanMAT"))
      }
      if(Use7z == 1){
        Sys.setenv(PATH = paste(Sys.getenv("PATH"), 
                                "C:\\Program Files\\7-Zip", 
                                sep = ";"))
        systCom = paste0('7z a ',paste0(pathZIPs, "\\",PatientCode,"_",wavename,"_cleanMAT.7z -o"),paste0(path, wavename, "_cleanMAT"), " ", paste0(path, wavename,"_cleanMAT\\*"))
        system(systCom)
      }
    }
  } else {
    print("No usable records - empty")
  }
}
