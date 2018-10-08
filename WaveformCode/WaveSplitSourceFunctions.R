unzipfile = function(outfile,filename,path){
  close( file( outfile, open="w" ) )
  print(readLines(outfile))
  ZipFile = paste0(path, filename, sep = "", collapse = "")
  ExtractTo = paste0(path, "temp_zip\\",filename, sep = "", collapse = "")
  
  ZipFile = gsub("/","\\", ZipFile, fixed = TRUE)
  ExtractTo = gsub("/","\\", ExtractTo, fixed = TRUE)
  
  write(paste0('ZipFile = "', ZipFile, '"',  sep = "", collapse =""), outfile, append = TRUE)
  write(paste0('ExtractTo = "', ExtractTo, '"',  sep = "", collapse =""), outfile, append = TRUE)
  
  write('Set fso = CreateObject("Scripting.FileSystemObject")', outfile, append = TRUE)
  write('If NOT fso.FolderExists(ExtractTo) Then', outfile, append = TRUE)
  write('fso.CreateFolder(ExtractTo)', outfile, append = TRUE)
  write('End If', outfile, append = TRUE)
  
  write('set objShell = CreateObject("Shell.Application")', outfile, append = TRUE)
  write('set FilesInZip=objShell.NameSpace(ZipFile).items', outfile, append = TRUE)
  write('objShell.NameSpace(ExtractTo).CopyHere(FilesInZip)', outfile, append = TRUE)
  write('Set fso = Nothing', outfile, append = TRUE)
  write('Set objShell = Nothing', outfile, append = TRUE)
  shell.exec(outfile)
}

unzipfileMatch = function(outfile,filename,path,strType){
  close( file( outfile, open="w" ) )
  print(readLines(outfile))
  ZipFile = paste0(path, filename, sep = "", collapse = "")
  ExtractTo = paste0(path, "temp_zip\\",filename, sep = "", collapse = "")
  
  ZipFile = gsub("/","\\", ZipFile, fixed = TRUE)
  ExtractTo = gsub("/","\\", ExtractTo, fixed = TRUE)
  
  write(paste0('ZipFile = "', ZipFile, '"',  sep = "", collapse =""), outfile, append = TRUE)
  write(paste0('ExtractTo = "', ExtractTo, '"',  sep = "", collapse =""), outfile, append = TRUE)
  
  write('Set fso = CreateObject("Scripting.FileSystemObject")', outfile, append = TRUE)
  write('If NOT fso.FolderExists(ExtractTo) Then', outfile, append = TRUE)
  write('fso.CreateFolder(ExtractTo)', outfile, append = TRUE)
  write('End If', outfile, append = TRUE)
  
  write('set objShell = CreateObject("Shell.Application")', outfile, append = TRUE)

  write('Sub ExtractWithSubfolders(PathToExtract, Pattern, OutputPath)', outfile, append = TRUE)
  write('For Each objItem in objShell.NameSpace(PathToExtract).items', outfile, append = TRUE)
  write('If objItem.Type = "File folder" Then', outfile, append = TRUE)
  write('ExtractWithSubfolders objItem.GetFolder, Pattern, OutputPath', outfile, append = TRUE)
  write('Else', outfile, append = TRUE)
  write('If Cdbl(InStr(objItem.Name,"ECG II.W"))>0 Then', outfile, append = TRUE)
  write('objShell.NameSpace(OutputPath).CopyHere(objItem)', outfile, append = TRUE)
  write('End If', outfile, append = TRUE)
  write('End If', outfile, append = TRUE)
  write('Next', outfile, append = TRUE)
  write('End Sub', outfile, append = TRUE)

  write('ExtractWithSubfolders ZipFile, "ECG II.W", ExtractTo', outfile, append = TRUE)

  write('Set fso = Nothing', outfile, append = TRUE)
  write('Set objShell = Nothing', outfile, append = TRUE)
  shell.exec(outfile)
}

zipfiles = function(outfile,
                    filename,
                    path,
                    pathSource,
                    pathExtraSource){
  close( file( outfile, open="w" ) )
  print(readLines(outfile))
  ZipFile = paste0(path, filename, sep = "", collapse = "")
  FolderToZip = pathSource
  
  write(paste0('ZipFile = "', ZipFile, '"',  sep = "", collapse =""), outfile, append = TRUE)
  write(paste0('FolderToZip = "', FolderToZip, '"',  sep = "", collapse =""), outfile, append = TRUE)
  
  extraText = readLines(paste0(pathExtraSource))
  for (testS in extraText){
    write(testS, outfile, append = TRUE)
  }
  shell.exec(outfile)
}

readAnovFiles = function(path,filename){
  con = file(paste0(path,filename,sep="",collapse=""))
  test1 = readLines(con)
  close(con)
  nchartest = nchar(test1)
  test1 = test1[nchartest == 35]
  testsplit = data.frame(VarName = as.character(trimws(gsub(",","",substr(test1,1,16)))),
                         Value = as.numeric(trimws(gsub(",","",substr(test1,17,25)))),
                         Date = trimws(gsub(",","",substr(test1,26,35))), stringsAsFactors = FALSE)
  return(testsplit)
}

cleanECG = function(filename){
  ECGfile = read.csv(paste0(pathCSV,filename),stringsAsFactors = FALSE, header = FALSE)
  names(ECGfile) = c("Date", "Value")
  ECGfile$Date = as.double(ECGfile$Date)
  ECGfile$Value = as.numeric(ECGfile$Value)
  curr_time = as.numeric(Sys.time())
  ECGfile$Date[ECGfile$Date>curr_time] = NA # Remove bad dates in future
  ECGfile$Value[ECGfile$Value == -32768] = NA 
  qcut = quantile(ECGfile$Value,c(0.005,0.995), na.rm = TRUE)
  ECGfile$Value[ECGfile$Value >= qcut[2] | ECGfile$Value <= qcut[1]] = NA
  ECGfile$Date = as.POSIXct(as.numeric(ECGfile$Date),origin = "1970-01-01")
  ECGfile = subset(ECGfile, !is.na(ECGfile$Value))
  ECGfile = subset(ECGfile, !is.na(ECGfile$Date))
  return(ECGfile)
}

cleanWave = function(filename){
  Wavefile = read.csv(paste0(pathCSV,filename),stringsAsFactors = FALSE, header = FALSE)
  names(Wavefile) = c("Date", "Value")
  Wavefile$Date = as.double(Wavefile$Date)
  Wavefile$Value = as.numeric(Wavefile$Value)
  curr_time = as.numeric(Sys.time())
  Wavefile$Date[Wavefile$Date>curr_time] = NA #Remove bad dates in future
  Wavefile$Value[Wavefile$Value == -32768] = NA
  Wavefile$Date = as.POSIXct(as.numeric(Wavefile$Date),origin = "1970-01-01")
  Wavefile = subset(Wavefile, !is.na(Wavefile$Value))
  Wavefile = subset(Wavefile, !is.na(Wavefile$Date))
  return(Wavefile)
}
