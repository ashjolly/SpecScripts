###
# function for trimming files for CM modelling
# Script is meant to take corrected EEM files and ensure that they have the same EX and EM ranges 
# so that file can be passed into CM model in matlab

# 8July2015
############

CMtrim <- function(directory, projectname, minex){
  
  # get directory with all of the corrected files
  setwd(directory) 
  filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")

  ######## Prepping files for Cory McKnight modelling in Matlab
  ########
  # Ensure that the files have the same ex and em lengths
  #myfiles = lapply(filelist_EEMScor, read.delim)
  #myfiles = do.call("rbind", lapply(filelist_EEMScor, 
  #    function(x) read.delim(x,stringsAsFactors = FALSE)))
  
  # CM - take out row and column names in first column and row and save in CM folder
  n = length(filelist_EEMScor)
  
  # get the ex and em for each of the files
  
  for (i in 1:n){
    temp.EEMS <- read.delim(filelist_EEMScor[i], header= TRUE, sep = ",")
  
    #trim so that exitation and emission goes from the same
    ex.temp <- colnames(temp.EEMS)
  
      if(ex.temp[1] != minex) {
      # if first value in ex.temp is not 240, trim 
      ex.length <- length(ex.temp)
      # find column where the exitation wavelength is 240 to cut from
      xmin = as.numeric(match(minex,names(temp.EEMS)))
      temp.EEMS <- temp.EEMS[,c(xmin:ex.length)]
  } 
  
  # cut out any columns containing Nas- this is 798 and 800 nm. Must cut last four rows of data from 20april2015
  #temp.EEMS.1 <- na.omit(temp.EEMS)
  #g <- length(temp.EEMS)
  #temp.EEMS.1 <- temp.EEMS[,c(1:(g-4))] #cut out the last four colomns manually
  
  #resave without the row and column names
  samplename <- strapplyc(filelist_EEMScor[i], "(.*)Prechlor", simplify = TRUE)
  corrpath <- file.path(directoryCM, paste(samplename,"CorrCM_",i,".csv", sep = ""))
  write.table(temp.EEMS.1, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")
  
}
}
