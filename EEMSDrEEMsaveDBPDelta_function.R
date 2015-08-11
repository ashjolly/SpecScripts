#
# Function to take corrected EEMS and compile the EEMS in a format for modelling through Dr EEMS or DOM Fluor package
# 6 Aug 2015
# Ashlee Jollymore's PhD Project
##################

DrEEM = function(filelist, project, exmin, filedirectory, ex) {

  # set working directory with sample files
  setwd(filedirectory) 
  
  # create sample ID list to populate
  sampleID = data.frame((0))
  
  n = length(filelist)
  for (i in 1:n){
    temp.EEMS <- read.delim(filelist[i], header= TRUE, sep = ",")
  
    #trim so that exitation and emission goes from the same
    ex.temp <- colnames(temp.EEMS)
  
    if(ex.temp[1] != exmin) {
      # if first value in ex.temp is not 240, trim 
      ex.length <- length(ex.temp)
      # find column where the exitation wavelength is 240 to cut from
      x240 = as.numeric(match(exmin,names(temp.EEMS)))
      temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
      } 
    
    #create em variables from trimmed EEMS
    em = row.names(temp.EEMS)
    
    # create a new dataset where the post-cut EEMS are compiled together by rows
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- temp.EEMS
      }
  
    # if the merged dataset does exist, append to it
    if (exists("dataset")){
     temp_dataset <-temp.EEMS
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
    
    # create list with sample IDs in it
    samplename <- strapplyc(filelist_EEMScor[i], paste("(.*)_", "DBPdelta", sep = ""), simplify = TRUE)
    sampleID[i] <- samplename
  }
  
  #seems to have doubled first dataset, remove?
  x <- length(em)
  y <- dim(dataset)[1]
  dataset.2 <- dataset[c((x+1):y),]
  y <- dim(dataset.2)[1]
  remove(x)
  remove(y)
  
  # Save files in correct locations
  corrpath <- file.path("/Users/ashlee/Documents/MATLAB/ExEmfiles", paste(project,"ex",".csv", sep = ""))

  write.table(dataset.2, file = file.path("/Users/ashlee/Documents/MATLAB/toolbox/DOMFluor", paste(project, "/fl.csv", sep = "")),
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

  #Ex file
  write.table(ex, file = file.path("/Users/ashlee/Documents/MATLAB/toolbox/DOMFluor", paste(project,"/Ex.csv", sep = "")),
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

  #Em
  write.table(em, file = file.path("/Users/ashlee/Documents/MATLAB/toolbox/DOMFluor", paste(project,"/Em.csv", sep = "")), 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

  #File containing sample names
  write.table(sampleID, file = file.path("/Users/ashlee/Documents/MATLAB/toolbox/DOMFluor", paste(project,"/01key.csv",sep = "")), 
            row.names = FALSE, col.names = FALSE, sep = ",") #saved in matlab folder

  return(dataset.2)
}