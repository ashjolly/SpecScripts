#
#
#
# Function to load and trim blank files from Aqualog
# 29june2015
########

BLANKtrim <- function(graphheadings, samplewd, loopnum) {
  
  #set wd where all of the blank files are located. Note in .dat format
  setwd(samplewd)
  
  #read in blank file
  # Blankfilename <- test2[i,2] # set blank file for the sample
  Blankfilename <- toString(graphheadings[loopnum,4]) 
  
  # Read in blank file according to name in graph headings 
  Blankfile <- read.delim(Blankfilename, header= FALSE, sep = "")
  
  #### Blank file trimming
  # take out the first two rows and first two columns of data in Blank - these contain text
  # assign column and row names
  
  #emission = y axis. get as row names
  y = nrow(Blankfile)
  em.b = as.numeric(t(data.frame(Blankfile[c(2:y), 1])))
  
  #excitation - x axis. Get as column names.
  x = ncol(Blankfile) 
  ex.b_initial = as.numeric(Blankfile[1, c(2:x)])
  ex.b = as.numeric((sort(Blankfile[1, c(2:x)], decreasing = FALSE)))
  
  # cut Blank files - first row and first column is ex and em wavelengths
  Blkcor <- Blankfile[c(2:y), c(2:x)]  
  
  #Assign col names as ex and row names as emission
  colnames(Blkcor) <- ex.b_initial
  rownames(Blkcor) <- em.b
  
  return(Blkcor)
}