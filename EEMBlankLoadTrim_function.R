#
#
#
# Function to load and trim blank files from Aqualog
# 29june2015
########

BLANKtrim <- function(graphheadings, samplewd) {
  
  #set wd where all of the blank files are located. Note in .dat format
  setwd(samplewd)
  
  #read in blank file
  # Blankfilename <- test2[i,2] # set blank file for the sample
  Blankfilename <- toString(graphheadings[i,4]) 
  
  # Read in blank file according to name in graph headings 
  Blankfile <- read.delim(Blankfilename, header= FALSE, sep = "")
  
  #### Blank file trimming
  # take out the first two rows and first two columns of data in Blank - these contain text
  # assign column and row names
  
  #emission = y axis. get as row names
  y = nrow(Blankfile)
  em = as.numeric(t(data.frame(Blankfile[c(3:y), 1])))
  em_all[i,] = em
  
  #excitation - x axis. Get as column names.
  x = ncol(Blankfile) 
  ex_initial = as.numeric(Blankfile[1, c(4:x)])
  ex = as.numeric((sort(Blankfile[1, c(4:x)], decreasing = FALSE)))
  ex_all[i,] = ex
  
  # cut Blank files
  Blkcor <- Blankfile[c(2:y), c(2:x)]  
  
  #Assign col names as ex and row names as emission
  colnames(Blkcor) <- ex_initial
  rownames(Blkcor) <- em
  
  return(Blkcor)
}