#
#
#
# Function to load and trim EEMS from Aqualog
# 29june2015
########

EEMtrim <- function(graphheadings, samplewd, loopnum) {
  
  #set wd where all of the EEM files are located. Note in .dat format
  setwd(samplewd)
  
  # Read in the EEM file
  EEMSfilename <- toString(graphheadings[loopnum,2]) # where the EEMS file are located in the graph headings file 
  
  # read in EEM file
  EEMSfile <- as.data.frame(read.delim(EEMSfilename, header= FALSE, sep = "", stringsAsFactors=FALSE))
  
  ############## trim and arrange the eem files
  # EEM file
  #set em and ex based on the columns (emission) and rows (excitation) of the raw file
  x = ncol(EEMSfile)
  ex_initial = as.numeric(EEMSfile[1, c(4:x)])
  ex = as.numeric((sort(EEMSfile[1, c(4:x)], decreasing = FALSE)))
  
  #emission = y axis
  y = nrow(EEMSfile)
  em = as.numeric(t(data.frame(EEMSfile[c(3:y), 1])))
  
  # take out the first three rows and first column of data in EEMS. These just contain text
  EEMScut <- EEMSfile[c(3:y), c(2:(x-2))]
  colnames(EEMScut) <- ex_initial
  rownames(EEMScut) <- em
  
  # first two columns are characters - return to numeric
  EEMScut[,c(1,2)] <- sapply(EEMScut[,c(1,2)],as.numeric)
  
  return(EEMScut)
}