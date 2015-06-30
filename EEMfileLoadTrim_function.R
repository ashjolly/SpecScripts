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
  EEMSfile <- read.delim(EEMSfilename, header= FALSE, sep = "")
  #samplename <- strapplyc(EEMSfilename, "001(.*)PEM", simplify = TRUE)
  samplename <- toString(graphheadings[loopnum,1]) #column with the sample IDs

  ############## trim and arrange the eem files
  # EEM file
  #set em and ex based on the columns (emission) and rows (excitation) of the raw file
  x = ncol(EEMSfile)
  ex_initial = as.numeric(EEMSfile[1, c(4:x)])
  ex = as.numeric((sort(EEMSfile[1, c(4:x)], decreasing = FALSE)))
  
  #ex_all is a variable that holds all of the ex wavelengths. It's a way of seeing if ex wavelengths are different between samples
  ex_all[loopnum,] = ex

  #emission = y axis
  y = nrow(EEMSfile)
  em = as.numeric(t(data.frame(EEMSfile[c(3:y), 1])))
  
  #em_all is a variable that holds all of the em wavelengths. It's a way of seeing if em wavelengths are different between samples
  em_all[loopnum,] = em

  # take out the first three rows and first column of data in EEMS. These just contain text
  EEMScut <- EEMSfile[c(3:y), c(2:x)]
  colnames(EEMScut) <- ex_initial
  rownames(EEMScut) <- em

  return(EEMScut)
}