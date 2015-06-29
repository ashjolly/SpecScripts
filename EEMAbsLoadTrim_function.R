#
#
#
# Function to load and trim abs from Aqualog
# 29june2015
########

ABStrim <- function(graphheadings, samplewd) {
  
  # change wd to file with all of the abs files in it
  setwd(samplewd)
  
  #read in abs file
  #Absfile <- test2[i,3] # set Abs file for the sample
  
  #Abs filename
  Absfilename <- toString(graphheadings[i,3]) 
  # Read in particular abs filename 
  Absfile <- read.delim(Absfilename, header= FALSE, sep = "")
  
  ############## trim and arrange the abs files
  ######## trim abs file
  # get excitation wavelengths. Contained in first column of data
  row.names(Absfile) <- Absfile[,1]
  
  #Abs <- data.frame(Absfile[,-1])
  Absfile1 = data.frame(t(Absfile))
  Absfile2 = Absfile1[-1,] #remove wavelengths from dataset
  
  return(Absfile2)
}