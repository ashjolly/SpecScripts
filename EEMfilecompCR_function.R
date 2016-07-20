# function for reading in Abs, blank and EEM files from the Aqualog so that they can be corrected and to calculate spectral parameters
# 22June2015 Ashlee
#############

EEMfilecompCR <- function(workdir, dil, EEMfiletype) {
  
  setwd(workdir)
  #above directory contains all blank, Abs and EEM files for correction and calculation from Aqualog
  
  #########
  #Blank files
  filelist_Blank <- list.files(pattern = "_BEM.dat$")
  
  #create column with sample ID - extracted from blank filename
  y = length(filelist_Blank)
  
  sample.ID <- 0 #create sample ID variable
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_Blank[i], "(.*)Blank_BEM", simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_Blank <- unique(cbind(filelist_Blank, sample.ID))
  
  ###########
  #Abs 
  filelist_Abs <- list.files(pattern = "_ABS.dat$")
  #create column with sample ID - extracted from ABS filename
  y = length(filelist_Abs)
  sample.ID <- 0 #create sample ID variable
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_Abs[i], "(.*)ABS_ABS", simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_Abs <- unique(cbind(filelist_Abs, sample.ID))
  
  #########
  # raw EEMS files - note that these are IFM and RM
  #filelist_EEMS <- list.files(pattern = "PEM.dat$")
  
  #below is raw eem without any corrections
  filelist_EEMS <- list.files(pattern = "_EEM")
  
  #create column with sample ID - extracted from EEMS filename
  y = length(filelist_EEMS)
  sample.ID <- 0 #create sample ID variable
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_EEMS[i], paste("(.*)IFERM", sep = ""), simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_EEMS <- unique(cbind(filelist_EEMS, sample.ID))
  
  #######
  # Merge blank, EEM, abs and dilution files according to sample ID
  #alter so that it mearges according to sample ID, which is contained
  data.3 <- Reduce(function(x, y) merge(x, y, by = "sample.ID", all=TRUE), list(filelist_EEMS, filelist_Abs, filelist_Blank, dil))
  
  return(data.3)
}