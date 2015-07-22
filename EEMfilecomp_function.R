# function for reading in Abs, blank and EEM files from the Aqualog so that they can be corrected and to calculate spectral parameters
# 22June2015 Ashlee
#############

EEMfilecomp <- function(workdir, dil, EEMfiletype) {
  
  setwd(workdir)
  #above directory contains all blank, Abs and EEM files for correction and calculation from Aqualog
  
  #########
  #Blank files
  filelist_Blank <- list.files(pattern = "BEM.dat$")
  
  #create column with sample ID - extracted from blank filename
  y = length(filelist_Blank)
  
  sample.ID <- 0 #create sample ID variable
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_Blank[i], "001(.*)BEM", simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_Blank <- cbind(filelist_Blank, sample.ID)
  
  ###########
  #Abs 
  filelist_Abs <- list.files(pattern = "ABS.dat$")
  #create column with sample ID - extracted from ABS filename
  y = length(filelist_Abs)
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_Abs[i], "001(.*)ABS", simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_Abs <- cbind(filelist_Abs, sample.ID)
  
  #########
  # raw EEMS files - note that these are IFM and RM
  #filelist_EEMS <- list.files(pattern = "PEM.dat$")
  
  #below is raw eem without any corrections
  filelist_EEMS <- list.files(pattern = EEMfiletype)
  
  #create column with sample ID - extracted from EEMS filename
  y = length(filelist_EEMS)
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_EEMS[i], paste("001(.*)", EEMfiletype, sep = ""), simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_EEMS <- cbind(filelist_EEMS, sample.ID)
  
  #######
  # Merge blank, EEM, abs and dilution files according to sample ID
  #alter so that it mearges according to sample ID, which is contained
  data.1 <- merge(filelist_EEMS, filelist_Abs,  by = "sample.ID", all = TRUE)
  data.2 <- merge(data.1, filelist_Blank, by = "sample.ID", all = TRUE)
  data.3 <- merge(data.2, dil, by = "sample.ID", all = TRUE)
  
  #Remove data.1 and data.2 - only need merged file
  remove(data.1, data.2)
  return(data.3)
}