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
  data.4 <- unique(data.3)
  
  ## Data massaging!
  # for samples 286 - 344 (8jan2014), use abs from Jan 20, 2014 (1339-1359)
  data.4[145:199,4] <- "CR1339Blank_BEM.dat"
  # for sample 1031, use blank 1031
  data.4[(which(data.4$sample.ID == "CR1031")), 4] <- "CR1032Blank_BEM.dat"
  # For sample755, use blank 756
  data.4[(which(data.4$sample.ID == "CR756")), 4] <- "CR755Blank_BEM.dat"
  # For sample 840- 845, use blank 839
  data.4[(which(data.4$sample.ID == "CR840"):which(data.4$sample.ID == "CR845")), 4] <- "CR839Blank_BEM.dat"
  # For sample 851, use blank 850
  data.4[(which(data.4$sample.ID == "CR851")), 4] <- "CR850Blank_BEM.dat"
  
  # only return the data that has an EEMs associatedwith it.
  data.5 <- data.4[complete.cases(data.4$filelist_EEMS),]
  
  return(data.5)
}