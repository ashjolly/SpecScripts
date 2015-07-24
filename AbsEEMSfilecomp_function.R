#
# function for creating master file for corrected absorbance and fluorescence files
# For passing these files into a file that will calculate absorbance and fluorescence indicies
# 24 July 2015
# AJ PhD project
##################

abseemfilecomp <- function(directoryAbsEEMs, projectname){
  
  setwd(directoryAbsEEMs)
  #create column with sample ID - extracted from corrected EEMS filename
  filelist_EEMScor <- list.files(pattern = "_Corrected.csv$")
  y = length(filelist_EEMScor)
  
  sample.ID <- 0 #create sample ID variable
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_EEMScor[i], paste("(.*)", "_", projectname, "_Corrected.csv", sep = ""), simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_EEM <- cbind(filelist_EEMScor , sample.ID)
  
  ############ create column with sample ID - extracted from corrected Abs filename
  #Abs 
  filelist_Abscor <- list.files(pattern = "_AbsCorrected.csv$")
  
  #create column with sample ID - extracted from ABS filename
  y = length(filelist_Abscor)
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_Abscor[i], paste("(.*)","_", projectname,"_AbsCorrected",".csv", sep = ""), simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_Abs <- cbind(filelist_Abscor, sample.ID)
  
  #######
  # Merge EEM and Abs filenames by sample ID to create file with all of the filenames
  data.1 <- merge(filelist_EEM, filelist_Abs,  by = "sample.ID", all = TRUE)
  
  return(data.1)
}
