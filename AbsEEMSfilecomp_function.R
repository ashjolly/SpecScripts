#
# function for creating master file for corrected absorbance and fluorescence files
# For passing these files into a file that will calculate absorbance and fluorescence indicies
# 24 July 2015
# AJ PhD project
##################

abseemfilecomp <- function(directoryAbsEEMs, directoryRaleigh, projectname, directorynoncorabs, filelist_EEMScor){
  
  setwd(directoryAbsEEMs)
  #create column with sample ID - extracted from corrected EEMS filename
  
  y = length(filelist_EEMScor)
  
  sample.ID <- 0 #create sample ID variable
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_EEMScor[i], paste("(.*)", "_", projectname, "_Raleighcorr.csv", sep = ""), simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist <- cbind(filelist_EEMScor , sample.ID)
  
  ###########
  #Abs - non-corrected
  setwd(directorynoncorabs)
  filelist_Abs_noncorr <- list.files(pattern = "ABS.dat$")
  #create column with sample ID - extracted from ABS filename
  y = length(filelist_Abs_noncorr)
  
  for (i in 1:y){
    sample.ID.temp <- strapplyc(filelist_Abs_noncorr[i], "001(.*)ABS", simplify = TRUE)
    sample.ID[i] <- sample.ID.temp
  }
  filelist_Abs_noncorr <- cbind(filelist_Abs_noncorr, sample.ID)
  
  ############ create column with sample ID - extracted from corrected Abs filename
  #Abs 
  setwd(directoryAbsEEMs)
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
  data.1 <- merge(filelist, filelist_Abs,  by = "sample.ID", all = TRUE)
  data.1 <- merge(data.1, filelist_Abs_noncorr, by = "sample.ID", all = TRUE)
  return(data.1)
}
