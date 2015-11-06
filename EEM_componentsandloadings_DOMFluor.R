##################################
# Components and loadings for results from DOM FLuor
# DOM Fluor results
# Results from DOM Fluor output as Fmax values

# Input

## set working directory
rm(list = ls())
ls()

#directoryDomFluor <-"/Users/ashlee/Documents/MATLAB/DOMFluor/DBP_pre" 
#directoryDomFluor <-"/Users/ashlee/Documents/MATLAB/DOMFluor/WL" 

directory_DBPpre <- "/Users/user/Documents/MATLAB/toolbox/PARAFACresults/DOMFluor_DREEMSHydbrid/DBPPre" 
setwd(directory_DBPpre)

####
# input Fmax file
fpath <- file.path(directory_DBPpre, paste("FMax.csv", sep = ""))
top <- c("C1", "C2", "C3", "C4", "C5", "C6")
FMax = read.csv(fpath, header= FALSE, sep = ",", col.names = top)

fpath <- file.path(directory_DBPpre, paste("01key.csv", sep = ""))
sample.ID <- read.csv(fpath, header= FALSE, sep = ",")

FMax <- cbind(sample.ID, FMax)

# sum rows and express all of the Fmax values as Fmax/sum(Fmax)

n <- dim(FMax)[1]
j <- dim(FMax)[2]

for (i in 1:n){
  Fsum <- sum(FMax[i,c(2:j)])
  #FMax[i,c(j+1)] <- Fsum
  for (k in 1:j){
    Fper <- FMax[i,k]/Fsum*100
    FMax[i,c(j+k)] <- Fper
  } 
}

#change 1 to 2 if you have sample IDs

# export Fmax/Fsum *100 values
len <- dim(FMax)[2]
#data.out <- FMax[,c(1,(j+2):len)]
data.out <- FMax[,c(1,(j+1):len)]

#save
write.table(data.out, file = paste(directory_DBPpre, "/percentloadings.csv", sep = ""), row.names = FALSE,col.names = FALSE, sep = ",")

################## end