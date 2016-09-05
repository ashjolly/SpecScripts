#######################
# Components and loadings script
# 21april2015
# Adapted from components and loadings spreadsheet from CM workshop

# Note that the CM PARAFAC code will output three text files, A B and C, containing PARAFAC results.
# This is in addition to residials in csv files.

# Eventually also add to the code to analyze these residuals to flag when residula EEMS are too high?
# Code will prodcue two files: 1) with FMax values; 2) one file with Fmax/sum(Fmax) across all components
# See Murphy, K. R., Stedmon, C. A., Graeber, D., & Bro, R. (2013). Fluorescence spectroscopy and multi-way techniques. PARAFAC. Analytical Methods, 5(23), 6557–11. http://doi.org/10.1039/c3ay41160e
# pg 6565
# Fmax = a*maxB*maxC
# a = model score
# C = loadings for excitation
# B = loadings for emission
#################################

## set working directory
rm(list = ls())
ls()

# First, save the A B and C text files from PARAFAC in a specific folder.
#Directory
#DBP
#directoryCMresults <-"/Users/user/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_CMPARAFAC" 
#directoryCMresults <-"/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_pre_CM_PARAFAC" 
#directoryCMresults <-"/Users/user/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_delta/DBP_delta_CMPARAFACresults"
#directoryCMresults <- "/Users/user/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_all_corrected/DBP_all_CMPARAFACresults"
#directoryCMresults <- "/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_preandpost/DBP_prepost_CMresults"
# WL
#directoryCMresults <-"/Users/user/Dropbox/PhD Work/PhD Data/WL_data/WL_Fluorescence/WL_CMParafac_Results"
# CR Results
#directoryCMresults <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR EEMS/CR_CMPARAFACResults"
# CR soil results
directoryCMresults <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_soilextracts/CR_soilextracts_fluorescence/CR_PARAFACCMresults"
# CR Soil Lysimeters
#directoryCMresults <- "/Users/user/Dropbox/PhD Work/PhD Data/CR_Data/CR_Lysimeter/CRLys_Fluorescence/CRLys_CMPARAFACResults"

setwd(directoryCMresults) 
#sample type
#sample.type <- "DBP"
#sample.type <- "WL"
#sample.type <- "DBPpre"
#sample.type <- "DBPpost"
#sample.type <- "DBPdelta"
#sample.type <- "DBPall"
#sample.type <- "DBPprepost"
#sample.type <- "CR"
sample.type <- "SoilExtracts"
#sample.type <- "CRLys"

# Read in the A file
#afile = "A_DBPpostandpre.txt"
Apath <- file.path(directoryCMresults, paste("A_", sample.type, ".txt", sep = ""))
A = as.data.frame(read.table(Apath, header= FALSE, sep = ""))

#Read in the B file
bpath <- file.path(directoryCMresults, paste("B_", sample.type, ".txt", sep = ""))
B = read.table(bpath, header= FALSE, sep = "")
# add in wavelengths for B (emission wavelengths for CM)
wave <- seq(350, 550, by = 2)
rownames(B) <- wave

#Read in the C file
cpath <- file.path(directoryCMresults, paste("C_", sample.type, ".txt", sep = ""))
C = read.table(cpath, header= FALSE, sep = "")
# add in wavelengths for the c - excitation for CM
wave <- seq(250,400, by = 5)
rownames(C) <- wave

#read in the graph headings file used to identify the samples in the A file
gH <- paste("GraphHeadings_", sample.type, ".txt", sep = "")
gHpath <- file.path("/Users/user/Documents/MATLAB/CM_graphheadings", paste(gH, sep = ""))
sample.ID <- read.table(gHpath, header= FALSE, sep = "")
sample.ID <- subset(sample.ID, !duplicated(sample.ID))

##########
# Find the fMax factors from the b and c files
# find max for each column (component) in c file
cmax = apply(C,2,max)
# find max for each column (component) in b file
bmax = apply(B,2,max)

#create empty row for Fmax calculation
temprow <- matrix(c(rep.int(NA,length(cmax))),nrow=1,ncol=length(cmax))
# make it a data.frame and give cols the same names as data
Fmax <- data.frame(temprow)
remove(temprow)

#Fmax factor = cmax*bmax for each of the 13 components
j = length(Fmax)

for (i in 1:j){
  Fmax.temp <- as.numeric(bmax[i])*as.numeric(cmax[i])
  #Fmax.temp <- as.numeric(Fmax.temp, stringsAsFactors=FALSE)
  Fmax[i] <- Fmax.temp
}

Fmax.all <- rbind(bmax, cmax, Fmax)
rownames(Fmax.all) <- c("bmax", "cmax", "Fmax")
Fmax.all <- as.data.frame(Fmax.all)

#################
# calculate loadings
# where loadings = A*Max factor for that component

#create empty loadings dataframe
temprow <- matrix(c(rep.int(NA,length(Fmax))),nrow=(dim(A)[1]),ncol=(dim(A)[2]))
# make it a data.frame and give cols the same names as data
loadings <- data.frame(temprow)
remove(temprow)

# Multiply each column by corresponding fmax file
n <- length(Fmax)
g = dim(A)[1]

for (i in 1:g){
  for (j in 1:n){
    loadings[i,j] <- A[i,j]*Fmax[j]
  }
}

# bind sample id
Fmax.CM <- cbind(sample.ID, loadings)
colnames(Fmax.CM) <- cbind("sample.ID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13")

# write file with Fmax values
PARAFACpath <- file.path(directoryCMresults, paste(sample.type,"_componentsandloadings_CM_Fmax", ".csv",sep = ""))
write.table(Fmax.CM, file = PARAFACpath, row.names = FALSE,col.names = TRUE, sep = ",")

### plot the ex/em loadings for each component and save within the PARAFAC folder
for (i in 1:13){
  png(filename=paste(directoryCMresults, "/CMC", i, "_Loadings.png", sep = ""))
  plot(rownames(C), C[,i], type = "l", main = paste("C", i, "- CM", sep = ""), xlim =c(250, 550))
  lines(rownames(B), B[,i], col = "blue")
  legend('topright', c('ex','em'), lty=c(1,1), lwd=c(2.5,2.5), col=c('black','blue'), cex = 0.3)
  dev.off()
}

################## Fmax/sumFmax (percent)
#take the sum of each row in the loadings
sums <- rowSums(loadings)

# find the percent of each, where percent = loading for each component/sum for each sample
num.samples <- dim(loadings)[1] #sample = row

per.loadings <- function(loading, sum) {
  percent <- loading / sum *100
  return (percent)
}
#create empty dataframe for new variable
temprow <- matrix(c(rep.int(NA,length(loadings))),nrow=(dim(loadings)[1]),ncol=(dim(loadings)[2]))
# make it a data.frame and give cols the same names as data
percent <- data.frame(temprow)
remove(temprow)

for (i in 1:num.samples){
    # i = sample ID
  percent.temp <- per.loadings(loadings[i,], sum = sums[i])
  percent[i,] <- percent.temp
    }
  
# bind sample id
PARAFAC.CM <- cbind(sample.ID, percent)
colnames(PARAFAC.CM) <- cbind("sample.ID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13")

#calculate the percent protein - C8 plus C13
PARAFAC.CM$perprotein <- PARAFAC.CM$C8 + PARAFAC.CM$C13

#calculate the redox index
PARAFAC.CM$redox <- (PARAFAC.CM$C4+PARAFAC.CM$C5+PARAFAC.CM$C7+PARAFAC.CM$C9)/(PARAFAC.CM$C2+PARAFAC.CM$C4+PARAFAC.CM$C5+PARAFAC.CM$C7+PARAFAC.CM$C9+PARAFAC.CM$C11+PARAFAC.CM$C12)

# Save file to directory with PARAFAC results
PARAFACpath <- file.path(directoryCMresults, paste(sample.type,"_componentsandloadings_CM", ".csv",sep = ""))
write.table(PARAFAC.CM, file = PARAFACpath, row.names = FALSE,col.names = TRUE, sep = ",")
