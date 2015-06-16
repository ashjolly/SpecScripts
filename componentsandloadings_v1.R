#######################
# Components and loadings script
# 21april2015
# Adapted from components and loadings spreadsheet from CM workshop

# Note that the CM PARAFAC code will output three text files, A B and C, containing PARAFAC results.
# This is in addition to residials in csv files.

# Eventually also add to the code to analyze these residuals to flag when residula EEMS are too high
#################################

## set working directory
rm(list = ls())
ls()

# First, save the A B and C text files ffrom PARAFAC in a specific folder.
#Directory
#DBP
directoryCMresults <-"/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_CM_Preandpostmodelled" 
#directoryCMresults <-"/Users/ashlee/Documents/UBC Data/WL_data/WL_Fluorescence/WL_CMParafac_Results/" 

setwd(directoryCMresults) 
#sample type
#sample.type <- "DBP"
#sample.type <- "WL"
sample.type <- "DBP_pre"

# Read in the A file
#afile = "A_DBPpostandpre.txt"
Apath <- file.path(directoryCMresults, paste("A_", sample.type, ".txt", sep = ""))
A = as.data.frame(read.table(Apath, header= FALSE, sep = ""))

#Read in the B file
#bfile = "B_DBPpostandpre.txt"
bpath <- file.path(directoryCMresults, paste("B_", sample.type, ".txt", sep = ""))
B = read.table(bpath, header= FALSE, sep = "")

#Read in the C file
#cfile = "C_DBPpostandpre.txt"
cpath <- file.path(directoryCMresults, paste("C_", sample.type, ".txt", sep = ""))
C = read.table(cpath, header= FALSE, sep = "")

#read in the graph headings file used to identify the samples in the A file
gH <- paste("GraphHeadings_", sample.type, ".txt", sep = "")
gHpath <- file.path("/Users/ashlee/Documents/MATLAB/CorrectedEEMS", paste(gH, sep = ""))
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


n <- length(Fmax)

for (i in 1:n){
  #must convert factors in A!
  loadings[,i] <-  as.numeric(levels(A[,i]))[A[,i]] * as.numeric(Fmax[i])
}

##################
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
