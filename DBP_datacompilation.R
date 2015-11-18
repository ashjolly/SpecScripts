# script for compiling data for analysis 
# through PCA
# DBP project
# #######

# clean up list
rm(list = ls())
ls()
################################################################################
# install libraries necessary for analysis
library('gsubfn')
library('abind')
library('zoo')

############ Fluorescence data
## Prechlorinated EEMS
pre.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_correctedEEMs'

## post Chlorinated EEMS
post.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_correctedEEMs'

# directory for saving data
save.directory <- '/Users/user/Dropbox/PhD Work/PhD Data/DBP_data/DBP_analysisdata'

##############################################################
# functions used in script
#function for transposing and compiling EEMS according to column
PCA.EEM <- function(EEM){
  # Compile and decompose EEM such that ex-em pairs are the columns and the sample ID is the row
  temp.PCA = data.frame(t(EEM))
  colnames(temp.PCA) = (paste(colnames(EEM), row.names(EEM), sep = '_'))
  
  # if the merged dataset does exist, append to it by column
  if (exists("PCA.data")){
    #temp_dataset <- temp.cut
    PCA.data<-cbind(PCA.data, temp.PCA)
    #rm(temp.cut)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("PCA.data")){
    PCA.data <- temp.PCA
  }
}

###############################
# Prechlorinated EEMS
# locate the prechlorinated corrected eems within the file
setwd(pre.directory) 
filelist_DBPpre <- list.files(pattern = "_Corrected.csv$")

# create graph heading variable
graphheadingspre = data.frame((0))

# compile all of the corrected pre chlorination EEMS and  correct from Raman and Raleigh scatter
n = length(filelist_DBPpre)
exmin = 'X240'
project = 'DBPPre'

# run loop over all files within the corrected file list
for (i in 1:n){
  # set working directory back to directory with sample ID + read in EEMs
  setwd(pre.directory) 
  temp.EEMS <- read.delim(filelist_DBPpre[i], header= TRUE, sep = ",")
  
  # ensure that the ex ranges are the same for all of the data - 240 to 200 nm in 2 nm increments
  ex.temp <- colnames(temp.EEMS)
  if(ex.temp[1] != exmin) {
    # if first value in ex.temp is not 240, trim 
    ex.length <- length(ex.temp)
    # find column where the exitation wavelength is 240 to cut from
    x240 = as.numeric(match(exmin,names(temp.EEMS)))
    temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
  } 
  
  # Correct for Raleigh scatter using function - interpolates for first and second Raleigh
  setwd("/Users/user/SpecScripts") 
  source("EEMRaleigh_function.R")
  temp.cut <- raleigh(eem = temp.EEMS, slitwidth1 = 25, slitwidth2 = 25, R1 = "no")
  
  # get the emission variables from the EEM
  em = row.names(temp.cut)
  # get the excitation variables from the cut EEMS
  ex = colnames(temp.cut)
  
  #save as a array
  # if the merged dataset does exist, append to it by column
  if (exists("EEM.dataset")){
    #temp_dataset <- temp.cut
    EEM.dataset <- abind(EEM.dataset, temp.cut, along = 3)
    #rm(temp.cut)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("EEM.dataset")){
    EEM.dataset <- temp.cut
  }
  
  # Create graph headings variable to identify samples
  samplename <- strapplyc(filelist_DBPpre[i], paste("(.*)_", project, "_Corrected", sep = ""), simplify = TRUE)
  graphheadingspre[i,] <-samplename
}

# Compile and decompose EEM in array such that ex-em pairs are the columns and the sample ID is the row prior to pCA
EEM.pre <- EEM.dataset
remove(EEM.dataset)

n = dim(EEM.pre)[3]

# create empty vector
EEM.row = data.frame(matrix(vector(), 5000, 200000))

# assemble using vectorized solution
for (i in 1:n){
  # get the sample data frame
  temp.EEM <- EEM.pre[,,i]
  
  # use function and apply to reorganize for PCA
  EEM.row[i,] <- as.data.frame(apply(temp.EEM, 2, PCA.EEM))
}

# Cut out the NAs inserted when you created a dataframe to put the PCA formatted data into
PCA.pre <- EEM.row[rowSums(is.na(EEM.row)) != ncol(EEM.row),]
PCA.pre <- PCA.pre[,colSums(is.na(PCA.pre)) != nrow(PCA.pre)]
#PCA.pre <- PCA.pre[,complete.cases(PCA.pre)]
#PCA.pre <- PCA.pre[,colSums(is.na(PCA.pre))<nrow(PCA.pre)]

# ensure that the rows are properly named
row.names(PCA.pre) = graphheadingspre[,1] #sample names

# create column names
ex = colnames(temp.cut)
em = row.names(temp.cut)

for (i in 1:length(ex)){
  ex.temp <- t(as.data.frame(paste(ex[i], em, sep="_")))
  
  # if the merged dataset exists, append to it by column
  if (exists("exem")){
    exem <- cbind(exem, ex.temp)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("exem")){
    exem <- ex.temp
  }
}

colnames(PCA.pre) = exem

PCA.pre.cut <- PCA.pre[,c(1:length(exem))]

#### save because vectorizing takes forever!
saveRDS(PCA.pre.cut, paste(save.directory, "/PCApre.rds", sep = ""))
saveRDS(EEM.pre, paste(save.directory, "/EEMpre.rds", sep = ""))
write.table(graphheadingspre, paste(save.directory, "/EEMpresampleID.csv", sep = ""))

######################################################################################################
# post chlorination
# read in the post chlorinated EEMS, correct for Raleigh and assemble for 
# locate the prechlorinated corrected eems within the file
setwd(post.directory) 
filelist_DBPpost <- list.files(pattern = "_Corrected.csv$")
project = 'DBPPost'

# create graph heading variable (or reset to 0s)
graphheadingspost = data.frame((0))

# compile all of the corrected pre chlorination EEMS and  correct from Raman and Raleigh scatter
n = length(filelist_DBPpost)
exmin = 'X240'

# run loop over all files within the corrected file list
for (i in 1:n){
  # set working directory back to directory with sample ID + read in EEMs
  setwd(post.directory) 
  temp.EEMS <- read.delim(filelist_DBPpost[i], header= TRUE, sep = ",")
  
  # ensure that the ex ranges are the same for all of the data - 240 to 200 nm in 2 nm increments
  ex.temp <- colnames(temp.EEMS)
  if(ex.temp[1] != exmin) {
    # if first value in ex.temp is not 240, trim 
    ex.length <- length(ex.temp)
    # find column where the exitation wavelength is 240 to cut from
    x240 = as.numeric(match(exmin,names(temp.EEMS)))
    temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
  } 
  
  # Correct for Raleigh scatter using function - interpolates for first and second Raleigh
  setwd("/Users/user/SpecScripts") 
  source("EEMRaleigh_function.R")
  temp.cut <- raleigh(eem = temp.EEMS, slitwidth1 = 25, slitwidth2 = 25, R1 = "no")
  
  # get the emission variables from the EEM
  em = row.names(temp.cut)
  # get the excitation variables from the cut EEMS
  ex = colnames(temp.cut)
  
  #save as a array
  # if the merged dataset does exist, append to it by column
  if (exists("EEM.dataset")){
    #temp_dataset <- temp.cut
    EEM.dataset <- abind(EEM.dataset, temp.cut, along = 3)
    #rm(temp.cut)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("EEM.dataset")){
    EEM.dataset <- temp.cut
  }
  
  # Create graph headings variable to identify samples
  samplename <- strapplyc(filelist_DBPpost[i], paste("(.*)_", project, "_Corrected", sep = ""), simplify = TRUE)
  graphheadingspost[i,] <-samplename
  
}

# Compile and decompose EEM in array such that ex-em pairs are the columns and the sample ID is the row prior to pCA
EEM.post <- EEM.dataset
remove(EEM.dataset)

n = dim(EEM.post)[3]

# create empty vector
EEM.row = data.frame(matrix(vector(), 5000, 200000))

# assemble using vectorized solution
for (i in 1:n){
  # get the sample data frame
  temp.EEM <- EEM.post[,,i]
  
  # use function and apply to reorganize for PCA
  EEM.row[i,] <- as.data.frame(apply(temp.EEM, 2, PCA.EEM))
}

# Cut out the NAs inserted when you created a dataframe to put the PCA formatted data into
PCA.post <- EEM.row[rowSums(is.na(EEM.row)) != ncol(EEM.row),]
PCA.post <- PCA.post[,colSums(is.na(PCA.post)) != nrow(PCA.post)]

# ensure that the column names are properly named
row.names(PCA.post) = graphheadingspost[,1] #asample names

# create column names
ex = colnames(temp.cut)
em = row.names(temp.cut)

remove(exem)

for (i in 1:length(ex)){
  ex.temp <- t(as.data.frame(paste(ex[i], em, sep="_")))
  
  # if the merged dataset exists, append to it by column
  if (exists("exem")){
    exem <- cbind(exem, ex.temp)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("exem")){
    exem <- ex.temp
  }
}

colnames(PCA.post) = exem
# cut out NA columns

PCA.post.cut <- PCA.post[,c(1:length(exem))]

#### save because vectorizing takes forever!
saveRDS(PCA.post.cut, paste(save.directory, "/PCApost.rds", sep = ""))
saveRDS(EEM.post, paste(save.directory, "/EEMpost.rds", sep = ""))
write.table(graphheadingspost, paste(save.directory, "/EEMpostsampleID.csv", sep = ""))
