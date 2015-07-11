#
# Script for creating delta DBP - from abs files
# The aim of this script is to get the abs of the prechlorination minus the post chlorination
# DBP project, ashlee's phd
# 10July2015
# ashlee jollymore
#############

## set working directory
rm(list = ls())
ls()

## necessary packages
library(reshape)
library(plyr)
library(gsubfn)

######
# load pre and post chlorination files, match according to sample ID, and then minus pre - post

# Where abs files are located
prechlor.files <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_all"
postchlor.files <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_all"

# directory where the delta abs indicies will be saved. No need to save spectra
delta.abs <- "/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_deltaabs"

# read in pre chlor abs files
# prechlor files
setwd(prechlor.files)
filelist_pre <- list.files(pattern = "ABS.dat$")

# Get sample ID from pre chlor files
sample.ID <- 0 #create sample ID variable
samplecode <- 0
n = length(filelist_pre)

for (i in 1:n){
  sampleID.pre.temp <- strapplyc(filelist_pre[i], "DBP(.*)ABS", simplify = TRUE)
  sample.ID[i] <- sampleID.pre.temp
  
  samplecode.pre.temp <- strapplyc(filelist_pre[i], "001(.*)ABS", simplify = TRUE) 
  samplecode[i] <- samplecode.pre.temp 
}

filelist.pre <- as.data.frame(cbind(filelist_pre, sample.ID, samplecode))

# get dilution factor
top = c("samplecode", "dilutionfactorpre")
dilution.pre <- as.data.frame(read.csv("/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_prechlorination/DBP_prechlor_Aqualogdilution.csv", 
                                       sep=",", header = TRUE, col.names = top))

filelist.pre <- merge(filelist.pre, dilution.pre, by = "samplecode")

# Get sample ID from post chlor files
#postchlor
setwd(postchlor.files)
filelist.post <- list.files(pattern = "ABS.dat$")

n = length(filelist_post)
for (i in 1:n){
  sampleID.post.temp <- strapplyc(filelist_post[i], "DBPChlor(.*)ABS", simplify = TRUE)
  sample.ID[i] <- sampleID.post.temp 
  
  samplecode.pre.temp <- strapplyc(filelist_post[i], "001(.*)ABS", simplify = TRUE) 
  samplecode[i] <- samplecode.pre.temp 
}

filelist.post <- as.data.frame(cbind(filelist.post, sample.ID, samplecode))

# need dilution factors
top = c("samplecode", "dilutionfactorpost")
dilution.post <-as.data.frame(read.csv("/Users/ashlee/Documents/UBC Data/DBP_data/DBP_fluorescence/DBP_postchlorination/DBP_postchlor_Aqualogdilution.csv", 
                                       sep=",", header = TRUE, col.names = top))
filelist.post <- merge(filelist.post, dilution.post, by = "samplecode")

# merge files by sample ID
all <- merge(filelist.post, filelist.pre, by = "sample.ID", all = FALSE)

# create abs ind variable for loop
abs.ind = data.frame(matrix(vector(), 0, 8))

## loop to calculate abs pre - abs post

n = dim(all)[1]
for (i in 1:n){
  # read in abs from pre and post files
  
  ############# pre files
  setwd(prechlor.files)
  pre.samplename <- toString(all[i,6])
  pre.temp <- as.data.frame(read.delim(pre.samplename, header= FALSE, sep = "", stringsAsFactors=FALSE))

  # account for dilution factor
  pre.dil <- as.data.frame(pre.temp[,2] * all[i,7])
  rownames(pre.dil) <- pre.temp[,1]
  
  ############# post files
  setwd(postchlor.files)
  post.samplename <- toString(all[i,3])
  post.temp <- as.data.frame(read.delim(post.samplename, header= FALSE, sep = "", stringsAsFactors=FALSE))
  
  # account for dilution factor
  post.dil <- as.data.frame(post.temp[,2] * all[i,4])
  rownames(post.dil) <- post.temp[,1]
  
  ############# trim
  #trim so that exitation and emission goes from the same
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("Abstrim_DBP_function.R")
  
  pre.trim <- abstrim(abs = pre.dil, minex = 240)
  post.trim <- abstrim(abs = post.dil, minex = 240)
  
  ############# calculate pre - post
  delta.abs <- pre.trim - post.trim
  row.names(delta.abs) <- paste("X",seq(800, minex, by = -2), sep = "")
  # convert to transpose
  delta.abs = as.data.frame(t(delta.abs))
  
  ############# calculate absorbance indicies from delta spectra
  # call function
  setwd("/Users/ashlee/SpecScripts") 
  source("Aqualog_Absindicies_v1.R")
  
  #call the function to calculate indicies
  Abs.ind <- Abs(absorbance = delta.abs)
  
  
