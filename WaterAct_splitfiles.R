
# Script for choosing random individual files for Water Act paper

###########

## set working directory
rm(list = ls())
ls()

library(reshape)
library(plyr)
library(gsubfn)

#######

# create a 

directory <-"/Users/ashlee/Documents/WaterAct Paper/test_splitfiles" #where all of the split PDFs are located

setwd(directory)

## Create a vector with all of the data files and the numbers they are associated with
filelist_pdf <- list.files(pattern = ".pdf$")

# for getting 10% of individual submissions
numind = length(filelist_pdf)

per.10 =round((numind*0.1),digits = 0)

#create a vector of randomly generated numbers that comprise 10% of the total number of submisisons, but spans to the number of individual submissions
random.10 <- sample(1:numind, per.10, replace=F)

#select file numbers that correspond to randomly generated vector
sample.ID <- 0 #create sample ID variable

for (i in 1:numind){
  sample.ID.temp <- strapplyc(filelist_pdf[i], "_(.*).pdf", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}

filelist_pdf.1 <- as.data.frame(cbind(filelist_pdf, sample.ID))


# take a random sample of 10% from a dataset mydata 
# sample without replacement
mysample <- filelist_pdf.1[sample(1:nrow(filelist_pdf.1), per.10,
                          replace=FALSE),]
