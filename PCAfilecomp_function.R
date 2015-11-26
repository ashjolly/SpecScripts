# function for compiling data for analysis through PCA
# DBP project
# Ashlee J's PhD 
# 25Nov2015
#########

PCA.eem <- function(filelist, directory){
  ################################################################################
  # install libraries necessary for analysis
  library('gsubfn')
  library('abind')
  library('zoo')
  library("EEM")
  library("reshape")
  library('plyr')
  # read all of the almost corrected EEMS into a dataframe, where the third dimention = the number of samples
  setwd(directory)
  EEM.IFE.rm = lapply(filelist, function(x) read.delim(x, header = TRUE, sep = ","))
  
  #################################
  # Compile and decompose EEM in array such that ex-em pairs are the columns and the sample ID is the row prior to pCA
  EEM.row = data.frame(matrix(vector(), 140500, length(EEM.IFE.rm)))
  
  #Solution using melt command, noting that the loop takes a long time.
  for (i in 1:length(EEM.IFE.rm)){
    
    test <- as.data.frame(EEM.IFE.rm[i])
    
    # transpose the data  
    eem.long <- melt(test, id.vars = NULL)
    
    # get variable for the row and column names
    Ex.nm <- colnames(test)
    Em.nm <- rownames(test)
    
    # create variable
    eem.long$exem <- paste(rep(Ex.nm,each=dim(test)[1]), rep(Em.nm,dim(test)[2]), sep = "_")
    
    #eem.long <- arrange(eem.long,Em.nm) arranges columns by the em wavelengths
    
    # append by rows to the dataframe that you've already created
    EEM.row[,i] <- as.numeric(eem.long[,2])
    rownames(EEM.row) <- eem.long$exem
  }
  
  # add the sample ID from the graph headings file as column names
  sample.ID <- lapply(strsplit(filelist, "_"), function(x) x[1])
  colnames(EEM.row) <- sample.ID
  # take transpose
  EEM.final <- as.data.frame(t(EEM.row))
  
  return(EEM.final)
}


