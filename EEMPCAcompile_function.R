# function for placing eems into a format suitable for PCA analysis
# Sample ID as rows, and exitation and emission pairs as columns
# Ashlee JOllymore's PhD
# 6Oct2015
#################

PCA.EEM <- function(EEM){
  
  # Compile and decompose EEM such that ex-em pairs are the columns and the sample ID is the row
  for (j in 1:as.numeric(dim(EEM)[2])){
    temp.PCA = data.frame(t(EEM[,j]))
    colnames(temp.PCA) = (paste(colnames(EEM[j]), row.names(EEM), sep = '_'))
    
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
    
    remove(temp.PCA) # remove the temporary PCA dataset
  }
  
  return(PCA.data)
}