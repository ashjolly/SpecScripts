#
#
# WSA Analysis script
# For analyzing results from NVivo
# 15Aug2015
# Kiely M and Ashlee J
# references: http://www.r-bloggers.com/search/heatmap
# http://sebastianraschka.com/Articles/heatmaps_in_r.html
################

## set working directory
rm(list = ls())
ls()
directory <- "/Users/ashlee/Documents/WaterAct Paper/"

# libraries
library('gplots')
library('plyr')
library("d3heatmap")


## load data
data <- as.data.frame(read.csv(paste(directory, 'test_heatmap.csv', sep = ''), sep = ",", header = TRUE))

# group by policy group + count total number of responses
total.responses <- aggregate(data[,3:18], by=list(Policy = data$Policy), FUN = sum)
total.responses$status <- 'total'
data1 <- rbind(data, total.responses)

policy.areas <- unique(data1$Policy)
# Normalize each entry by the total number of responses
test <- subset(data1, data1$Policy == policy.areas[1])
normalize <- test[c(1:3), c(3:18)]/test[4, c(3:18)]

# normalize by ddply
#f <- function {
#  # get total responses for each
#  aggregate(data[,3:18], by = by=list(Policy = data$Policy))
  
#}
#normalize <- apply(data[,3:18], 2, function)

# merge first two columns into third column
data$Policy_response <- paste(data$Policy, data$status, sep="_")

###### preprocessing
# find row names
rnames <- data$Policy_response
mat.data <- data.matrix(data[,3:18])          # convert data to matrix
rownames(mat.data) <- rnames                  # assign row names 

# do heat map
# create colour palette
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# creates a 5 x 5 inch image
png(paste(directory, "heatmaps_in_r.png", sep = ""),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "test- WSA paper 
          (non-normalized data)", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(8,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 0.6,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering

####### second method
# interactive heatmap
d3heatmap(mat.data, scale = "column", 
          dendrogram = "row",
          col=my_palette,
          margins =c(8,50), 
          cexRow=0.6, 
          cexCol = 0.35,
          key = TRUE
          )

# categorize according to the
# allocate weighting to each categorization

