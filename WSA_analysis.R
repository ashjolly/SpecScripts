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

## Normalize data by the total number of responses within a group

# group by policy group + count total number of responses
total.responses <- aggregate(data[,3:18], by=list(Policy = data$Policy), FUN = sum)
total.responses$status <- 'total'

data$Policy_response <- paste(data$Policy, data$status, sep="_") # merge first two columns into third column

policy.areas <- unique(data$Policy)      # find the total number of unique policy areas
stakeholders <- unique(colnames(data[,3:18]))

#############################################
# normalize by the total number of responses within a category
# create a new dataframe with the normalized responses

j = length(policy.areas)

for(i in 1:j){
  
  polarea.temp <- policy.areas[i]   # get specific policy area
  
  #get sum of each column for the policy area
  temp.sum <- total.responses[total.responses$Policy == polarea.temp,]

  # take subset of data according to policy area
  temp.subset <- subset(data, data$Policy == polarea.temp)
  
  # normalize subset across columns in subset of data
  temp.normalize <- as.data.frame(sweep(as.matrix(temp.subset[,3:18]), 2, as.matrix(temp.sum[,2:17]), "/"))
  
  # bring back in Policy and status columns
  temp.normalize$Policy <- temp.subset$Policy
  temp.normalize$status <- temp.subset$status
  temp.normalize$Policy_response <- paste(temp.normalize$Policy, temp.normalize$status, sep = "_")
  
  # create a new dataset where the normalized data together by rows
  # if the merged dataset doesn't exist, create it
  # if the normalized dataset does exist, append to it
  
  if (exists("normalized")){
    normalized<-rbind(normalized, temp.normalize)
    rm(temp.normalize)
  }
  
  if (!exists("normalized")){
    normalized <- temp.normalize
  }
  
} # end of loop

# replace Nas with 0's
normalized[is.na(normalized)] = 0

####################### Graphing
###### preprocessing for all data
# find row names
rnames <- data$Policy_response
mat.data <- data.matrix(data[,3:18])          # convert data to matrix
rownames(mat.data) <- rnames                  # assign row names 

# do heat map
# create colour palette
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# creates a 5 x 5 inch image
png(paste(directory, "WSA_heatmap_policyfocus_col.png", sep = ""),    # create PNG for the heat map        
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

#################
# graph according to normalized data
# preprocessing

# find row names
rnames <- normalized$Policy_response
rounded.norm<- format(round(normalized[,1:16], 1), nsmall = 1)      # ensure that only have 1 decimal places
mat.data2 <- data.matrix(rounded.norm)          # convert data to matrix
rownames(mat.data2) <- rnames                         # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normalized.png", sep = ""),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data2,
          cellnote = mat.data2,  # same data set for cell labels
          main = "test- WSA paper 
          (normalized data)", # heat map title
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
d3heatmap(mat.data2, scale = "column", 
          dendrogram = "row",
          col=my_palette,
          margins =c(8,50), 
          cexRow=0.6, 
          cexCol = 0.35,
          key = TRUE
)

#####################
# Categorize according to allocate weighting to each categorization
# group according to each policy area. Multiply the number of respondents by..

# final score (for stakeholder group for policy area) = 
# (# respondents for increase/total number of respondents)*3 
# + (# respondents for maintain/total number of respondents)*2
# + (# respondents for decrease/total number of respondents)*1

#multiply each 'increase' segment by 3 in normalized

increase = subset(normalized, normalized$status == 'Increase, strengthen regulation')
increase.factor = data.frame(increase[,1:16]*3, increase[,17:19])

decrease = subset(normalized, normalized$status == 'Decrease, weaken regulation')
decrease.factor= data.frame(decrease[,1:16]*1, decrease[,17:19])

maintain = subset(normalized, normalized$status == 'Maintain status quo, minor changes')
maintain.factor= data.frame(maintain[,1:16]*2, maintain[,17:19])

factor <- rbind(increase.factor, decrease.factor, maintain.factor) # bind all three subsets back together

# sum all three groupings together by policy area
total.factors<- data.frame(aggregate(factor[,1:16], by=list(Policy = factor$Policy), FUN = sum))

#### graph above
# find row names
rnames <- total.factors$Policy
rounded.factors<- format(round(total.factors[,2:17], 1), nsmall = 1)      # ensure that only have 2 decimal places
mat.data3 <- data.matrix(rounded.factors)            # convert data to matrix
rownames(mat.data3) <- rnames                             # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normfactors.png", sep = ""),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data3,
          cellnote = mat.data3,  # same data set for cell labels
          main = "test- Policy Areas Weighting Factors", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(8,16),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 0.6,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering

####### second method
# interactive heatmap
d3heatmap(mat.data3, scale = "column", 
          dendrogram = "column",
          col=my_palette,
          margins =c(8,50), 
          cexRow=0.6, 
          cexCol = 0.5,
          key = TRUE,
          main = "test- Policy Areas Weighting Factors"
)


##### box plots - to answer   question of whether stakeholder groups called for more regulation uniformly along policy areas
box.data <- data.frame(total.factors[,2:17])
boxplot.matrix(as.matrix(box.data), use.cols = TRUE)

