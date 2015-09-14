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
data.original <- as.data.frame(read.csv(paste(directory, 'codestats_heatmap.csv', sep = ''), sep = ",", header = TRUE))

# take out the 'all' categories from the data
data <- subset(data.original, data.original$Status != "All")

# Make two new columns that contain the Main, Sub and status of the policy areas
data$Policy_response <- paste(data$Main_Policy, data$Sub_Policy, data$Status, sep="_") # merge first three columns into third column
data$Policy_only <- paste(data$Main_Policy, data$Sub_Policy, sep="_")

## Normalize data by the total number of responses within a policy area
# First, by the main policy area
# group by policy group + count total number of responses
# According to the seven main policy areas
total.responses <- aggregate(data[,5:23], by=list(Policy = data$Main_Policy), FUN = sum)

# According to the sub policy areas within the larger policy area
subpolicy.responses <- aggregate(data[,5:23], by=list(Policy = data$Policy_only), FUN = sum)

main.policyareas <- unique(data$Main_Policy)    # find the total number of unique main policy areas
sub.policyareas <- unique(data$Policy_only)     # find the total number of unique sub piolicy areas

stakeholders <- unique(colnames(data[,5:23]))   # find the total number of stakeholder groups

#############################################
# normalize by the total number of responses within a category
# create a new dataframe with the normalized responses
# Note that you want to normalize according to the sub policy area, to get what percent wanted a certain level of legislation

j = length(sub.policyareas)

for(i in 1:j){
  
  polarea.temp <- sub.policyareas[i]   # get specific policy area
  
  #get sum of each column for the policy area
  temp.sum <- subpolicy.responses[subpolicy.responses$Policy == polarea.temp,]

  # take subset of data according to policy area
  temp.subset <- subset(data, data$Policy_only == polarea.temp)
  
  # normalize subset across columns in subset of data
  temp.normalize <- as.data.frame(sweep(as.matrix(temp.subset[,5:23]), 2, as.matrix(temp.sum[,2:20]), "/"))
  
  # bring back in Policy and status columns
  temp.normalize$Policy_only <- temp.subset$Policy_only
  temp.normalize$Main_Policy <- temp.subset$Main_Policy 
  temp.normalize$Sub_Policy <- temp.subset$Sub_Policy 
  temp.normalize$Status <- temp.subset$Status 
  temp.normalize$Policy_response <- temp.subset$Policy_response
    
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
###### preprocessing for all data. 
# This is to create a heat map from all of the policy groups according to the number of submissions

# find row names
rnames <- data$Policy_response
mat.data <- data.matrix(data[,5:22])          # convert data to matrix
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
rounded.norm<- format(round(normalized[,1:19], 1), nsmall = 1)      # ensure that only have 1 decimal places
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

increase = subset(normalized, normalized$Status == 'More, stronger regulation')
increase.factor = data.frame(increase[,1:19]*101, increase[,20:24])

maintain = subset(normalized, normalized$Status == 'Moderate regulation')
maintain.factor= data.frame(maintain[,1:19]*51, maintain[,20:24])

decrease = subset(normalized, normalized$Status == 'Less, weaker regulation')
decrease.factor= data.frame(decrease[,1:19]*1, decrease[,20:24])

factor <- rbind(increase.factor, decrease.factor, maintain.factor) # bind all three subsets back together

# sum all three groupings together by policy area
total.factors<- data.frame(aggregate(factor[,1:19], by=list(Policy = factor$Main_Policy), FUN = sum))

#### graph above
# find row names
rnames <- total.factors$Policy
rounded.factors<- format(round(total.factors[,2:19], 1), nsmall = 1)      # ensure that only have 2 decimal places
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
box.data <- data.frame(total.factors[,2:19])
png(paste(directory, "WSA_normbox.png", sep = ""),    # create PNG for the heat map        
    width = 7*300,        # 5 x 300 pixels
    height = 7*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

boxplot.matrix(as.matrix(box.data), use.cols = TRUE, 
               main=toupper("Distribution - Weighted Response Factor"), font.main=10, 
               cex.main=1.2, xlab="Stakeholder", ylab="Response factor", font.lab=10, 
               col="darkgreen", las=3, margins =c(8,50)
               )

# overlay where the act itself is in the boxplots?

############### influence maps
# Compare what was placed in the Act to pespectives of different groups

##### Decision prediction
# Try and use decision prediction
# Assumes that decision made is a combination of different inputs

