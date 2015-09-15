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
col_breaks = c(seq(-1,0,length=100),  # This is for uneven breaking... WILL NEED TO CHANGE..for red
               seq(0,0.8,length=100),              # for yellow
               seq(0.8,1,length=100))              # for green

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
mat.data <- data.matrix(rounded.norm)          # convert data to matrix
rownames(mat.data) <- rnames                         # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normalized.png", sep = ""),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
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

## Second normalization method- normalize data to the total number of respondents within a stakeholder group
# Row at bottom of CSV where Policy_general = Number of submissions
# Aim is to normalize all of the responses within a column by the corresponding total number of submissions within a stakeholder group

# Choose the row with the total nnumber of submissions
no.submissions <- data[117,]

# divide each row by the number of submissions
norm.2 <- data.frame(lapply(data[, 5:22], function(x) x/tail(x,1)))
# add back the row names
norm.2 <- cbind(data[,c(1:3,24:25)], norm.2)

## plot the second way of normalizing data
# find row names
rnames <- norm.2$Policy_response
rounded.norm<- format(round(norm.2[,6:23], 1), nsmall = 1)      # ensure that only have 1 decimal places
mat.data <- data.matrix(rounded.norm)          # convert data to matrix
rownames(mat.data) <- rnames                         # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normalized2.png", sep = ""),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "WSA paper 
          (normalized by total submissions)", # heat map title
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

#########################################################################################################
# Categorize according to allocate weighting to each categorization
# group according to each policy area. Multiply the number of respondents by..

# final score (for stakeholder group for policy area) = 
# (# respondents for increase/total number of respondents)*scalingfactor 1
# + (# respondents for maintain/total number of respondents)*scalingfactor 2
# + (# respondents for decrease/total number of respondents)*scalingfactor 3

#multiply each 'increase' segment by 101 in normalized (according to policy area)

increase = subset(normalized, normalized$Status == 'More, stronger regulation')
increase.factor = data.frame(increase[,1:19]*101, increase[,20:24])

maintain = subset(normalized, normalized$Status == 'Moderate regulation')
maintain.factor= data.frame(maintain[,1:19]*51, maintain[,20:24])

decrease = subset(normalized, normalized$Status == 'Less, weaker regulation')
decrease.factor= data.frame(decrease[,1:19]*1, decrease[,20:24])

factor <- rbind(increase.factor, decrease.factor, maintain.factor) # bind all three subsets back together

#### Graph heatmaps for the first normalization method (by number of responses within a specific policy focus)
# First, graph according to sub policy areas
total.factors.1<- data.frame(aggregate(factor[,1:19], by=list(Policy = factor$Sub_Policy), FUN = sum))

# find row names
rnames <- total.factors.1$Policy
rounded.factors<- format(round(total.factors.1[,2:19], 1), nsmall = 1)      # ensure that only have 2 decimal places
mat.data <- data.matrix(rounded.factors)                          # convert data to matrix
rownames(mat.data) <- rnames                                      # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normfactor1_Sub.png", sep = ""),    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "Response to policy areas - Sub policy areas (normalized by # respondents)", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 1,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv=TRUE)            # turn off column clustering

###### Try where the stakeholder groups are the rows
# find row names
rnames <- colnames(total.factors.1)
rounded.factors<- format(round(total.factors.1[,2:19], 1), nsmall = 1)      # ensure that only have 2 decimal places
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix
colnames(mat.data) <- total.factors.1$Policy                                     # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normfactor1_Sub_b.png", sep = ""),    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "Response to policy areas - Sub policy areas (normalized by # respondents)", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 1,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv=TRUE)            # turn off column clustering


##### Graph by the main policy areas (Seven in total)
# sum all three groupings together by policy area
total.factors.main1<- data.frame(aggregate(factor[,1:19], by=list(Policy = factor$Main_Policy), FUN = sum))

# find row names
rnames <- total.factors.main1$Policy
rounded.factors<- format(round(total.factors.main1[,2:19], 1), nsmall = 1)      # ensure that only have 2 decimal places
mat.data <- data.matrix(rounded.factors)            # convert data to matrix
rownames(mat.data) <- rnames                             # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normfactorsmain.png", sep = ""),    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "Response to policy areas - main policy areas", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,10),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 0.6,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv=TRUE)            # turn off column clustering

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

################
##### Influence maps - normalized according to the number of respondents (method 2)
# multiply each 'increase' segment by 101 in normalized

increase = subset(norm.2, norm.2$Status == 'More, stronger regulation')
increase.factor = data.frame(increase[,6:23]*101, increase[,1:5])

maintain = subset(norm.2, norm.2$Status == 'Moderate regulation')
maintain.factor= data.frame(maintain[,6:23]*51, maintain[,1:5])

decrease = subset(norm.2, norm.2$Status == 'Less, weaker regulation')
decrease.factor= data.frame(decrease[,6:23]*1, decrease[,1:5])

factor <- rbind(increase.factor, decrease.factor, maintain.factor) # bind all three subsets back together

### Graph heatmaps according to specific policy approach - second normalization approach
# First, graph according to sub policy areas
total.factors.2<- data.frame(aggregate(factor[,1:18], by=list(Policy = factor$Sub_Policy), FUN = sum))

# find row names
rnames <- total.factors.2$Policy
rounded.factors<- format(round(total.factors.2[,2:19], 1), nsmall = 1)      # ensure that only have 2 decimal places
mat.data <- data.matrix(rounded.factors)                          # convert data to matrix
rownames(mat.data) <- rnames                                      # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normfactor2_Sub.png", sep = ""),    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "Response to policy areas - Sub policy areas (normalized by # stakeholders)", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 1,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv=TRUE)            # turn off column clustering

###### Try where the stakeholder groups are the rows - norm method 2
# find row names
rnames <- colnames(total.factors.2)
rounded.factors<- format(round(total.factors.2[,2:19], 1), nsmall = 1) # ensure that only have 2 decimal places
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix
colnames(mat.data) <- total.factors.2$Policy                           # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_normfactor2_Sub_b.png", sep = ""),    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "Response to policy areas - Sub policy areas (normalized by # stakeholders)", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(10,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 1,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv=TRUE)            # turn off column clustering



############################################
##### box plots - to answer   question of whether stakeholder groups called for more regulation uniformly along policy areas
box.data.1 <- data.frame(total.factors.1[,2:19])
png(paste(directory, "WSA_norm1_boxplot.png", sep = ""),    # create PNG for the heat map        
    width = 7*300,        # 5 x 300 pixels
    height = 7*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

boxplot.matrix(as.matrix(box.data.1), use.cols = TRUE, 
               main=toupper("Distribution - Weighted Response Factor (normalization method 1"), font.main=10, 
               cex.main=1.2, xlab="Stakeholder", ylab="Response factor", font.lab=10, 
               col="darkgreen", las=3, margins =c(8,50)
               )

# second normalization method
box.data.2 <- data.frame(total.factors.2[,2:19])
png(paste(directory, "WSA_norm2_boxplot.png", sep = ""),    # create PNG for the heat map        
    width = 7*300,        # 5 x 300 pixels
    height = 7*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

boxplot.matrix(as.matrix(box.data.2), use.cols = TRUE, 
               main=toupper("Distribution - Weighted Response Factor (normalization method 2)", font.main=10, 
               cex.main=1.2, xlab="Stakeholder", ylab="Response factor", font.lab=10, 
               col="darkgreen", las=3, margins =c(8,50)))

# overlay where the act itself is in the boxplots?

############### influence maps
# Compare what was placed in the Act to pespectives of different groups
# Calculate according to:
# Influence = What the group wanted for that policy area - what was in the WSA



##### Decision prediction
# Try and use decision prediction
# Assumes that decision made is a combination of different inputs

############### 
# Heat maps of specific policy areas - normalized data
# Regulate and Protect Groundwater Use - norm method 1
groundwater.1 <- subset(normalized, normalized$Main_Policy == "Regulate and Protect Groundwater Use")

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(groundwater.1[,1:18])
rounded.factors<- format(round(groundwater.1[,1:18], 1), nsmall = 1) # ensure that only have 2 decimal places
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix
colnames(mat.data) <- paste(groundwater.1$Sub_Policy, groundwater.1$Status, sep="_") 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_groundwater_norm1.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "Response to groundwater - Norm method 1", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 1,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv=TRUE)            # turn off column clustering

###
# Regulate Groundwater - Normalized by the number of stakeholders within a group (norm method 2)
groundwater.2 <- subset(norm.2, norm.2$Main_Policy == "Regulate and Protect Groundwater Use")

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(groundwater.2[,6:23])
rounded.factors<- format(round(groundwater.2[,6:23], 1), nsmall = 1) # ensure that only have 2 decimal places
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix
colnames(mat.data) <- paste(groundwater.2$Sub_Policy, groundwater.2$Status, sep="_")                   # assign row names 

# creates a 5 x 5 inch image
png(paste(directory, "WSA_groundwater_norm2.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          cellnote = mat.data,  # same data set for cell labels
          main = "Response to groundwater - Norm method 2", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(20,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 1,          # decrease row font size to fit
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv=TRUE)            # turn off column clustering
