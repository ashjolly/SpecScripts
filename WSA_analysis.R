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
directory <- "/Users/user/Dropbox/PhD Work/WaterAct Paper/"

# libraries
library('gplots')
library('plyr')
library("d3heatmap")
library('RColorBrewer')
library('scales')
library('ggplot2')

## load data
data.original <- as.data.frame(read.csv(paste(directory, 'codestats_heatmap.csv', sep = ''), sep = ",", header = TRUE))

################# Data Preprocessing
# take out the 'all' categories from the data
data <- subset(data.original, data.original$Status != "All")

# Make two new columns that contain the Main, Sub and status of the policy areas
data$Policy_response <- paste(data$Main_Policy, data$Sub_Policy, data$Status, sep="_") # merge first three columns into third column
data$Policy_only <- paste(data$Main_Policy, data$Sub_Policy, sep="_")

# combine 'health' with community groups (only one respondent in 'health')
data$Community.Groups <- data$Community.Groups + data$Health
data$Health <- NULL

# combine the two individual groups
data$Individuals <- data$Individual.form.submissions.total..10..+data$Individual.non.form.submissions.total..10..

# Get rid of mining?
##################### Data normalization
## Normalize data by the total number of responses within a policy area
# First, by the main policy area
# group by policy group + count total number of responses
# According to the seven main policy areas
total.responses <- aggregate(data[,c(7:21, 26)], by=list(Policy = data$Main_Policy), FUN = sum)

# According to the sub policy areas within the larger policy area
subpolicy.responses <- aggregate(data[,c(7:21, 26)], by=list(Policy = data$Policy_only), FUN = sum)

main.policyareas <- unique(data$Main_Policy)    # find the total number of unique main policy areas
sub.policyareas <- unique(data$Policy_only)     # find the total number of unique sub policy areas

stakeholders <- unique(colnames(data[,c(7:21, 26)]))   # find the total number of stakeholder groups

########
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
  temp.normalize <- as.data.frame(sweep(as.matrix(temp.subset[,c(7:21, 26)]), 2, as.matrix(temp.sum[,2:17]), "/"))
  
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
# normalized[is.na(normalized)] = 0

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
rounded.norm<- format(round(normalized[,1:16], 1), nsmall = 1)      # ensure that only have 1 decimal places
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
d3heatmap(mat.data, scale = "column", 
          dendrogram = "row",
          col=my_palette,
          margins =c(8,50), 
          cexRow=0.6, 
          cexCol = 0.35,
          key = TRUE
)

#########################################################################################################
# Categorize according to allocate weighting to each categorization
# group according to each policy area. Multiply the number of respondents by..

# final score (for stakeholder group for policy area) = 
# (# respondents for increase/total number of respondents)*scalingfactor 1
# + (# respondents for maintain/total number of respondents)*scalingfactor 2
# + (# respondents for decrease/total number of respondents)*scalingfactor 3

#multiply each 'increase' segment by 101 in normalized (according to policy area)

increase = subset(normalized, normalized$Status == 'More, stronger regulation')
increase.factor = data.frame(increase[,1:16]*101, increase[,17:21])

maintain = subset(normalized, normalized$Status == 'Moderate regulation')
maintain.factor= data.frame(maintain[,1:16]*51, maintain[,17:21])

decrease = subset(normalized, normalized$Status == 'Less, weaker regulation')
decrease.factor= data.frame(decrease[,1:16]*1, decrease[,17:21])

factor <- rbind(increase.factor, decrease.factor, maintain.factor) # bind all three subsets back together
remove(increase, increase.factor, maintain, maintain.factor, decrease, decrease.factor)

#### Graph heatmaps for the first normalization method (by number of responses within a specific policy focus)

# Graph according to sub policy areas
total.factors.1<- data.frame(aggregate(factor[,1:16], by=list(Policy = factor$Policy_only), FUN = sum))
# create column with sub policy area
total.factors.1$Sub_policy <- gsub("^.*?_","", total.factors.1$Policy)
total.factors.1$Main_Policy <- sub("(.*?)_.*","\\1", total.factors.1$Policy)

#sort by Main Policy area
total.factors.1 <- total.factors.1[with(total.factors.1, order(total.factors.1$Main_Policy)), ]

# take away column with main policy - leave only sub policy
total.factors <- total.factors.1[,c(-1,-19)]

######  Heat map
#Where the stakeholder groups are the rows (Thus, grouping by stakeholdr groups)
# find row names
rnames <- colnames(total.factors)
cnames <- total.factors$Sub_policy 
rounded.factors<- format(round(total.factors[,1:16], 1), nsmall = 1)      # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames                                # assign row names
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix
                                        
# change colour pallet
palette <- colorRampPalette(c("#E5F5F9", "#005824"))(n = 299)

# creates a 5 x 5 inch image
png(paste(directory, "WSA_response_subpolicy.png", sep = ""),    # create PNG for the heat map        
    width = 11*300,        # 5 x 300 pixels
    height = 11*300,
    res = 600,            # 300 pixels per inch
    pointsize = 5)        # smaller font size

heatmap.2(mat.data,
          # Change the data within the heat map boxes
          cellnote = mat.data,  # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          # labels
          main = "Stakeholder Response to Specific Policy Areas (Normalized)", # heat map title
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="histogram",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.5, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=TRUE)            # turn off column clustering
dev.off()

##### Graph by the main policy areas (Seven in total). Group by stakeholder groups
# sum all three groupings together by policy area
total.factors.main1<- data.frame(aggregate(factor[,1:16], by=list(Policy = factor$Main_Policy), FUN = sum))

# find row names
rnames <- colnames(total.factors.main1)
rounded.factors<- format(round(total.factors.main1[,2:17], 1), nsmall = 1)      # ensure that only have 2 decimal places
mat.data <- t(data.matrix(rounded.factors))            # convert data to matrix
colnames(mat.data) <- total.factors.main1$Policy                             # assign row names

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
dev.off()

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


############################################
##### box plots - to answer   question of whether stakeholder groups called for more regulation uniformly along policy areas
box.data.1 <- data.frame(total.factors.1[,2:17])
png(paste(directory, "WSA_boxplot.png", sep = ""),    # create PNG for the heat map        
    width = 7*300,        # 5 x 300 pixels
    height = 7*300,
    res = 300,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

op <- par(mar = c(14,4,4,4) + 2) # increase margins, especialy x axis margins

plot <-boxplot.matrix(as.matrix(box.data.1), use.cols = TRUE, 
               main=toupper("Distribution - Weighted Response Factor"), font.main=12, 
               cex.main=1.2, 
               xlab="", 
               ylab = "Percent Response", cex.axis=1.5,cex.lab = 2,
               col="#005824", las=3, margins =c(8,50), outline=FALSE
               )
par(op)
dev.off()
# add x labels

# overlay where the act itself is in the boxplots?

############### influence maps
# Compare what was placed in the Act to pespectives of different groups
# Calculate according to:

# Influence = What the group wanted for that policy area - what was in the WSA

 # try #1 - influence = WSA * response factor (calculated previously)
# Step 1: Calculate the WSA response factor

WSA.increase = subset(data, data$Status == 'More, stronger regulation')
WSA.increasefactor = data.frame(WSA.increase[,23]*101, WSA.increase[,c(1:3, 24,25)])
colnames(WSA.increasefactor)[1] <- "Response"

WSA.maintain = subset(data, data$Status == 'Moderate regulation')
WSA.maintainfactor= data.frame(WSA.maintain[,23]*51, WSA.maintain[,c(1:3, 24,25)])
colnames(WSA.maintainfactor)[1] <- "Response"

WSA.decrease = subset(data, data$Status == 'Less, weaker regulation')
WSA.decreasefactor= data.frame(WSA.decrease[,23]*1, WSA.decrease[,c(1:3, 24,25)])
colnames(WSA.decreasefactor)[1] <- "Response"

WSA.factors <- rbind(WSA.increasefactor, WSA.maintainfactor, WSA.decreasefactor) # bind all three subsets back together
remove(WSA.maintain, WSA.increase, WSA.increasefactor, WSA.maintainfactor, WSA.decrease, WSA.decreasefactor)

# Aggregate WSA response according to sub policy area
WSA.totalfactors <- data.frame(aggregate(WSA.factors[,1], by=list(Policy = WSA.factors$Policy_only), FUN = sum))

# create column with sub policy area
WSA.totalfactors$Sub_policy <- gsub("^.*?_","", WSA.totalfactors$Policy)
WSA.totalfactors$Main_Policy <- sub("(.*?)_.*","\\1", WSA.totalfactors$Policy)

#sort by Main Policy area
WSA.totalfactors <- WSA.totalfactors[with(WSA.totalfactors, order(WSA.totalfactors$Main_Policy)), ]
# insert column name for WSA response
colnames(WSA.totalfactors)[2] <- "WSA_Response"

## calculate Influence, where Influence = WSAresponse - Stakeholder response
# Response factor for stakeholder in total.factors.1
response.all <- merge(total.factors.1, WSA.totalfactors, by = "Policy", all = TRUE)

# calculate influence for each stakeholder group (columns in data, columnes 2-17)
n = dim(response.all)[1]
remove(Influence)

for (i in 1:n){
    temp.influence <- data.frame(t(apply(response.all[i,2:17], 2, function(x) response.all[i,20] - x)))
    
    # if the merged dataset  exists, append to it by row
    if (exists("Influence")){
      Influence <- rbind(Influence, temp.influence)
      #rm(temp.cut)
    }
    
    # if the merged dataset doesn't exist, create it
    if (!exists("Influence")){
      Influence <- temp.influence
    }
}
remove(temp.influence)

Influence <- cbind(response.all$Sub_policy.x, Influence)
             
### Do heat map of the influence factors
# find row names
rnames <- colnames(Influence)
cnames <- Influence[,1]
rounded.factors<- format(round(Influence[,2:17], 1), nsmall = 0)      # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames      # assign column names to data matrix
mat.data <- t(data.matrix(rounded.factors))            # convert data to matrix

# creates a 5 x 5 inch image
png(paste(directory, "WSA_InfluenceMap_columngrouping.png", sep = ""),    # create PNG for the heat map        
    width = 10*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,    # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",        # change font color of cell labels to black
          
          # labels
          main = "Influence Map - WSA Sub-Policy Areas", # heat map title
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="histogram",  # turns on histogram plot inside color legend
          trace="none",         # turns off trace lines across both rows and columns
          
          # appearance
          margins =c(15,12),     # widens margins around plot
          col= my_palette,       # use on color palette defined earlier 
          cexCol=1, 
          cexRow = 1,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=TRUE)            # turn off column clustering


##### Decision prediction
# Try and use decision prediction
# Assumes that decision made is a combination of different inputs

############################## ############### ############### ############### 
# Heat maps of specific policy areas - normalized data
# Regulate and Protect Groundwater Use 
groundwater.1 <- subset(normalized, normalized$Sub_Policy == "Regulate Groundwater Extraction and Use")

# express as percent
groundwater = cbind(groundwater.1[,1:16] *100, groundwater.1$Sub_Policy, groundwater.1$Status)
colnames(groundwater)[18] <- "Status"
colnames(groundwater)[17] <- "Sub_Policy"

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(groundwater[,1:16])
cnames <- groundwater$Status
rounded.factors<- format(round(groundwater[,1:16], 1), nsmall = 1) # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix

# creates a 5 x 5 inch image
png(paste(directory, "WSA_groundwater.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,  # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Percent of Responses to Regulate Groundwater Extraction and Use", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="density",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.5, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=FALSE)            # turn off column clustering
dev.off()

#################################
 #Second Focus Area - FITFIR ("General' main policy area)
general <- subset(normalized, normalized$Sub_Policy == "Allocation system")

# express as percent
general = cbind(general[,1:16] *100, general$Sub_Policy, general$Status)
colnames(general)[18] <- "Status"
colnames(general)[17] <- "Sub_Policy"

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(general[,1:16])
cnames <- general$Status 
rounded.factors<- format(round(general[,1:16], 1), nsmall = 1) # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix

# creates a 5 x 5 inch image
png(paste(directory, "WSA_allocation.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,  # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Percent of Responses to FITFIR ('Allocation System' Sub Policy Area)", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="density",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.5, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=FALSE)            # turn off column clustering
dev.off()

#################
# grouping #3 - Environmental Flow Needs
EFN <- subset(normalized, normalized$Sub_Policy == "Environmental Flow Needs")

# express as percent
EFN = cbind(EFN[,1:16] *100, EFN$Sub_Policy, EFN$Status)
colnames(EFN)[18] <- "Status"
colnames(EFN)[17] <- "Sub_Policy"

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(EFN[,1:16])
cnames <- EFN$Status
rounded.factors<- format(round(EFN[,1:16], 1), nsmall = 1) # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix

# creates a 5 x 5 inch image
png(paste(directory, "WSA_EFN.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,  # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Percent of Responses to Environmental Flow Needs", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="density",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.4, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=FALSE)            # turn off column clustering
dev.off()

#################
# Grouping #4 - Enabling Shared Governance
ERGA <- subset(normalized, normalized$Sub_Policy == "Enable a Range of Governance Approaches")

# express as percent
ERGA = cbind(ERGA[,1:16] *100, ERGA$Sub_Policy, ERGA$Status)
colnames(ERGA)[18] <- "Status"
colnames(ERGA)[17] <- "Sub_Policy"

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(ERGA[,1:16])
cnames <- ERGA$Status
rounded.factors<- format(round(ERGA[,1:16], 1), nsmall = 1) # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix

# creates a 5 x 5 inch image
png(paste(directory, "WSA_ERGA.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,  # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Percent of Responses to Enabling a Range of Governance Approaches", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="density",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.5, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=FALSE)            # turn off column clustering
dev.off()

####################################################################
# percent of stakeholders that responsed to different policy areas

# Bar graph of the number of respondents by stakeholder group
sum.total <- subset(data, data$Status == "Submitter Category")

# creates a 5 x 5 inch image
png(paste(directory, "WSA_number of submissions.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

op <- par(mar = c(14,4,4,2) + 1) # increase magins, especialy x axis margins

plot <- barplot(as.matrix(sum.total[,c(7:21, 26)]), main = "Number of Submissions Per Stakeholder Group", cex.main = 2, 
                cex.axis = 2,
         col = "#005824", axes = FALSE, axisnames = FALSE)
par(op)        
text(plot, par("usr")[3]-.5, labels = (colnames(sum.total)[c(7:21, 26)]),srt = 45,
     adj = c(1.3,1.3), xpd = TRUE, cex=1.3)
axis(2)

# add to top of graph the total number for each stakeholder group
text(x = plot, y = (sum.total[,c(7:21, 26)])+12, labels=(sum.total[,c(7:21, 26)]), xpd=TRUE)

title(main = NULL, sub = NULL, 
      xlab = NULL, ylab = "Number of Responses per Stakeholder Group", cex.lab = 1.6)
dev.off()

#########################################################################################
### Calculate the percent of respondents per sub policy area
# Where percent = number of respondents in sub policy area/ total number of respondents in stakeholder grou
# combine 'health' with community groups (only one respondent in 'health')
data.original$Community.Groups <- data.original$Community.Groups + data$Health
data.original$Health <- NULL

# combine the two individual groups
data.original$Individuals <- data.original$Individual.form.submissions.total..10..+data.original$Individual.non.form.submissions.total..10..
data.original$Total <- NULL

num.submitters <- subset(data.original, data.original$Status == 'Submitter Category')

# Find the 'all' category for each subgroup in the data.original dataframe

status.all <- subset(data.original, data.original$Status == "All")

# for sub policy, remove the instances where the 'Sub_policy' heading == ALL
sub.all <- subset(status.all, status.all$Sub_Policy != "All")

# calculate the percent that respondent to each sub policy by dividing by the total number of submissions
sub.all <- data.frame(rbind(sub.all, num.submitters))

# loop over all rows to divide by number of submitters
remove(per.subpolicy)

for (i in c(7:21,23)) {
  temp.subpolicy <- data.frame(t(apply(data.frame(sub.all[,i]), 1, function(x) x/sub.all[27,i]*100)))
  
  # if the merged dataset  exists, append to it by row
  if (exists("per.subpolicy")){
    per.subpolicy <- rbind(per.subpolicy, temp.subpolicy)
    #rm(temp.cut)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("per.subpolicy")){
    per.subpolicy <- temp.subpolicy
  }
}
per.subpolicy <- t(per.subpolicy)

# Assign row and column names

row.names(per.subpolicy) <- sub.all$Sub_Policy # Note that row names are sub policy area
colnames(per.subpolicy) <- colnames(sub.all)[c(7:21,23)]

# Remove the row with the number of submissions per stakeholder group

per.subpolicy <- per.subpolicy[c(-27,-28),]

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(per.subpolicy)
cnames <- row.names(per.subpolicy)
rounded.factors<- round(per.subpolicy, digits = 0) # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix

# creates a 5 x 5 inch image
png(paste(directory, "WSA_PercentRespondents_colgroupings.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,  # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Percent of Respondents Per Policy Focus Across Stakeholder Groups", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="density",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.5, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=TRUE)            # turn off column clustering
dev.off()

## Do the same for main policy areas

# for sub policy, remove the instances where the 'Sub_policy' heading == ALL
main.all <- subset(status.all, status.all$Sub_Policy == "All")

# calculate the percent that respondent to each sub policy by dividing by the total number of submissions
main.all <- data.frame(rbind(main.all, num.submitters))

# loop over all rows to divide by number of submitters
remove(per.mainpolicy)

for (i in c(7:21,23)) {
  temp.mainpolicy <- data.frame(t(apply(data.frame(main.all[,i]), 1, function(x) x/main.all[8,i]*100)))
  
  # if the merged dataset  exists, append to it by row
  if (exists("per.mainpolicy")){
    per.mainpolicy <- rbind(per.mainpolicy, temp.mainpolicy)
    #rm(temp.cut)
  }
  
  # if the merged dataset doesn't exist, create it
  if (!exists("per.mainpolicy")){
    per.mainpolicy <- temp.mainpolicy
  }
  remove(temp.mainpolicy)
}

per.mainpolicy <- t(per.mainpolicy)

# Assign row and column names

row.names(per.mainpolicy) <- main.all$Main_Policy # Note that row names are sub policy area
colnames(per.mainpolicy) <- colnames(main.all)[c(7:21,23)]

# Remove the row with the number of submissions per stakeholder group

per.mainpolicy <- per.mainpolicy[c(-8),]

# do heat map
# find row names - where the stakeholder groups are the rows
rnames <- colnames(per.mainpolicy)
cnames <- row.names(per.mainpolicy)
rounded.factors<- round(per.mainpolicy, digits = 0) # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix

# creates a 5 x 5 inch image
png(paste(directory, "WSA_PercentRespondents_main.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,  # same data set for cell labels
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Percent of Respondents Per Main Policy Focus Across Stakeholder Groups", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("none"),     # only draw a row dendrogram
          density.info="density",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.2, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = TRUE,
          Colv=TRUE)            # turn off column clustering
dev.off()


################# What is the response to consultation?
# Show a bar graph split into the percent of stakeholders who expressed negative, positive and suggestions for change to consultation process
# read in csv with the data

consult.response <- as.data.frame(read.csv(paste(directory, 'ResponsetoConsultation.csv', sep = ''), sep = ",", header = TRUE))
consult.response <- consult.response[,c(1:20)] # trim

# combine 'health' with community groups (only one respondent in 'health')
consult.response$Community.Groups <- consult.response$Community.Groups + consult.response$Health
consult.response$Health <- NULL

# combine the two individual groups
consult.response$Individuals <- consult.response$Individual.form.submissions.total..10..+consult.response$Individual.non.form.submissions.total..10..

# get the percent who responded
consult.response$Total <- NULL
percent.response <- percent(as.matrix(consult.response[1,4:19]/consult.response[9,4:19]))
# round to nearest percent
percent.response.2 <- round(((consult.response[1,4:19]/consult.response[9,4:19])*100), digits = 0)

# create percentages for positive, neutral, negative
neg <- consult.response[2,2:19]/consult.response[1,2:19]*100
pos <- consult.response[3,2:19]/consult.response[1,2:19]*100
suggestions <- consult.response[4,2:19]/consult.response[1,2:19]*100

#neutral <- 100-pos-neg

percent.happy <- rbind(pos, neutral, suggestions)
row.names(percent.happy) <- c( "Positive Response", "Negative Response", "Offered Suggestion to Process")

############# Stacked Bar Plot with Colors and Legend - splitting into pos, neutral, negative
# creates a 5 x 5 inch image
png(paste(directory, "WSA_ResponsetoConsult_bargraph.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

op <- par(mar = c(10,2,4,12) + 2) # increase margins, especialy x axis margins

plot <- barplot(as.matrix(percent.happy[,3:18]), main = "Response to Consultation Process",
                col = (colorRampPalette(c("dark blue", "light blue"))(n = 3)), 
                #bty='L',
                axes = FALSE, axisnames = FALSE)

par(op)        # resize area
# add legend in top left outside of plot area
par(xpd=TRUE)
legend(19.5,100, legend = c("Negative", "Offered Suggestion", "Positive"), 
       fill = (colorRampPalette(c("light blue", "dark blue"))(n = 3)), title="Response")

# add x labels
text(plot, par("usr")[3]-.5, labels = (colnames(percent.happy)[c(3:18)]),srt = 45,
     adj = c(1.3,1.3), xpd = TRUE, cex=1)
axis(2)

# Add y axis label
title(main = NULL, sub = NULL, 
      xlab = NULL, ylab = "Percent of Respondents in Stakeholder Group", cex = 1)

# Add in numbers to top of bar graph showing percent in stakeholder group that responded
text(x = plot, y = 100+2, labels=(percent.response), xpd=TRUE)
dev.off()

############ express as heat map
# bind the percent that expressed an opinion to the percent happy
percent.happy.2 <- rbind(percent.response.2, percent.happy[,3:18])
row.names(percent.happy.2)[1] <- 'Percent That Commented on Process'

rnames <- colnames(percent.happy.2)
cnames <- row.names(percent.happy.2)
rounded.factors<- round(percent.happy.2, digits = 0) # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix
colnames(rounded.factors) <- cnames

# creates a 5 x 5 inch image
png(paste(directory, "WSA_Comments_heatmap_per.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size#

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = mat.data,  # number of respondents
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Response to Consultation Process", # heat map title
          
  #       dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("row"),     # only draw a row dendrogram
          density.info="density",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
        margins =c(20,20),     # widens margins around plot
         col= palette,       # use on color palette defined earlier 
         cexCol=1.5, 
         cexRow = 1.5,          # decrease row font size to fit
         srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
        na.color = 'white',   # colour of NA blocks
         keysize = 1,         # size of the colour key
         Rowv = TRUE,
        Colv=FALSE)            # turn off column clustering

dev.off()

#######
#### Percent that asked for changes to be made
percent.asked <- (as.matrix(consult.response[4,4:19]/consult.response[9,4:19]*100)) #show on top of bar graph

bettercom <- (as.matrix(consult.response[5,4:19]/consult.response[4,4:19]*100))
extend <- (as.matrix(consult.response[6,4:19]/consult.response[4,4:19]*100)) 
increase <- (as.matrix(consult.response[7,4:19]/consult.response[4,4:19]*100)) 
FN <- (as.matrix(consult.response[8,4:19]/consult.response[4,4:19]*100)) 

comments <- rbind(percent.asked, bettercom, extend, increase, FN)
colnames(comments) <- colnames(consult.response)[4:19]
comments <- round(comments, digits = 0)
row.names(comments) <- c("Expressed Suggestion", "Better communication regarding consultation",
                      "Extend comment period", "Increase opportunities for input", "Meaningful FN participation required")

# Express as a heat map...
# do heat map
# find row names - where the stakeholder groups are the rows
comments.n <- t(as.matrix(consult.response[4:8,4:19])) # raw number of respondents

rnames <- colnames(comments)
cnames <- row.names(comments)
rounded.factors<- comments # ensure that only have 2 decimal places
row.names(rounded.factors) <- cnames
mat.data <- t(data.matrix(rounded.factors))                          # convert data to matrix
colnames(comments.n) <- cnames

# creates a 5 x 5 inch image
png(paste(directory, "WSA_Comments_heatmap_per.png", sep = ""),         # create PNG for the heat map        
    width = 20*300,        # 5 x 300 pixels
    height = 10*300,
    res = 600,            # 300 pixels per inch
    pointsize = 6)        # smaller font size

heatmap.2(mat.data,
          
          # Change the data within the heat map boxes
          cellnote = comments.n,  # number of respondents
          notecex = 0.8,          # Change the font size of the data labels
          notecol="black",      # change font color of cell labels to black
          
          main = "Suggestions for Consultation Improvement", # heat map title
          
          # dendorgram and groupings
          #breaks=col_breaks,    # enable color transition at specified limits
          dendrogram=c("none"),     # only draw a row dendrogram
          density.info="none",  # turns on density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          
          # appearance
          margins =c(20,20),     # widens margins around plot
          col= palette,       # use on color palette defined earlier 
          cexCol=1.2, 
          cexRow = 1.5,          # decrease row font size to fit
          srtCol=45,           # rotate the x labels at 45 deg so they fit
          #axisnames = FALSE,
          
          na.color = 'white',   # colour of NA blocks
          keysize = 1,         # size of the colour key
          Rowv = FALSE,
          Colv=FALSE)            # turn off column clustering
dev.off()
