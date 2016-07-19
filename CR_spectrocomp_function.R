#
# Function for compiling par file
# for use within CR analysis
#
########


#input .fp files
spec <- function(spectro.direct){
setwd(spectro.direct)
file_list <- list.files(pattern = ".fp$")
dataset <- ldply(file_list, read.delim, header=FALSE, skip=2, na.strings = "NaN",col.names=top)


############### 2a. Arrange merged fp files according to date

dataset = data.frame(dataset)

#Optional - export dataset to see format
#write.table(dataset, file = "C:/Users/Ashlee/Documents/compiledfp.csv",
#sep = ",", col.names = NA, qmethod = "double")

#get rid of duplicate entries
uniquedataset = unique(dataset)  

#sort according to date (2009-2013)
uniquedataset$DateTime <- as.POSIXct(strptime(uniquedataset$DateTime, format = "%Y.%m.%d  %H:%M:%S", tz="GMT"))
sort1.fp <- uniquedataset[order((uniquedataset$DateTime)),]    # sort by date/time

#check 
#write.table(sort1.dataset, file = "/Users/ashlee/Dropbox/par and fp compilation/fpcomp/fp_unique.csv",
#sep = ",", col.names = NA, qmethod = "double")

################# 1b. Read in PAR files into one master file

#column names for par files, RIVERCOL 
toppar = c('DateTime', 'Status',	'Turb.FTU',	'Turb.status',	'NO3.N',	'NO3.status',	'TOC',	'TOC.status',	'DOC',	'DOC.status','SAC254',	'SAC254.status',	'Level [m]NaN-NaN_2',	'[Level_0.0_1.0_0.0_0.0]',	'Temp [°C]NaN-NaN_2',	'[Temp_0.0_1.0_0.0_0.0]',	'analogINNaN-NaN_2',	'[analogIN_0.0_1.0_0.0_0.0]')

#input .par files
file_listpar <- list.files(pattern = ".par$")
datasetpar <- ldply(file_listpar, read.delim, header=FALSE, skip=2,na.strings = "NaN",col.names=toppar)

#Optional - export dataset to see format - ok
#write.table(dataset, file = "C:/Users/Ashlee/Documents/compiledfp.csv",
#sep = ",", col.names = NA, qmethod = "double")

############### 2b. Arrange merged par files according to date

datasetpar = data.frame(datasetpar)

#Optional - export dataset to see format
#write.table(dataset, file = "C:/Users/Ashlee/Documents/compiledfp.csv",
#sep = ",", col.names = NA, qmethod = "double")

#get rid of duplicate entries
uniquedatasetpar = unique(datasetpar)  

#sort according to date (2013-2014 files)
uniquedatasetpar$DateTime <- as.POSIXct(strptime(uniquedatasetpar$DateTime, format = "%Y.%m.%d  %H:%M:%S", tz="GMT"))
sort1.datasetpar <- uniquedatasetpar[order((uniquedatasetpar$DateTime)),]    

#check 
#write.table(sort1.dataset, file = "/Users/ashlee/Dropbox/par and fp compilation/fpcomp/fp_unique.csv",
#sep = ",", col.names = NA, qmethod = "double")

############# 
###3. Clean par file

############ 3a. Subset of status = OK
parOK = subset(sort1.datasetpar, sort1.datasetpar$Status == 'Ok')
parOK = data.frame(parOK)

############ 3b. Creating a subsets: Threshold conditions
DOCcleaned1 = subset(parOK, as.numeric((parOK$DOC)[parOK$DOC]) >= 0.950 )  #create a subset of data that contains only the data above 0.95.  mainly to get rid of the blanks/error
DOCcleaned2 = subset(DOCcleaned1, as.numeric((Turb.FTU)[Turb.FTU]) < 170.000) # subsets data where turbidity <=170.  Note that a turbidity above this likely indicative of when the spectro was buried, thus data collected during this time likely not good, but the threshold point chosen on the basis of the limit for 

############ 3c. DOC: Correct for the shimadzu (lab corrected value)
#  Streamwater correction
#If spectro serial number =  9210014

DOCcleaned2$DOCcorr = as.numeric(DOCcleaned2$DOC)*1.0271-0.1998  
correctedpar = DOCcleaned2

############# 3d. Calculating SUVA
correctedpar$SUVA = correctedpar$SAC254/correctedpar$DOCcorr

#####################
#### 4. Merging corrected par files with fp files
# 4a. merge files
# note that all - FALSE meanse that only those entries present in cleaned par file will be merged
mergeddata = merge(correctedpar, sort1.fp, by="DateTime", all = FALSE)

# 4b. rename DateTime column for openair function
mergeddata$DateTime <- as.POSIXct(strptime(mergeddata$DateTime, format = "%Y-%m-%d  %H:%M:%S", tz="GMT"))
colnames(mergeddata)[1] <- "date"

############## 5. Apply openair to create 30 min dataset
#install.packages("openair")
library(openair)
CR.ave <- timeAverage(mergeddata, avg.time = "30 min", start.date = "2013-10-20", fill = TRUE)
#note for above - may have to check start.date

#check - 
write.table(CR.ave, file = "/Users/ashlee/Documents/CR_Data/openairtest.csv",
            sep = ",", col.names = NA, qmethod = "double")



}