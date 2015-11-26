###
# function for trimming files for CM modelling
# Script is meant to take corrected EEM files and ensure that they have the same EX and EM ranges 
# so that file can be passed into CM model in matlab

# 8July2015
############

CMtrim <- function(filedirectory, filelist, project, exmin, directoryCM, ex){
  
  # create graph heading variable
  graphheadings = data.frame((0))
  
  # set working directory with filder that has corrected files
  setwd(filedirectory)
  
  n = length(filelist)
  for (i in 1:n){
    temp.EEMS <- read.delim(filelist[i], header= TRUE, sep = ",")
    
    #trim so that exitation and emission goes from the same
    ex.temp <- colnames(temp.EEMS)
    em = row.names(temp.EEMS)
    
    if(ex.temp[1] != exmin) {
      # if first value in ex.temp is not 240, trim 
      ex.length <- length(ex.temp)
      # find column where the exitation wavelength is 240 to cut from
      x240 = as.numeric(match(exmin,names(temp.EEMS)))
      temp.EEMS <- temp.EEMS[,c(x240:ex.length)]
    } 
    
    # Create graph headings variable. Also insert "_i" to use in CM modelling
    samplename <- strapplyc(filelist_EEMScor[i], paste("(.*)_", project, "_Corrected", sep = ""), simplify = TRUE)
    graphheadings[i,] <-paste(samplename, project,"CorrCM_",i, sep = "")
    
    # Resave trimmed EEM without the row and column names
    corrpath <- file.path(directoryCM, paste(samplename, project,"CorrCM_",i,".csv", sep = ""))
    write.table(temp.EEMS, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")
    
  }
  
  # save ex and em in two separate files, to make it easier to read into CM PARAFAC files
  corrpath <- file.path("/Users/user/Documents/MATLAB/ExEmfiles", paste(project,"em",".csv", sep = ""))
  write.table(em, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")
  
  corrpath <- file.path("/Users/user/Documents/MATLAB/ExEmfiles", paste(project,"ex",".csv", sep = ""))
  write.table(ex.PARAFAC, file = corrpath, row.names = FALSE,col.names = FALSE, sep = ",")
  
  # write graph headings file
  corrpath <- file.path("/Users/user/Documents/MATLAB/CM_graphheadings", paste("GraphHeadings_", project,".txt", sep = ""))
  write.table(graphheadings, file = corrpath, row.names= FALSE, col.names = FALSE, sep= ",")
  
  return(graphheadings)
}