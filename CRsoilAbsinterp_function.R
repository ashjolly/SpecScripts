
abs.interpret <- function(filelist_EEMScor, directoryCorrectedEEMS, directoryCorrectedAbs){

  interp.abs <- function(abs){
    colnames(abs) <- seq(800,239, by = -3)
    # interpolate for all excitation to regular emission from 248-831 nm, 1 nm incrementes
    # if the excitation doesn't go in 2 nm increments, interpolate so that it does 
    if (as.numeric(colnames(abs)[1]) - as.numeric(colnames(abs)[2]) != 2){
      ex.seq <- seq(800, 240, by = -2)
      abs.inter <- as.data.frame(t(approx(x = as.numeric(colnames(abs)), y = abs[1,], xout = ex.seq)$y))
      colnames(abs.inter) <- ex.seq
    } 
    return(abs.inter)
  }

  # do for all of the abs file and return
  n = dim(filelist_EEMScor)[1]
  for (i in 1:n){
  
    samplename <- toString(filelist_EEMScor[i,1]) #column with the sample IDs
  
    ###########
    # load the Abs file
    setwd(directoryCorrectedAbs)
    abs.temp.corr <-as.data.frame(read.delim(as.character(filelist_EEMScor[i,3]), 
                                           header= TRUE, sep = ",", stringsAsFactors=FALSE))
    # apply interpolation function
    abs.interp <- interp.abs(abs = abs.temp.corr)
    # write interpolated file
    corrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_CorrInterp",".csv", sep = ""))
    write.table(abs.interp, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
  }
}

