# Script for running corrections loops for EEMS
# This function will basically run different functions for running corrections
# Basically will take all uncorrected EEMS within a folder, and apply functions for corrections
# Avoids having the loop within the script - idea is to make corrections more adaptable between different projects and easier to apply
# 7aug2015 Started
# Ashlee Jollymore's PhD project
################

EEMcorrection = function(data.3, directoryall, directoryCorrectedEEMS, slitwidth1, slitwidth2, 
                         ex.wavelengths, em.wavelengths){
  
  # set working directory where all uncorrected EEMS, blank and absorbance files are
  setwd(directoryall)
  
  # Run corrections on all files contained within the data.3 vector
  # insert empty variables for populating with ex an em vectors
  n = nrow(data.3) #number of files you are going to correct in the file. Double check ths prior to proceeding
  
  ex_all = data.frame(matrix(vector(), 5000, n))
  em_all = data.frame(matrix(vector(), 5000, n))  
  
  #ex and em from blank - create dataframe
  ex_blank = data.frame(matrix(vector(), 5000, n))
  em_blank = data.frame(matrix(vector(), 5000, n)) 
  
  # 5000 is a placeholder for the number of rows. Note that this speeds the loop in comparison to leaving it as zero.
  
  ### set up loop to correct files and calculate indicies from all files in folder
  
  for (i in 1:n){
    
    #### load and trim files: functions to load, trim and correct EEMS data - blank, EEM and files
    #sample name
    samplename <- toString(data.3[i,1]) #column with the sample IDs
    
    #### load and trim EEMS
    #call function. Note that wd will change to sample WD in function
    setwd("/Users/user/SpecScripts") 
    source("EEMfileLoadTrim_function.R")
    
    EEM <- EEMtrim(graphheadings = data.3, samplewd = directoryall, loopnum = i)
    
    #ex and em wavelengths
    ex = as.numeric((sort(colnames(EEM), decreasing = TRUE)))
    em = as.numeric((sort(rownames(EEM), decreasing = FALSE)))
    
    #em_all is a variable that holds all of the em wavelengths. It's a way of seeing if em wavelengths are different between samples
    em_all = rbind(em_all, rownames(EEM))
    #ex_all is a variable that holds all of the ex wavelengths. It's a way of seeing if ex wavelengths are different between samples
    ex_all= rbind(ex_all, colnames(EEM))
    
    #### load and trim ABS
    # call function. Note that wd will change to sample WD in function
    setwd("/Users/user/SpecScripts") 
    source("EEMAbsLoadTrim_function.R")
    
    Abstrim <- ABStrim(graphheadings = data.3, samplewd = directoryall, loopnum = i)
    
    #### load and trim blank
    setwd("/Users/user/SpecScripts")
    source("EEMBlankLoadTrim_function.R")
    
    Blktrim <- BLANKtrim(graphheadings = data.3, samplewd = directoryall, loopnum = i)
    
    #ex and em_all are variables that output the complete ex and em. Make sure that all of your samples are the same
    ex_blank = rbind(ex_blank, colnames(Blktrim))
    em_blank = rbind(em_blank, rownames(Blktrim))
    
    #### identify dilution factor in master file
    # Dilution = column 5 in data.3
    dil = data.3[i,5]
    
    ################################## Corrections
    ########### IFE: Correct raw EEM for IFE if sample has not been corrected for this
    # Determine if needs IFE done if EEM samplename is SYM.dat
    # Otherwise, if sample name is PEM, IFE has been done in the software
    EEMsampletype <- strapplyc(as.character(data.3[i,2]), paste(samplename, "(.*).dat", sep = ""), simplify = TRUE)
    
    if (EEMsampletype == "SYM") {
      # call function
      setwd("/Users/user/SpecScripts") 
      source("EEMIFECorr_function.R")
      
      EEM.IFC <- innerfilter(eem = EEM, abs = Abstrim, em = em, ex = ex)
      # note that IFE should be between 0.4 and  0.98 according to McKnight 2001 (doi: 10.4319/lo.2001.46.1.0038)
    }
    
    # if sample type = PEM, do NOT do IFE correction - already done in file
    if (EEMsampletype == "PEM"){
      EEM.IFC <- EEM
    }
    
    ########### Normalize IFE EEM and blank file according to area under Raman peak
    # Need to do for all EEMs exported from Aqualog
    # call function
    setwd("/Users/user/SpecScripts") 
    source("Ramancorrect_v1.R")
    
    # wavelengths for Raman corrections
    # tell R where em = 375 nm, em = 430 nm; ex = 350 nm
    em375 <- as.numeric(grep(em.wavelengths[c("em.375"),], rownames(EEM)))
    em430 <- as.numeric(grep(em.wavelengths[c("em.430"),], rownames(EEM)))
    ex350 <- as.numeric(match(ex.wavelengths[c("ex.350"),], colnames(EEM)))
    
    # get the Raman correction file from the Raman function stored
    Raman.area <- Ramancor(blank = Blktrim) 
    
    # normalize the EEM file for the area underneath the raman curve calculated above 
    # Raman normalize the raw EEM
    EEM.ram = EEM.IFC/Raman.area 
    
    # Raman normalize the blank file (important for blank correction)
    blankram = Blktrim/Raman.area
    
    ##################################
    ########### Correct for Blank
    # Only if have not been done in software
    # corrected EEM = Sample - blank. Both SYM and PEM have been instrument corrected
    
    #if (EEMsampletype == "") {
    #  EEM.blk <- EEM.ram - blankram
    #}
    
    # if sample type = SYM and PEM, do NOT do blank correction
    #elseif {
    #  EEM.blk <- EEM.ram
    #}
    ##### note that SYM + PEM files from Aqualog have already been blank corrected
    
    ###########################
    ##### Apply dilution factor to EEM and to Abs file
    EEM.dil = EEM.ram*dil 
    Absdil = Abstrim*dil #second column = absorbance readings. Assumes Beer's Law applies (c~abs)
    
    ##################################
    ########### Correct for Raleigh Masking
    # Done via interpolation software avalaible for Matlab script - 'eemscat.m'
    # Access on 23July2015 from http://www.models.life.ku.dk/EEM_correction
    # Also see DOI http://dx.doi.org/10.1016/j.chemolab.2015.01.017 explaination of why interpolating
    # through the Raleigh scatter is the best option in terms of PARAFAC modelling
    
    # R script emulates the interpolation script found here (interpolation according to the matlab file above)
    # Don't know why I can;t get that gawd damned script to work. Oh well ;)
    
    # this portion of the script runs the eemscat.m function to get the interpolated spectra. Last correction before saving
    
    #if (EEMsampletype == "SYM"){
      
    #  # call function
     # setwd("/Users/ashlee/SpecScripts") 
    #  source("EEMRaleigh_function.R")
    #  # note that this will gap fill the second order Raleigh scatter with na.spline function in zoo
    #  EEM.rm <- raleigh(eem = EEM.dil, slitwidth1, slitwidth2)
   # }
    
    # if Raleigh has already been done in Aqualog software (inserted 0's, not the best option)
   # if(EEMsampletype == "PEM"){
     # EEM.rm <- EEM.dil
   # }
    
    EEM.rm <- EEM.dil
    ##### Apply correction factor for Fe concentration
    ##### TO DO!!!!!
    
    # change the corrected EEM so that excitation goes from 200-800, rather than 800-200
    EEMcorr <- EEM.rm[,sort(names(EEM.rm), decreasing = FALSE)]
    
    ###########
    ##### Save the corrected EEM
    corrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_Corrected",".csv", sep = ""))
    write.table(EEMcorr, file = corrpath, row.names = TRUE,col.names = TRUE, sep = ",")
    
    ###########
    ##### Save the corrected absorbance file (trimmed and dilution factor accounted for). In same folder as EEMS
    abscorrpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project,"_AbsCorrected",".csv", sep = ""))
    write.table(Absdil, file = abscorrpath, row.names = TRUE,col.names = TRUE, sep = ",")
    
    ##########
    # last - plot corrected eems as a contour plot
    #variables to change
    zmax = max(EEMcorr,na.rm=TRUE) # put the max intensity of that you want to graph
    #EEMmax[i] <- zmax #to show the maximum fluorescence for all files
    xlimit <- range(300, 700, finite=TRUE)
    ylimit <- range(240, 450, finite = TRUE)
    
    numcont = 100 # number of contour levels you want: Change if you want
    
    ##### contour plotting function
    # call contour plot function
    setwd("/Users/user/SpecScripts") 
    source("EEM_contour_v1.R")
    
    #Plot contours and save in correction file
    plotpath <- file.path(directoryCorrectedEEMS, paste(samplename,"_", project, "_Contour",".jpeg", sep = ""))
    
    g <- length(EEMcorr)
    EEMplot <- EEMcorr # not cutting out the last two columns
    EEMplot[EEMplot<0]=0 # have line so that if fluorescence is less than 0, becomes 0.
    
    explot = as.numeric(colnames(as.matrix(EEMplot)))
    emplot = as.numeric(row.names(EEMplot))
    
    jpeg(file=plotpath)
    contour.plots(eems = as.matrix(EEMplot), Title = paste(samplename, project, sep = "_"), ex = explot, em = emplot, 
                  zmax = zmax, zmin = 0, numcont = numcont)  
    dev.off()

    # note that the above is meant to be a crude graphing - better graphing done in matlab once
    # you figure out the max emission for your dataset (normalize all of the plots to this)
  }
  
  #### End of corrections loop!
  wavelengths = cbind(ex_all, em_all, ex_blank, em_blank)
  return(wavelengths)
  
}