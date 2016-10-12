# 
# CR Event Analysis
# Calculate stats for events for CR - max, min, mean, sd
# - for grab samples
#######

event.grab <- function(data){
  #1. identify events
  start <- as.Date(c("2009-01-06 17:30", "2009-02-25 8:00", "2009-03-01 19:30", "2009-05-04 18:00", "2009-07-06 13:30",
                     "2009-08-01 0:00", "2009-08-04 8:00", "2009-08-06 6:30", "2010-01-10 21:00", "2010-01-14 15:00",
                     "2010-01-17 21:00", "2010-01-19 12:30", "2010-01-24 22:30", "2010-02-11 23:30", "2010-02-26 14:30",
                     "2010-03-02 23:30", "2010-03-11 6:30", "2010-03-29 7:30", "2010-04-01 0:00", "2010-04-27 7:30",
                     "2010-04-28 3:00", "2010-05-01 0:00", "2010-05-03 2:30", "2010-05-25 23:00", "2010-05-27 0:30",
                     "2010-05-31 9:00", "2010-06-15 2:30", "2010-08-07 13:00", "2010-08-27 20:30", "2010-08-31 14:00", 
                     "2010-09-21 1:00", "2010-09-26 2:00", "2010-10-24 5:00", "2010-11-07 2:00", "2010-11-09 17:00", "2010-11-14 2:00",
                     "2010-11-17 12:00","2010-11-30 3:00","2010-12-07 11:30","2011-01-07 5:00","2011-01-14 22:00","2011-02-03 17:00",
                     "2011-02-12 2:30","2011-02-20 19:30","2011-03-02 19:00","2011-03-09 12:30","2011-04-01 0:00","2011-05-13 3:00",
                     "2011-05-15 9:30","2011-05-19 17:00","2011-05-22 20:00","2011-05-27 20:30","2011-06-01 5:00","2011-06-15 0:00",
                     "2011-06-24 21:30","2011-07-01 0:00","2011-07-12 10:30","2011-07-13 22:00","2011-08-01 0:00","2011-09-22 1:30",
                     "2011-10-03 10:30","2011-11-22 3:00","2011-12-28 4:30","2012-01-03 19:30","2012-01-20 22:00","2012-01-22 19:30",
                     "2012-01-28 22:30","2012-02-08 12:00","2012-02-16 4:00","2012-02-17 5:00","2012-03-03 16:00","2012-03-09 7:00",
                     "2012-04-16 5:00","2012-04-25 10:30","2012-05-01 0:00","2012-06-07 6:30","2012-06-12 16:00","2012-06-19 3:00","2012-06-22 19:00",
                     "2012-07-16 2:30","2012-08-01 0:00","2012-08-13 4:00","2012-08-14 1:30","2012-08-19 13:00","2012-08-23 2:00","2012-08-24 2:30",
                     "2012-08-26 21:00","2012-08-29 1:00","2012-10-14 7:00","2012-10-18 20:00","2012-10-21 20:00","2012-10-23 6:00","2012-11-04 15:30",
                     "2012-11-06 6:30","2012-11-16 21:00","2012-11-28 15:30","2012-12-23 16:00","2012-12-25 6:00","2012-12-28 11:30","2013-01-05 4:30",
                     "2013-01-24 19:00","2013-01-25 16:00","2013-01-29 17:30","2013-02-24 19:00","2013-02-28 9:00","2013-03-19 23:30","2013-04-04 20:00",
                     "2013-04-10 4:00","2013-05-23 7:00","2013-05-26 15:00","2013-06-11 7:00","2013-06-18 2:00","2013-08-11 17:30","2013-08-27 10:00",
                     "2013-09-20 11:30","2013-09-22 10:30","2013-11-02 1:00","2012-11-06 19:30","2013-11-12 1:00","2013-11-15 18:00","2013-12-23 14:00",
                     "2013-12-27 7:30","2013-12-29 18:30","2014-01-02 12:30","2014-01-10 18:00","2014-02-13 4:30","2014-02-14 2:30","2014-02-28 16:00",
                     "2014-03-03 16:30","2014-03-05 21:00","2014-03-25 8:30","2014-04-17 16:30","2014-04-24 1:00","2014-05-01 0:00","2014-05-17 16:30",
                     "2014-05-18 12:30","2014-05-23 13:30","2014-05-25 2:00","2014-06-03 6:00","2014-07-06 8:00","2014-08-30 15:00","2014-09-23 23:30",
                     "2014-10-17 12:30","2014-10-19 15:00","2014-11-03 16:30","2014-11-06 7:00","2014-11-25 13:00","2014-12-05 18:00","2014-12-27 9:30"))
  
  end <-  as.Date(c("2009-01-22 15:00", "2009-02-25 12:00", "2009-03-05 9:00","2009-05-17 13:00","2009-07-11 14:30","2009-08-04 1:30",
                    "2009-08-06 2:00","2009-08-11 23:30","2010-01-13 0:00","2010-01-16 5:00","2010-01-19 6:30","2010-01-20 10:00","2010-01-26 3:30",
                    "2010-02-19 0:00","2010-02-28 18:00","2010-03-03 3:00","2010-03-19 16:00","2010-03-31 9:30","2010-04-12 4:00","2010-04-27 17:30",
                    "2010-04-29 8:00","2010-05-01 12:30","2010-05-04 2:30","2010-05-26 6:00","2010-05-28 2:00","2010-06-14 18:30","2010-06-17 5:30",
                    "2010-08-09 0:30","2010-08-28 5:00","2010-08-31 23:30","2010-09-21 6:00","2010-09-26 14:00","2010-11-05 11:00","2010-11-08 17:30",
                    "2010-11-13 13:30","2010-11-16 18:30","2010-11-19 6:00","2010-12-03 10:30","2010-12-31 15:00","2011-01-10 2:30","2011-01-31 23:30",
                    "2011-02-10 6:30","2011-02-20 9:30","2011-02-21 6:00","2011-03-03 11:30","2011-03-27 7:00","2011-04-05 21:30","2011-05-13 23:00",
                    "2011-05-19 11:00","2011-05-20 2:30","2011-05-23 11:00","2011-05-28 22:30","2011-06-12 3:30","2011-06-17 2:30","2011-06-26 8:00",
                    "2011-07-08 3:30","2011-07-13 8:00","2011-07-29 2:30","2011-08-31 23:30","2011-09-30 23:30","2011-10-14 10:00","2011-12-02 21:30",
                    "2012-01-01 2:30","2012-01-07 8:00","2012-01-21 11:00","2012-01-27 7:30","2012-02-05 4:00","2012-02-15 18:30","2012-02-16 20:00",
                    "2012-02-26 20:00","2012-03-07 23:00","2012-04-14 18:00","2012-04-23 7:30","2012-04-26 7:30","2012-05-06 14:00","2012-06-11 6:30","2012-06-14 16:00",
                    "2012-06-19 19:00","2012-07-15 17:30","2012-07-16 15:30","2012-08-12 15:30","2012-08-13 15:00","2012-08-14 14:00","2012-08-20 16:30",
                    "2012-08-23 13:00","2012-08-24 12:30","2012-08-27 16:30","2012-08-29 12:00","2012-10-14 15:30","2012-10-19 6:30","2012-10-22 19:00",
                    "2012-11-03 17:00","2012-11-04 23:00","2012-11-07 13:30","2012-11-27 0:00","2012-12-08 21:00","2012-12-23 22:30","2012-12-27 5:30",
                    "2012-12-31 13:00","2013-01-11 9:00","2013-01-25 5:30","2013-01-28 11:30", "2013-02-23 10:30","2013-02-27 21:00","2013-03-04 17:00",
                    "2013-03-23 23:30","2013-04-05 10:00","2013-04-10 17:30","2013-05-23 16:00","2013-06-07 14:00","2013-06-15 10:30",
                    "2013-07-09 15:30","2013-08-20 14:00","2013-08-31 23:30","2013-09-21 6:30","2013-10-03 19:30","2013-11-03 21:30",
                    "2013-11-10 21:30","2013-11-14 22:30","2013-11-16 6:00","2013-12-23 19:00","2013-12-28 1:00","2013-12-31 15:30",
                    "2014-01-03 18:00","2014-01-15 19:00","2014-02-13 7:00","2014-02-21 4:30","2014-03-02 5:00","2014-03-03 21:00",
                    "2014-03-21 5:30", "2014-04-08 17:00","2014-04-18 3:00", "2014-04-30 15:30","2014-05-15 17:30","2014-05-18 2:00","2014-05-18 19:30",
                    "2014-05-24 16:30","2014-06-02 17:00","2014-06-03 11:00","2014-07-07 14:00","2014-08-31 3:30","2014-09-24 5:30","2014-10-18 3:30",
                    "2014-11-02 5:30","2014-11-04 21:30","2014-11-10 5:30","2014-11-30 15:30","2014-12-26 17:00","2014-12-27 15:00"))
  mean.out <- data.frame(matrix(vector(), length(start),32)) #creating an empty vector
  sum.out <- data.frame(matrix(vector(), length(start),32)) #creating an empty vector 
  max.out <- data.frame(matrix(vector(), length(start),32)) #creating an empty vector
  min.out <- data.frame(matrix(vector(), length(start),32)) #creating an empty vector 
  sd.out <- data.frame(matrix(vector(), length(start),32)) #creating an empty vector 
  duration <- data.frame(matrix(vector(), length(start),32)) #creating an empty vector
  DWM <- data.frame(matrix(vector(), length(start),32)) #creating an empty vector
  
  for (i in 1:length(start)) {
    temp.start <- start[i]
    temp.end <- end[i]
    
    #temp.event <- subset(data, as.Date(data$date) <= temp.start & data$date >= temp.end)
    myfunc <- function(x,y){data[as.Date(data$date) >= x & as.Date(data$date) <= y,]}
    temp.event <- myfunc(temp.start,temp.end)  
    temp.event <- temp.event[,c(2,3,7:11,15:29,49:52,55,58,59,61,63,68)]
    top <- colnames(temp.event)
    
    # calculate stats
    mean.out[i,] <- apply(temp.event, 2, function(x) mean(x, na.rm = TRUE))
    colnames(mean.out) <- paste(top, "mean", sep = "_")
    #sum.out[i,] <- apply(temp.event, 2, function(x) sum(x, na.rm = TRUE))
    #colnames(sum.out) <- paste(top, "sum", sep = "_")
    max.out[i,] <- apply(temp.event, 2, function(x) max(x, na.rm = TRUE))
    colnames(max.out) <- paste(top, "max", sep = "_")
    min.out[i,] <- apply(temp.event, 2, function(x) min(x, na.rm = TRUE))
    colnames(min.out) <- paste(top, "min", sep = "_")
    sd.out[i,] <- apply(temp.event, 2, function(x) sd(x, na.rm = TRUE))
    colnames(sd.out) <- paste(top, "sd", sep = "_")
    duration[i,] <- temp.end - temp.start
    
    ## Calculate the flow weighted mean - as per {Inamdar:2011bo}
    # Cw = sum(c*Q)/sum(Q)
    DWM[i,1] <- sum((temp.event$FI*temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,2] <- sum((temp.event$HIX_ohno_area*temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,3] <- sum((temp.event$FrI*temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,4] <- sum((temp.event$peakA *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,5] <- sum((temp.event$peakC *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,6] <- sum((temp.event$peakB *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,7] <- sum((temp.event$peakT *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,8] <- sum((temp.event$C1.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,9] <- sum((temp.event$C2.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,10] <- sum((temp.event$C3.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,11] <- sum((temp.event$C4.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,12] <- sum((temp.event$C5.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,13] <- sum((temp.event$C6.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,14] <- sum((temp.event$C7.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,15] <- sum((temp.event$C8.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,16] <- sum((temp.event$C9.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,17] <- sum((temp.event$C10.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,18] <- sum((temp.event$C11.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,19] <- sum((temp.event$C12.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,20] <- sum((temp.event$C13.x *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,21] <- sum((temp.event$perprotein *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,22] <- sum((temp.event$redox *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,23] <- sum((temp.event$C1_per *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,24] <- sum((temp.event$C2_per *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,25] <- sum((temp.event$C3_per *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,26] <- sum((temp.event$C4_per *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,27] <- sum((temp.event$Cl_mgL *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,28] <- sum((temp.event$NO2_N_mgL *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,29] <- sum((temp.event$NO3_N_mgL *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    DWM[i,30] <- sum((temp.event$SO4_S_mgL *temp.event$Q.L.s), na.rm = TRUE)/sum(temp.event$Q.L.s, na.rm = TRUE)
    colnames(DWM) <- c("DWM_Fi","DWM_HIX", "DWM_FrI", "DWM_peakA", "DWM_peakC", "DWM_peakB","DWM_peakT","DWM_CMC1", "DWM_CMC2",
                       "DWM_CMC3", "DWM_CMC4","DWM_CMC5", "DWM_CMC6","DWM_CMC7", "DWM_CMC8","DWM_CMC9", "DWM_CMC10",
                       "DWM_CMC11", "DWM_CMC12","DWM_CMC13", "DWM_perprotein","DWM_redox", "DWM_C1per","DWM_C2per",
                       "DWM_C3per","DWM_C4per","DWM_Cl","DWM_NO2","DWM_NO3","DWM_SO4")
    # hysterisis calculations?
  }
  event.stats.grab <- cbind(mean.out, sum.out, max.out, min.out, sd.out, DWM ,duration)
  event.stats.grab$start <- start
  event.stats.grab$end <- end
  
  return(event.stats)
}