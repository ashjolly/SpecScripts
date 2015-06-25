#
#
#
# Function for taking files from a list of files and saving into another directory
# note that this takes the files from your set wd
# 25 june2015, from water act paper
#######

save.wc <- function(listwc, sub.wc) {
  y = length(listwc)
  for (i in 1:y){
    filename <- paste(toString(listwc[i]), ".txt", sep = "")
    filepath.temp <- file.path(dname, paste("corpus_", sub.wc, "/", filename, sep = ""))
    temp <- file.copy(filename, filepath.temp)
  }
  return(temp)
}