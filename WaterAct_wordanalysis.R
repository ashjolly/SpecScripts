## Word frequency analysis
# aim is to take the submissions, look at the most frequenct words, and try and use this to filter submissions that are the same

# Make a list of submissions whose top 10 most frequent words are the same?
# references: http://onepager.togaware.com/TextMiningO.pdf
# https://gist.github.com/benmarwick/11333467
# 22June2015
################

## set working directory
rm(list = ls())
ls()

Selected_individual <- "/Users/ashlee/Documents/WaterAct Paper/Selected_individual"

setwd(Selected_individual)

#install.packages('tm')
library(tm)

# read PDFs from file with randomly selected files
length(dir(Selected_individual))

# make a vector of PDF file names
myfiles <- list.files(path = Selected_individual, pattern = "pdf",  full.names = TRUE)

# convert each PDF file that is named in the vector into a text file 
# text file is created in the same directory as the PDFs

#lapply(myfiles, function(i) system(paste('"/Applications/xpdf/bin64/pdftotext.exe"', 
#                                         paste0('"', i, '"')), wait = FALSE) )

# read PDF files into one 
#docs <- Corpus(DirSource(Selected_individual), readerControl=list(reader=readPDF))

test <- lapply(myfiles, function(i) system(paste('"usr/bin/pdftotext.exe"', paste0('"', i, '"')), wait = FALSE) )

# convert all upper case to lower case
docs <- tm_map(docs, content_transformer(tolower))

# remove numbers
docs <- tm_map(docs, removeNumbers)

# remove punctuation
docs <- tm_map(docs, removePunctuation)

# remove stop words like for, very, and, of, are...
docs <- tm_map(docs, removeWords, stopwords("english"))

# add own stop words here
#docs <- tm_map(docs, removeWords, c("OWN STOP WORDS!!"))

# strip white spaces from 
docs <- tm_map(docs, stripWhitespace)

#create document term matrix
dtm <- DocumentTermMatrix(docs)


### word count


