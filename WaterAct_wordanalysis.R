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

########## Convert PDF to txt files and input as a tm object (corpus)
# Necessary packages
library(tm)
library(SnowballC)

# read PDFs from file with randomly selected files
length(dir(Selected_individual))

# make a vector of PDF file names
myfiles <- list.files(path = Selected_individual, pattern = "pdf",  full.names = TRUE)

# convert each PDF file that is named in the vector into a text file 
# text file is created in the same directory as the PDFs
#puts a bunch of text files into the destination

#lapply(myfiles, function(i) system(paste('"/usr/local/bin/pdftotext"', paste0('"', i, '"')), wait = FALSE) )

# get text files that you just created in that directory. Not that .txt files are moved to 'Corpus' 
# folder within the same folder

#dname <- file.path(".", "corpus")
dname <-  "/Users/ashlee/Documents/WaterAct Paper/Selected_individual/corpus"

length(dir(dname)) #should be number of text files that

# input the text files into a corpus
YourCorpus <- Corpus(DirSource(dname), readerControl = list(language="lat"))

# inspect corpus to make sure that the documents are input
# Should pull up one of the submissions
inspect(YourCorpus[2])

MyCorpus <- tm_map(YourCorpus,
                   content_transformer(function(x) iconv(x, to='UTF-8-MAC', sub='byte')),
                   mc.cores=1)

######## Preprocessing text data prior to nalysis
# To see transformations possible within the tm package  -> getTransformations()
# using two cores! tell r to nly use one core using lazy = TRUE

# get rid of weird punctuation - change it to a space
docs <- tm_map(MyCorpus, toSpace, "/|@|\\|***|", lazy=TRUE)

# convert all upper case to lower case
docs <- tm_map(docs, content_transformer(tolower), lazy=TRUE)

# remove numbers
docs <- tm_map(docs, removeNumbers, lazy=TRUE)

# remove punctuation
docs <- tm_map(docs, removePunctuation, lazy=TRUE)

# remove stop words like for, very, and, of, are...
docs <- tm_map(docs, removeWords, stopwords("english"), lazy=TRUE)

# add own stop words here
# ***PERSONAL IDENTIFIERS REMOVED***
docs <- tm_map(docs, removeWords, c("***personal identifiers removed***", "wsa", "water sustainability act"), lazy = TRUE)

# strip white spaces from 
docs <- tm_map(docs, stripWhitespace, lazy = TRUE)

########## Create document term matrix
# A document term matrix is simply a matrix with documents as the rows and terms as the columns and a count of the frequency 
# of words as the cells of the matrix. We use DocumentTermMatrix() to create the matrix

#create document term matrix
dtm <- DocumentTermMatrix(docs)


### word count


