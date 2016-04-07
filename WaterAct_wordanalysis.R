## Word frequency analysis
# aim is to take the submissions, look at the most frequenct words, and try and use this to filter submissions that are the same

# Make a list of submissions whose top 10 most frequent words are the same?
# references: http://onepager.togaware.com/TextMiningO.pdf
# https://gist.github.com/benmarwick/11333467
# https://jhuria.wordpress.com/2012/07/01/text-mining-in-r/
# 22June2015
# Ashlee Jollymore
################

## set working directory
rm(list = ls())
ls()

Selected_individual <- "/Users/user/Dropbox/PhD Work/WaterAct Paper/Selected_individual"
dname <-  "/Users/user/Dropbox/PhD Work/WaterAct Paper/Selected_individual/corpus"
# where txt files will be stored

setwd(Selected_individual)

########## Convert PDF to txt files and input as a tm object (corpus)
# Necessary packages
library(tm)
library(wordcloud)
library(SnowballC)
library(reshape)
library(plyr)
library(gsubfn)

############### Convert PDFs to txt files
# Necessary to do this prior to doing text analysis. 
# Not necessary if your files are already in .txt files
# read PDFs from file with randomly selected files
length(dir(Selected_individual))

# make a vector of PDF file names
myfiles <- list.files(path = Selected_individual, pattern = "pdf",  full.names = TRUE)

# convert each PDF file that is named in the vector into a text file 
# text file is created in the same directory as the PDFs
# puts a bunch of text files into the destination

#lapply(myfiles, function(i) system(paste('"/usr/local/bin/pdftotext"', paste0('"', i, '"')), wait = FALSE) )

# get text files that you just created in that directory. Not that .txt files are moved to 'Corpus' 
# folder within the same folder

# get file names
setwd(dname)
filelist_txt <- list.files(pattern = ".txt$")
y = length(filelist_txt)

sample.ID <- 0 #create sample ID variable

for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_txt[i], "(.*).txt", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}

############## Creating corpus file
# First step of text analysis
# input the text files into a corpus
corpus.i <- Corpus(DirSource(dname), readerControl = list(language="lat"))

# inspect corpus to make sure that the documents are input
# Should pull up one of the submissions
inspect(corpus.i[2])

#MyCorpus <- tm_map(YourCorpus,
#                   content_transformer(function(x) iconv(x, to='UTF-8-MAC', sub='byte')),
#                   mc.cores=1)

#test <- tm_map(YourCorpus, stemDocument, lazy = TRUE)

######## Preprocessing text data prior to analysis
# To see transformations possible within the tm package  -> getTransformations()
# using two cores! tell r to nly use one core using lazy = TRUE
# this was converted into a function so that could use this later for partitioned txt files

text.pro <- function(txtdoc){
  
  corpus <- tm_map(txtdoc, content_transformer(tolower))

  # get rid of weird punctuation - change it to a space
  toSpace <- content_transformer(function(x, pattern) gsub(pattern, " ", x))
  docs <- tm_map(corpus, toSpace, "/|@|\\|***|", lazy = TRUE)

  # convert all upper case to lower case
  docs <- tm_map(corpus, content_transformer(tolower))

  # remove numbers
  docs <- tm_map(docs, content_transformer(removeNumbers))

  # remove punctuation
  docs <- tm_map(docs, content_transformer(removePunctuation))

  # remove stop words like for, very, and, of, are, plus personal stop words
  docs <- tm_map(docs, removeWords, c(stopwords("english"), 
     "personal","identifiers","removed","water","wsa","sustainability","act", "proposal", "need"),
      lazy = TRUE)

  docs <- tm_map(docs, removeWords, c(stopwords("english"),"my","custom","words")) 

  # strip white spaces  
  docs <- tm_map(docs, stripWhitespace, lazy = TRUE)

  # stem document - remove common word endings
  docs <- tm_map(docs, stemDocument, lazy = TRUE)

  # convert back to plan text document
  docs <- tm_map(docs, PlainTextDocument, lazy = TRUE)

  return(docs)
}
# use this function to clean corpus from all individual txt submissions
docs <- text.pro(txtdoc = corpus.i)

########## Create document term matrix
# A document term matrix is simply a matrix with documents as the rows and terms as 
# the columns and a count of the frequency 
# of words as the cells of the matrix. We use DocumentTermMatrix() to create the matrix

#create document term matrix
dtm <- DocumentTermMatrix(docs)
tdm <- TermDocumentMatrix(docs)

########### Frequency Analysis on all data!
#Frequent Terms and Associations
#freq.terms.1 <- findFreqTerms(dtm, lowfreq=4)
freq.terms <- findFreqTerms(tdm, lowfreq=100)

freq <- colSums(as.matrix(dtm))

ord <- order(freq)

# Least frequent terms
least <- freq[head(ord)]

# Most frequent terms
most <- freq[tail(ord)]

# calculate the frequency of words
wordfreq <- sort(rowSums(as.matrix(tdm)), decreasing=TRUE)

# #which words are associated with what? can find context of words
# according to frequency- most frequent
top.associations <- findAssocs(dtm, c("protect", "use", "must", "public", "new", "resourc"), 
                               c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4))
env.associations <- findAssocs(dtm, c("environment"), 0.4)
                               
# save

dput(top.associations, file = file.path(paste(dname, "/top_associations.txt", sep=""))) 
dput(env.associations, file = file.path(paste(dname, "/env_associations.csv", sep="")))
    #clustering :: k-means clustering
cluster <- kmeans(tdm, 10)
#colnames(cluster) <- sample.ID

############### Subsetting form data
########### Word Count Analysis
# sort submissions by the word count, on the assumption that forms will have similar word count
# After preliminary sorted by word count, will confirm that it is a form using word frequency analysis
# first sort, and then find qord frequencies of groups. 
# then compare the word frequency of inidvidual submissions to the groups

# Word count per document
wordc <- data.frame(rowSums(as.matrix(dtm)))
row.names(wordc) <- sample.ID

# plot to see if there are distributions of word counts to tease out forms
hist(wordc[,1], breaks = 200)

# get most frequent word counts
num.groups <- 20 #number of word count groups you want to create
top.wc <- sort(table(wordc[,1]),decreasing=TRUE)[1:num.groups]
top.wcn <- as.numeric(rownames(top.wc))

# get the submissions sorted according to word count 
# adaptable to choose the number of groups

###########
# save into different folders - function
save.wc <- function(listwc, sub.wc) {
  y = length(listwc)
  for (i in 1:y){
    filename <- paste(toString(listwc[i]), ".txt", sep = "")
    filepath.temp <- file.path(dname, paste("corpus_", sub.wc, "/", filename, sep = ""))
    temp <- file.copy(filename, filepath.temp)
  }
  return(temp)
}

###########
# use to sort and save the files into separate folders

temp <- seq(1,num.groups, by = 1)
n = length(temp)

for (i in 1:n) {
  # partition files
  wc.temp <- rownames(subset(wordc, wordc[,1] >= (top.wcn[i]-5) & wordc[,1] <= (top.wcn[i]+5)))
  #name each of the divisions -wc.[the partition you are looking for]
  assign(paste("wc.", i, sep = ""), wc.temp) 
  
  ## save files
  #create folder if it doesn't exist
  folder.name <- file.path(dname, paste("corpus_", i,sep = ""))
  
  # if the folder doesn't exist, create it.
  if (file.exists(folder.name) == FALSE) {
    dir.create(folder.name, showWarnings = FALSE)
  }
  
  #save files in correct folder                         
  save.wc(listwc = wc.temp, sub.wc = i)
}

############### Word frequency analysis
######## Test that the word count partitioned the files well

wfa.matrix <- function(parentdir, section){
    
    filepath.1 <- file.path(parentdir, paste("corpus_", section, "/", sep = ""))
    corpus.1 <- Corpus(DirSource(filepath.1), readerControl = list(language="lat"))
    # massage data
    corpus.1 <- text.pro(txtdoc = corpus.1 )

    #create document term matrix
    dtm.1 <- DocumentTermMatrix(corpus.1)
    tdm.1 <- TermDocumentMatrix(corpus.1)
    freq.terms <- findFreqTerms(tdm.1, lowfreq=100)

    # calculate the frequency of words
    # names of the top 50 most common names within the subset
    if (length(rownames(as.matrix(sort(rowSums(as.matrix(tdm.1)), decreasing=TRUE)))) >= 50) {
      wordfreq.50 <- rownames(as.matrix(sort(rowSums(as.matrix(tdm.1)), decreasing=TRUE)))[1:50]
    }
    
    if (length(rownames(as.matrix(sort(rowSums(as.matrix(tdm.1)), decreasing=TRUE)))) < 50) {
        wordfreq.50 <- rownames(as.matrix(sort(rowSums(as.matrix(tdm.1)), decreasing=TRUE)))
    }
    
    # convert tdm to matrix
    term.matrix <- as.matrix(tdm.1)

    #select only the rows with the 50 most common names in the subset to simplify the tdm matrix
    test <- as.matrix(term.matrix[wordfreq.50,])

    #attach filenames as the column names of the matrix
    setwd(filepath.1)
    sample.ID <- 0 #reset sample.ID variable
    filelist_txt <- list.files(pattern = ".txt$")
    y <- length(dir(filepath.1))
      for (i in 1:y){
        sample.ID.temp <- strapplyc(filelist_txt[i], "(.*).txt", simplify = TRUE)
        sample.ID[i] <- sample.ID.temp
      }
    colnames(test) <- sample.ID

    # write matrix to compare
    directory.matrix <- file.path(Selected_individual, paste("form_splitresults/", section, ".csv", sep = ""))
    write.table(test, file = directory.matrix, row.names = TRUE, col.names = TRUE, sep = ",")
}

### Use function to write matrix of most frequent works to separate csv files
for (t in 1:n) {
  wfa.matrix(parentdir = dname , section = t)
}
