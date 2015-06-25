## Word frequency analysis
# aim is to take the submissions, look at the most frequenct words, and try and use this to filter submissions that are the same

# Make a list of submissions whose top 10 most frequent words are the same?
# references: http://onepager.togaware.com/TextMiningO.pdf
# https://gist.github.com/benmarwick/11333467
# https://jhuria.wordpress.com/2012/07/01/text-mining-in-r/
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
library(wordcloud)
library(SnowballC)
library(reshape)
library(plyr)
library(gsubfn)

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

dname <-  "/Users/ashlee/Documents/WaterAct Paper/Selected_individual/corpus"

y <- length(dir(dname)) #should be number of text files that you pasted here

# get file names
setwd(dname)
filelist_txt <- list.files(pattern = ".txt$")

sample.ID <- 0 #create sample ID variable

for (i in 1:y){
  sample.ID.temp <- strapplyc(filelist_txt[i], "(.*).txt", simplify = TRUE)
  sample.ID[i] <- sample.ID.temp
}


############## Creating corpus file
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

#which words are associated with what? can find context of words
#protect.a <- findAssocs(tdm, "protect",.8)

freq <- colSums(as.matrix(dtm))

ord <- order(freq)

# Least frequent terms
least <- freq[head(ord)]

# Most frequent terms
most <- freq[tail(ord)]

# calculate the frequency of words
wordfreq <- sort(rowSums(as.matrix(tdm)), decreasing=TRUE)

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
top.wc <- sort(table(wordc[,1]),decreasing=TRUE)[1:5]
top.wcn <- as.numeric(rownames(top.wc))

# get the top 5 groups of submissions sorted according to word count 
wc.1 <- rownames(subset(wordc, wordc[,1] >= (top.wcn[1]-5) & wordc[,1] <= (top.wcn[1]+5)))
wc.2 <- rownames(subset(wordc, wordc[,1] >= (top.wcn[2]-5) & wordc[,1] <= (top.wcn[2]+5)))
wc.3 <- rownames(subset(wordc, wordc[,1] >= (top.wcn[3]-5) & wordc[,1] <= (top.wcn[3]+5)))
wc.4 <- rownames(subset(wordc, wordc[,1] >= (top.wcn[4]-5) & wordc[,1] <= (top.wcn[4]+5)))
wc.5 <- rownames(subset(wordc, wordc[,1] >= (top.wcn[5]-5) & wordc[,1] <= (top.wcn[5]+5)))

#total <- cbind(wc.1, wc.2, wc.3, wc.4, wc.5)

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

###### use function to save files into different folders
test.1 <- save.wc(listwc = wc.1, sub.wc = 1) #copy and do it for everyone

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
    wordfreq.50 <- rownames(as.matrix(sort(rowSums(as.matrix(tdm.1)), decreasing=TRUE)))[1:50]
    
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

## use function to write matrix of most frequent works  to separate csv files

num1 <- wfa.matrix(parentdir = dname , section = 5)