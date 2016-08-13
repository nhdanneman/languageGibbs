

# Read in english words
en <- as.character(read.delim("/Users/nathan.danneman/Documents/ndgit/languageGibbs/englishWords.txt", header=F)[,1])
fr <- as.character(read.delim("/Users/nathan.danneman/Documents/ndgit/languageGibbs/frenchWords.txt", header=F)[,1])

# Change all letters to unaccented ones:
unwantedArray <- list(    'Š'='S', 'š'='s', 'Ž'='Z', 'ž'='z', 'À'='A', 'Á'='A', 'Â'='A', 'Ã'='A', 'Ä'='A', 'Å'='A', 'Æ'='A', 'Ç'='C', 'È'='E', 'É'='E',
                          'Ê'='E', 'Ë'='E', 'Ì'='I', 'Í'='I', 'Î'='I', 'Ï'='I', 'Ñ'='N', 'Ò'='O', 'Ó'='O', 'Ô'='O', 'Õ'='O', 'Ö'='O', 'Ø'='O', 'Ù'='U',
                          'Ú'='U', 'Û'='U', 'Ü'='U', 'Ý'='Y', 'Þ'='B', 'ß'='Ss', 'à'='a', 'á'='a', 'â'='a', 'ã'='a', 'ä'='a', 'å'='a', 'æ'='a', 'ç'='c',
                          'è'='e', 'é'='e', 'ê'='e', 'ë'='e', 'ì'='i', 'í'='i', 'î'='i', 'ï'='i', 'ð'='o', 'ñ'='n', 'ò'='o', 'ó'='o', 'ô'='o', 'õ'='o',
                          'ö'='o', 'ø'='o', 'ù'='u', 'ú'='u', 'û'='u', 'ý'='y', 'ý'='y', 'þ'='b', 'ÿ'='y' )
out <- NULL
for(j in 1:length(fr)){
  word <- fr[j]
  for(i in seq_along(unwantedArray)){
     word <- gsub(names(unwantedArray)[i],unwantedArray[i],word)
  }
out <- c(out, word)
}

fr <- out

# keep only those with at least 4 letters:
enFinal <- en[sapply(en,nchar) >= 4]
frFinal <- fr[sapply(fr,nchar) >= 4]

# balanced data sets, not necessary but nice
enFinal <- enFinal[1:length(frFinal)]

# from here, pretend we don't know which words are from which distribution
allwords <- c(enFinal, frFinal)

# need a set of all possible bigrams, including "stop"
allbigrams <- NULL
for(i in 1:length(allwords)){
  word <- allwords[i]
  for(j in 2:nchar(word)-1){
    allbigrams <- c(allbigrams, substring(word, j, j+1))
  }
  allbigrams <- c(allbigrams, paste(substring(word, nchar(word), nchar(word)+1), "END", sep=""))
}
allPossibleBigrams <- unique(allbigrams)  # 368 of them out of 702 logically possible

## Gibbs sampler:

iterations <- 200

# place to hold results:
samples <- matrix(0, nrow=iterations, ncol=length(allwords) )

# initialize the first 'type' of each word randomly
samples[1,] <-  rbinom(length(allwords), 1, 0.5)

# represent each word as a vector of its bigrams, all in a list
allWordList <- list()
for(i in 1:length(allwords)){
  word <- allwords[i]
  holder <- NULL
  for(j in 2:nchar(word)-1){
    holder <- c(holder, substring(word, j, j+1))
  }
  holder <- c(holder, paste(substring(word, nchar(word), nchar(word)+1), "END", sep=""))
  allWordList[[i]] <- holder
}

# actual start of sampler
for(iter in 2:iterations){
  if(iter / 20 == floor(iter/20)){print(paste("Iteration", iter, sep=" " ))}

  # get counts of every bigram of each type:
  bigrams0 <- unlist(allWordList[which(samples[iter-1,]==0)])
  bigrams1 <- unlist(allWordList[which(samples[iter-1,]==1)])
  
  # get distribution
  bigDistr0 <- table(bigrams0)/length(bigrams0)
  bigDistr1 <- table(bigrams1)/length(bigrams1)
  
  # for each word, get its probability of being each type (prod of grams), and sample:
  for(i in 1:length(allWordList)){
    grams <- allWordList[[i]]
    pr0perGram <-  NULL
    pr1perGram <- NULL
    names0 <- names(bigDistr0)
    names1 <- names(bigDistr1)
    for(j in 1:length(grams)){
      if(grams[[j]] %in% names0 ){
        pr0perGram <- c(pr0perGram, bigDistr0[grams[j]])
      } else {
        pr0perGram <- c(pr0perGram, min(bigDistr0)*0.1)
      }
      if(grams[[j]] %in% names1 ){
        pr1perGram <- c(pr1perGram, bigDistr1[grams[j]])
      } else {
        pr1perGram <- c(pr1perGram, min(bigDistr1)*0.1)
      }
    }
    pr0 <- prod(pr0perGram)
    pr1 <- prod(pr1perGram)
    samples[iter, i] <- rbinom(1,1, (pr1/(pr0+pr1)))
  }

} # end of sampler

est <- samples[iterations,]
truth <- c(rep(0, length(enFinal)), rep(1, length(frFinal)))
table(est, truth)


s <- samples[150:iterations,]
diff <- which(apply(s, 2, mean) > 0.4 & apply(s, 2, mean) < 0.5)
allwords[sample(diff, 25, F)]

zeros <- which(apply(s, 2, mean) < 0.05)
allwords[sample(zeros, 25, F)]

ones <- which(apply(s, 2, mean) > 0.95)
allwords[sample(ones, 25, F)]
