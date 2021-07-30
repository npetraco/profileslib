#-----------------------------------------------------
#Determine the word indices for each chapter in a dictionary
#-----------------------------------------------------
parse.dictionary <- function(dictionary.mat, overlapping.lins=F, parse.tolerance=0.001){
  
  line.info.list <- rep(list(NULL),nrow(dictionary.mat))
  line.totals <- rep(NA,nrow(dictionary.mat))
  
  #Loop over the chapters in the dictionary and determine how many/indices of the words
  for(i in 1:nrow(dictionary.mat)){
    
    line.info.list[[i]] <- find.line.idxs( as.numeric( na.omit( dictionary.mat[i,] ) ) , tol=parse.tolerance, plotQ=F, printQ=F)
    
    ref.line.idxs <- line.info.list[[i]][,1]
    
    if(overlapping.lins == T){
      num.lines <- length(ref.line.idxs)-2
    }
    
    if(overlapping.lins == F){
      num.lines <- floor((length(ref.line.idxs)-1)/2) #IS THIS CORRECT ???????
    }
    line.totals[i] <- num.lines
  }
  
  #Indices of the words in the entire dictionary
  #Rows are the chapter numbers
  #Cols list the index span of the words in each chapter
  score.intervals <- cbind(c(1,cumsum(line.totals)+1)[-(length(cumsum(line.totals))+1)], cumsum(line.totals))
  
  dictionary.words.info <- list(line.info.list, line.totals, score.intervals)
  names(dictionary.words.info) <- c("available word info","# lines in chapters", "word indices for chapters")
  
  return(dictionary.words.info)
  
}


#-----------------------------------------------------
#Generate a Biasotti-Murdock feature vector for a signal
#-----------------------------------------------------
gen.bm.fv.ref <- function(signal, dictionary.mat, dictionary.words.info, overlapped.linesQ=F) {
  
  line.info.list <- dictionary.words.info[[1]]  # "available word info", i.e. critical values and extrema types
  line.totals <- dictionary.words.info[[2]]     # "# lines in chapters"
  #score.intervals <- dictionary.words.info[[3]] # "word indices for chapters"
  
  #fvec <- rep(-10,sum(line.totals))
  fvec <- NULL
  
  tmp.count<-0
  for(i in 1:nrow(dictionary.mat) ){ #Loop over the dictionary
    
    rv <- as.numeric( na.omit( dictionary.mat[i,] ) )
    qv <- as.numeric( na.omit( signal ) )             #Do this a better way?
    ref.line.idxs<-line.info.list[[i]]
        
    #-------------------------------------- Generate scores on the chapter
    if(!(names(dev.cur()) == "null device")){
      dev.off()
    }
    chapter.score.vec <- merge(rv, qv, 
                              ref.line.idxs, 
                              score.tolerance=c(0.8,50,0.1), 
                              max.lag=2000, 
                              rackQ=TRUE, 
                              renormalize.linesQ=FALSE, 
                              overlapped.linesQ=overlapped.linesQ, 
                              plotQ=F)[,4]
    #print(chapter.score.vec)
    print(paste("Chapter ", i, " done."))
    #--------------------------------------    
    #fvec[score.intervals[i,1]:score.intervals[i,2]] <- chapter.score.vec #This should only work with overlapping lines
    fvec <- c(fvec,chapter.score.vec)
  }
  
  return(fvec)
  
}

#-----------------------------------------------------
#Generate a Biasotti-Murdock feature vector for a signal. Includes plot option
#-----------------------------------------------------
gen.bm.fv.ref2 <- function(signal, dictionary.mat, dictionary.words.info, overlapped.linesQ=F, plotQ=FALSE) {
  
  line.info.list <- dictionary.words.info[[1]]  # "available word info", i.e. critical values and extrema types
  line.totals <- dictionary.words.info[[2]]     # "# lines in chapters"
  #score.intervals <- dictionary.words.info[[3]] # "word indices for chapters"
  
  #fvec <- rep(-10,sum(line.totals))
  fvec <- NULL
  
  tmp.count<-0
  for(i in 1:nrow(dictionary.mat) ){ #Loop over the dictionary
    
    rv <- as.numeric( na.omit( dictionary.mat[i,] ) )
    qv <- as.numeric( na.omit( signal ) )             #Do this a better way?
    ref.line.idxs<-line.info.list[[i]]
    
    #-------------------------------------- Generate scores on the chapter
    if(!(names(dev.cur()) == "null device")){
      dev.off()
    }
    chapter.score.vec <- merge(rv, qv, 
                               ref.line.idxs, 
                               score.tolerance=c(0.8,50,0.1), 
                               max.lag=2000, 
                               rackQ=TRUE, 
                               renormalize.linesQ=FALSE, 
                               overlapped.linesQ=overlapped.linesQ, 
                               plotQ=plotQ)[,4]
    #print(chapter.score.vec)
    print(paste("Chapter ", i, " done."))
    #--------------------------------------    
    #fvec[score.intervals[i,1]:score.intervals[i,2]] <- chapter.score.vec #This should only work with overlapping lines
    fvec <- c(fvec,chapter.score.vec)
  }
  
  return(fvec)
  
}

#-----------------------------------------------------
#Generate a Biasotti-Murdock feature vector for a signal. Optimized merge function.
#-----------------------------------------------------
gen.bm.fv.fast <- function(signal, dictionary.mat, dictionary.words.info, overlapped.linesQ=F, score.tolerances) {
  
  line.info.list <- dictionary.words.info[[1]]  # "available word info", i.e. critical values and extrema types
  line.totals <- dictionary.words.info[[2]]     # "# lines in chapters"
  #score.intervals <- dictionary.words.info[[3]] # "word indices for chapters"
  
  #fvec <- rep(-10,sum(line.totals))
  fvec <- NULL
  
  tmp.count<-0
  for(i in 1:nrow(dictionary.mat) ){ #Loop over the dictionary
    
    rv <- as.numeric( na.omit( dictionary.mat[i,] ) )
    qv <- as.numeric( na.omit( signal ) )             #Do this a better way?
    ref.line.idxs<-line.info.list[[i]]
    
    #-------------------------------------- Generate scores on the chapter
    chapter.score.vec <- merge5a(rv, qv, 
                                ref.line.idxs, 
                                score.tolerance=score.tolerances, 
                                max.lag=2000, 
                                rackQ=TRUE, 
                                renormalize.linesQ=FALSE, 
                                overlapped.linesQ=overlapped.linesQ)[,5]
    #print(paste("Chapter ", i, " done."))
    #--------------------------------------    
    #fvec[score.intervals[i,1]:score.intervals[i,2]] <- chapter.score.vec #This should only work with overlapping lines
    fvec <- c(fvec,chapter.score.vec)
  }
  
  return(fvec)
  
}