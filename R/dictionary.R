#
#
#
compute.dictionary.line.info<-function(dictionary.list, parser.tolerance, overlapping.linesQ=T){
  
  line.info.list<-rep(list(NULL),length(dictionary.list))
  line.totals<-rep(NA,length(dictionary.list))
  
  for(i in 1:length(dictionary.list)){
    
    line.info.list[[i]]<-find.line.idxs(dictionary.list[[i]],tol=parser.tolerance,plotQ=F,printQ=F) #Holds the points of the signals defining the lines
    ref.line.idxs<-line.info.list[[i]][,1]
    
    if(overlapping.linesQ<-T){
      num.lines<-length(ref.line.idxs)-2
    }
    if(overlapping.lins<-F){
      num.lines<-floor(length(ref.line.idxs)/2)
    }
    line.totals[i]<-num.lines
  }
  score.intervals<-cbind(c(1,cumsum(line.totals)+1)[-(length(cumsum(line.totals))+1)], cumsum(line.totals))
  
  grp.ids.for.lines<-rep(0,sum(line.totals))
  for(i in 1:length(line.totals)){
    tmp<-rep(i,line.totals[i])
    grp.ids.for.lines[score.intervals[i,1]:score.intervals[i,2]]<-tmp
  }
  
  all.line.info<-list(line.info.list, line.totals, score.intervals, grp.ids.for.lines)
  return(all.line.info)
  
}

#
#
#
biasotti.murdock.fv<-function(dictionary.list, query.signal, line.score.tolerances=c(0.8,50,0.1), line.info, 
                              renormalize.linesQ=F, max.lag=2000, rackQ=T, overlapped.linesQ=T, printQ=F){
  
  line.extrema.list<-line.info[[1]] #Peak and valley indices for all the signals in the dictionary
  line.totals<-line.info[[2]]       #Number of lines in each group (i.e. word, element of the dictionary)
  score.intervals<-line.info[[3]]   #Start and stop index of each element in the doctionary
  
  total.mean.scores.vec<-rep(-10,sum(line.totals)) #Initalize FV with impossible values. Check for them at the end
  
  #Assumes dictionary is in list form!
  for(i in 1:length(dictionary.list)){
    
    dict.vec<-dictionary.list[[i]] ##Grab one of the group references (words) from the BM-dictionary. There should be no NAs in these
    ref.line.idxs<-line.extrema.list[[i]] #Info for what kind of extreama each of the critical points are. Needed for merge
    
    if(names(dev.cur()) != "null device") { #Hack to help prevent fail on the merge function
      dev.off()
    }
    
    if(printQ==TRUE){
      print(paste("Starting query vec<->dictionary element",i,"computation."))
    }
    line.scores.mat<-merge(dict.vec, query.signal, ref.line.idxs, score.tols, renormalize.linesQ=renormalize.linesQ, max.lag=max.lag, rackQ=rackQ, overlapped.linesQ=overlapped.linesQ, plotQ=F)
    mean.scores.vec<-line.scores.mat[,4] #FV portion for qv<->dict j
    total.mean.scores.vec[score.intervals[i,1]:score.intervals[i,2]]<-mean.scores.vec    
  }
  
  return(total.mean.scores.vec)
}