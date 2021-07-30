#---------------------------------------------------------------
#Generate an IBIS like "hit-list" of possible signal I.D.s using the ccf
#---------------------------------------------------------------
get.ccf.hit.list<-function(signal.dictionary, query.signal, lagmax, cut.off.val, printQ=FALSE){
  
  maxccf.scores<-rep(NA,nrow(signal.dictionary))
  
  for(i in 1:nrow(signal.dictionary)){
    rv <- as.numeric(na.omit(signal.dictionary[i,]))
    qv <- as.numeric(na.omit(query.signal))
    maxccf.score <- ccfmax.distance(qv,rv,maxlag=lagmax,printQ=F)[1]
    maxccf.scores[i] <- maxccf.score
    
    if(printQ==TRUE){
      print(paste("Grp:",i,"MCCF:",maxccf.score))
    }
  }
  #sort(maxccf.scores,decreasing=T)[1:10]
  
  poss.id.idxs<-which(maxccf.scores>cut.off.val)
  
  return(poss.id.idxs)
  
}