#-------------------------------------------------------------------------------------
#"Project" lines of a query profile onto reference profile (a chapter of the  dictionary)
#-------------------------------------------------------------------------------------
merge5<-function(reference.signal, query.signal, reference.line.info.mat, score.tolerances, max.lag=10, rackQ=TRUE, renormalize.linesQ=FALSE, overlapped.linesQ=FALSE, plotQ=FALSE){
  
  #Align signals and determine shift to the line indices that will be needed:
  if(rackQ==TRUE){ #Rack the signals over each other with ccf to determine the "best" alignment
    aligned.signals<-align.mov.to.ref.pad.both2(reference.signal,query.signal,lagmax=max.lag,padtyp="NAs",printQ=FALSE)
  }
  if(rackQ==FALSE){ #Manually specify the alignment shift:
    aligned.signals<-manual.align.mov.to.ref.pad.both(reference.signal,query.signal,specified.lag=max.lag,padtyp="NAs",printQ=FALSE)
  }
  
  #Determine how the reference signal indices have shifted with any padding
  reference.idxs<-which(is.na(aligned.signals[[1]])==FALSE)
  
  #Loop over the line indices of the reference and pull out those index spans
  #from both the reference and the query.
  #First correct the line indices due to alignment/padding. Only needs to be done for the reference:
  orig.reference.line.idxs<-reference.line.info.mat[,1] #Line indices in the reference before alignment/padding
  aug.extrema.typ<-reference.line.info.mat[,2]   #For later in plots
  
  padded.reference.line.idxs<-rep(NA,length(orig.reference.line.idxs))
  
  line.widths <- NULL
  corr.vals<-NULL
  dtw.lb_keogh.vals<-NULL
  dtw.lb_improved.vals<-NULL
  dtwDist.vals<-NULL
  dtwNDist.vals<-NULL
  ssqd.vals<-NULL
  mean.scores<-NULL
  norm.scores<-NULL
  mid.line.idxs<-NULL #Used for plotting the summary plot
  
  for(i in 1:length(orig.reference.line.idxs)){
    
    orig.line.idx <- orig.reference.line.idxs[i] #Take into account any padding in the reference
    padded.reference.line.idx <- reference.idxs[orig.line.idx]
    padded.reference.line.idxs[i] <- padded.reference.line.idx
    
  }
  
  padded.reference.signal<-aligned.signals[[1]]
  padded.query.signal<-aligned.signals[[2]]
  if(overlapped.linesQ==TRUE){
    loop.idxs<-seq(1,length(padded.reference.line.idxs),1)  #Idxs for overlapping lines
  } else{
    loop.idxs<-seq(1,length(padded.reference.line.idxs),2) #Idxs for adjacent lines
  }
  
  for(i in 1:(length(loop.idxs))){
    
    left.idx<-padded.reference.line.idxs[loop.idxs[i]]    #Start of the line
    if(overlapped.linesQ==TRUE) {
      right.idx<-padded.reference.line.idxs[loop.idxs[i+2]] #End of the line for for overlapping lines
    } else {
      right.idx<-padded.reference.line.idxs[loop.idxs[i+1]] #End of the line for adjacent lines
    }
    #print("Comming up:")
    #print(loop.idxs[i])
    #print(loop.idxs[i+2])
    #print(left.idx)
    #print(right.idx)
    
    #This is a jerry-rig not to kill the plot window at the last line:
    if(is.na(right.idx)){
      if(plotQ==TRUE){
        inp.end<-readline(prompt="END. Hit any key to quit.") 
      }
      break()
    } else {
      mid.line.idx<-floor(mean(c(right.idx,left.idx))) #This idx will be used for plotting the summary plot below
      mid.line.idxs<-c(mid.line.idxs,mid.line.idx)
    }
    
    reference.line <- padded.reference.signal[left.idx:right.idx] #Line from the reference
    query.line <- padded.query.signal[left.idx:right.idx]         #Line from the query
    
    #Re-normalize lines if desired:
    if(renormalize.linesQ==TRUE){
      reference.line <- norm.profile2(reference.line)
      query.line <- norm.profile2(query.line)
    }
    
    #Replace any NA padding with 0s. Necessary to compare them:
    reference.line[which(is.na(reference.line)==TRUE)] <- 0
    query.line[which(is.na(query.line)==TRUE)] <- 0
    
    line.widths <- c(line.widths,length(reference.line))
    
    #Correlation:
    if(sd(reference.line)==0 | sd(query.line)==0) {
      corr.val<-NA
    } else {
      corr.val<-cor(reference.line,query.line)
    }
    
    if(is.na(corr.val)){
      corr.val<-(-2) #This should trip the flag not to accept this line and not mess up the scoring system below
      #NOTE though, it is NOT a valid correlation value and should not be treated as such
    }
    corr.vals<-c(corr.vals,corr.val)
 
    #DTW:
    wind.leng <- floor(0.05*length(query.line)) #Window size is ~5% of query/reference length
    #print(paste("Line ",i, "  Ref leng: ",length(reference.line), " Qurery leng: ", length(query.line)))
    #print(paste("Window width:",wind.leng))
    if(wind.leng >=5) {
      
      lemire.envelope <- computeEnvelope(query.line,wind.leng)
      query.Upper.lim <- lemire.envelope[,1]
      query.Lower.lim <- lemire.envelope[,2]
      #print(query.Upper.lim)
      #print(query.Lower.lim)
      dtw.lb_keogh.val <- LB_Keogh(reference.line,query.Upper.lim,query.Lower.lim) #Rcpp imp seems faster. Does lots of copying though
      dtw.lb_improved.val <- LB_Improved(reference.line, query.line, query.Upper.lim, query.Lower.lim, wind.leng)
      #print(paste("LB_Keogh: ", dtw.lb_keogh.val))
      #print(paste("LB_Improved: ", dtw.lb_improved.val))
      
      
    } else {
      #print("Skipping lower bound computation. Window too small.")
      dtw.lb_keogh.val <- 1e6
      dtw.lb_improved.val <- 1e6
    }
    
    dtw.lb_keogh.vals <- c(dtw.lb_keogh.vals, dtw.lb_keogh.val)
    dtw.lb_improved.vals <- c(dtw.lb_improved.vals, dtw.lb_improved.val)
    
    #Sum square deviations
    ssqd.val<-sum((query.line-reference.line)^2)
    ssqd.vals<-c(ssqd.val,ssqd.vals)
    ##Sum abs deviations
    #sabd.val<-sum(abs(query.line-reference.line))
    #sabd.vals<-c(sabd.vals,sabd.val)
    ##Mutual Info:
    #mi.val<-1-MIdist(rbind(reference.line,query.line))[[1]]
    #mi.vals<-c(mi.vals,mi.val)
    
    
    #Check scores:
    sc.ck.val<-sum(c(corr.val>=score.tolerances[1], 
                     dtw.lb_keogh.val<=score.tolerances[2], 
                     ssqd.val<=score.tolerances[3], 
                     length(reference.line)>=score.tolerances[4]))
    #print(sc.ck.val)
    #print(sign(corr.val)==1)
    #print(sc.ck.val>=2)
    if(sign(corr.val)==1 & sc.ck.val>=3){ #Corr must be positive and 3 out of 4 preliminary must be to tolerance
      
      #Do full DTW if you made it this far:
      #dtw.info<-dtw(query.line,reference.line)
      #dtwDist.val<-dtw.info$distance
      #dtwDist.vals<-c(dtwDist.vals,dtwDist.val)
      #dtwNDist.val<-dtw.info$normalizedDistance
      #dtwNDist.vals<-c(dtwNDist.vals,dtwNDist.val)
      
      dtw.lb_improved.val <- LB_Improved(reference.line, query.line, query.Upper.lim, query.Lower.lim, wind.leng)
      dtwDist.val <- dtw.lb_improved.val
      dtwNDist.val <- dtw.lb_improved.val
      dtwDist.vals<-c(dtwDist.vals,dtw.lb_improved.val)
      dtwNDist.vals<-c(dtwNDist.vals,dtw.lb_improved.val)
      
      #Check scores again:
      sc.ck.val2<-sum(c(corr.val>=score.tolerances[1], dtwDist.val<=score.tolerances[4], dtwNDist.val<=score.tolerances[5])) 
      if(sign(corr.val)==1 & sc.ck.val2>=2) {
        
        cscr<-(corr.val) #Assumes only positive corr vals mad it through the check
        
        dscr<-(score.tolerances[5]-dtwDist.val)/score.tolerances[5]
        if(dscr<0){ #If DTW dist exceeds the tolerance just set the score to 0 
          dscr<-0   #MAYBE CHANGE THIS TO DOUBLE PENALIZING THE SCORE INSTEAD??
        }
        
        dnscr<-(score.tolerances[6]-dtwNDist.val)/score.tolerances[6]
        if(dnscr<0){ #If Normalized DTW dist exceeds the tolerance just set the score to 0. CHECK THEORY ON NDTW!!!!!!
          dnscr<-0
        }
        
        mnscr<-mean(c(cscr,dscr,dnscr))
        normscr<-sqrt(sum(c(cscr,dscr,dnscr)^2))/sqrt(3)
        
      
      } else {
        
        cscr<-0
        dscr<-0
        dnscr<-0
        mnscr<-0
        normscr<-0
        
        #For the score matrix output
        #DTW:
        dtwDist.val<-NA
        dtwDist.vals<-c(dtwDist.vals,dtwDist.val)
        dtwNDist.val<-NA
        dtwNDist.vals<-c(dtwNDist.vals,dtwNDist.val)
        
      }
            
    } else {
      cscr<-0
      dscr<-0
      dnscr<-0
      mnscr<-0
      normscr<-0
      
      #For the score matrix output
      #DTW:
      dtwDist.val<-NA
      dtwDist.vals<-c(dtwDist.vals,dtwDist.val)
      dtwNDist.val<-NA
      dtwNDist.vals<-c(dtwNDist.vals,dtwNDist.val)
      
    }
    mean.scores<-c(mean.scores,mnscr)
    norm.scores<-c(norm.scores,normscr)
    
    #Plotting section:
    if(plotQ==TRUE){
      
      if(names(dev.cur())=="null device") {
        dev.new() 
      }
      
      inp<-readline(prompt="Hit enter for next plot set, q to quit... ")
      if(inp=="q") {
        print("Quitting Loop")
        break()
      }
      
      print(paste("Line#:",i,"Start:",left.idx,"Stop:",right.idx,"Width:",right.idx-left.idx))
      
      par(mfrow=c(2,1)) #To plot the full signals and the zoomed in lines below
      
      color.vec<-rep(NA,length(aug.extrema.typ))
      for(j in 1:length(color.vec)){
        if(aug.extrema.typ[j]==(-1)){
          color.vec[j]<-"blue"
        }
        if(aug.extrema.typ[j]==1){
          color.vec[j]<-"red"
        }
        if(aug.extrema.typ[j]==0){
          color.vec[j]<-"black"
        } 
      }
      
      #Plot the aligned signals:
      ymax.tot<-max(c(padded.reference.signal, padded.query.signal),na.rm=T) #For the full profiles
      ymin.tot<-min(c(padded.reference.signal, padded.query.signal),na.rm=T)
      ymax<-max(c(reference.line,query.line),na.rm=T) #For the lines
      ymin<-min(c(reference.line,query.line),na.rm=T)
      
      plot(1:length(padded.reference.signal),padded.reference.signal,typ="l",xlab="index",ylab="",
           xlim=c(1,length(padded.reference.signal)),ylim=c(ymin.tot,ymax.tot))
      par(new=T)
      plot(padded.reference.line.idxs,
           padded.reference.signal[padded.reference.line.idxs],col=color.vec,xlab="index",ylab="",
           xlim=c(1,length(padded.reference.signal)),ylim=c(ymin.tot,ymax.tot))
      abline(v=left.idx) #These mark off the line
      abline(v=right.idx)
      par(new=T)
      plot(1:length(padded.query.signal),padded.query.signal,typ="l",col="red",xlab="index",ylab="")      
      
      plot(left.idx:right.idx,reference.line,typ="l",ylim=c(ymin,ymax),
           main=paste("Line#:",i,"Comp val:",round(corr.val,3),round(dtwDist.val,3),round(dtwNDist.val,4),
                      round(mnscr,2),round(normscr,2)),
           xlab="",ylab="")
      par(new=T)
      plot(left.idx:right.idx,query.line,typ="l",col="red",ylim=c(ymin,ymax),xlab="index",ylab="")
      
      print(paste("Corr:      ",round(corr.val,4), "  ", "Corr score:",cscr))
      #print(paste("MI:        ",mi.val))
      print(paste("DTW dist:  ",round(dtwDist.val,4), "  ","DTW score: ",dscr))
      print(paste("N-DTW dist:",round(dtwNDist.val,4), "  ","NDTW score:",dnscr))
      #print(paste("Sum Sq:    ",ssqd.val))
      #print(paste("Sum Abs:   ",sabd.val))
      print(paste("Mean score:",mnscr))
      print(paste("Norm score:",normscr))
      
    }#End plot loop
    
  }
  
  #If parsed graphics device is still active, shut it off. It should be it plot==TRUE. Then plot the summary plot
  if(names(dev.cur()) != "null device") { 
    dev.off()
    
    ymax.tot<-max(c(padded.reference.signal, padded.query.signal),na.rm=T) #For the full profiles
    ymin.tot<-min(c(padded.reference.signal, padded.query.signal),na.rm=T)
    #print(ymin.tot)
    #Summary plot:
    #Scores (on scale 0-1)
    plot(mid.line.idxs, mean.scores,typ="h",lwd=4,col="blue",xlim=c(1,length(padded.reference.signal)),ylim=c(ymin.tot,ymax.tot),xlab="",ylab="")
    points(mid.line.idxs,norm.scores,pch=16,col="blue")
    par(new=T)
    plot(1:length(padded.reference.signal),padded.reference.signal,typ="l",xlab="index",ylab="",
         xlim=c(1,length(padded.reference.signal)),ylim=c(ymin.tot,ymax.tot))
    par(new=T)
    plot(padded.reference.line.idxs,
         padded.reference.signal[padded.reference.line.idxs],col=color.vec,xlab="index",ylab="",
         xlim=c(1,length(padded.reference.signal)),ylim=c(ymin.tot,ymax.tot))
    par(new=T)
    plot(1:length(padded.query.signal),padded.query.signal,typ="l",col="red",xlab="index",ylab="")
  }
  
  
  scores.mat<-cbind(corr.vals,dtwDist.vals,dtwNDist.vals,mean.scores,norm.scores,
                    dtw.lb_keogh.vals,dtw.lb_improved.vals,ssqd.vals,line.widths)
  #scores.mat<-cbind(corr.vals,dtwDist.vals,dtwNDist.vals,mean.scores,norm.scores)
  
  return(scores.mat)
  
}

#----------------------
#No plot
#LB_Keogh only for DTW
#----------------------
merge5a<-function(reference.signal, query.signal, reference.line.info.mat, score.tolerances, max.lag=10, rackQ=TRUE, renormalize.linesQ=FALSE, overlapped.linesQ=FALSE){
  
  #Align signals and determine shift to the line indices that will be needed:
  if(rackQ==TRUE){ #Rack the signals over each other with ccf to determine the "best" alignment
    aligned.signals<-align.mov.to.ref.pad.both2(reference.signal,query.signal,lagmax=max.lag,padtyp="NAs",printQ=FALSE)
  }
  if(rackQ==FALSE){ #Manually specify the alignment shift:
    aligned.signals<-manual.align.mov.to.ref.pad.both(reference.signal,query.signal,specified.lag=max.lag,padtyp="NAs",printQ=FALSE)
  }
  
  #Determine how the reference signal indices have shifted with any padding
  reference.idxs<-which(is.na(aligned.signals[[1]])==FALSE)
  
  #Loop over the line indices of the reference and pull out those index spans
  #from both the reference and the query.
  #First correct the line indices due to alignment/padding. Only needs to be done for the reference:
  orig.reference.line.idxs<-reference.line.info.mat[,1] #Line indices in the reference before alignment/padding
  aug.extrema.typ<-reference.line.info.mat[,2]   #For later in plots
  
  padded.reference.line.idxs<-rep(NA,length(orig.reference.line.idxs))
  
  line.widths <- NULL
  corr.vals<-NULL
  dtw.lb_keogh.vals<-NULL
  ssqd.vals<-NULL
  mean.scores<-NULL
  
  for(i in 1:length(orig.reference.line.idxs)){
    
    orig.line.idx <- orig.reference.line.idxs[i] #Take into account any padding in the reference
    padded.reference.line.idx <- reference.idxs[orig.line.idx]
    padded.reference.line.idxs[i] <- padded.reference.line.idx
    
  }
  
  padded.reference.signal<-aligned.signals[[1]]
  padded.query.signal<-aligned.signals[[2]]
  if(overlapped.linesQ==TRUE) {
    loop.idxs<-seq(1,length(padded.reference.line.idxs),1)  #Idxs for overlapping lines
  } else {
    loop.idxs<-seq(1,length(padded.reference.line.idxs),2) #Idxs for adjacent lines
  }
  
  for(i in 1:(length(loop.idxs))){
    
    left.idx<-padded.reference.line.idxs[loop.idxs[i]]    #Start of the line
    if(overlapped.linesQ==TRUE) {
      right.idx<-padded.reference.line.idxs[loop.idxs[i+2]] #End of the line for for overlapping lines
    } else {
      right.idx<-padded.reference.line.idxs[loop.idxs[i+1]] #End of the line for adjacent lines
    }
    
    #This is a jerry-rig not to kill the plot window at the last line:
    if(is.na(right.idx)){
      break()
    }
    
    reference.line <- padded.reference.signal[left.idx:right.idx] #Line from the reference
    query.line <- padded.query.signal[left.idx:right.idx]         #Line from the query
    
    #Re-normalize lines if desired:
    if(renormalize.linesQ==TRUE){
      reference.line <- norm.profile2(reference.line)
      query.line <- norm.profile2(query.line)
    }
    
    #Replace any NA padding with 0s. Necessary to compare them:
    reference.line[which(is.na(reference.line)==TRUE)] <- 0
    query.line[which(is.na(query.line)==TRUE)] <- 0
    
    #Line length (width):
    line.leng <- length(reference.line)
    line.widths <- c(line.widths, line.leng)
    
    #Correlation:
    if(sd(reference.line)==0 | sd(query.line)==0) {
#      corr.val<-NA
      corr.val<- -2
    } else {
      corr.val<-cor(reference.line,query.line)
    }
    
#     if(is.na(corr.val)){
#       corr.val<-(-2) #This should trip the flag not to accept this line and not mess up the scoring system below
#       #NOTE though, it is NOT a valid correlation value and should not be treated as such
#     }
    corr.vals<-c(corr.vals,corr.val)
    
    #Approx DTW:
    wind.leng <- floor(0.05*length(query.line)) #Window size is ~5% of query/reference length
    if(wind.leng >=5) {
      
      lemire.envelope <- computeEnvelope(query.line,wind.leng)
      query.Upper.lim <- lemire.envelope[,1]
      query.Lower.lim <- lemire.envelope[,2]
      dtw.lb_keogh.val <- LB_Keogh(reference.line,query.Upper.lim,query.Lower.lim) #Rcpp imp seems faster. Does lots of copying though      
      
    } else {
      #print("Skipping lower bound computation. Window too small.")
      dtw.lb_keogh.val <- 1e6     #Can't used NA for score checking later. Use big number instead
    }
    
    dtw.lb_keogh.vals <- c(dtw.lb_keogh.vals, dtw.lb_keogh.val)
    
    #Sum square deviations (Euclidian distances)
    ssqd.val<-sum((query.line-reference.line)^2)
    ssqd.vals<-c(ssqd.val,ssqd.vals)

    ##Sum abs deviations
    #sabd.val<-sum(abs(query.line-reference.line))
    #sabd.vals<-c(sabd.vals,sabd.val)    
    
    #Check scores:
    sc.ck.val<-sum(c(corr.val>=score.tolerances[1], 
                     dtw.lb_keogh.val<=score.tolerances[2], 
                     ssqd.val<=score.tolerances[3], 
                     line.leng>=score.tolerances[4]))

    if(sign(corr.val)==1 & sc.ck.val>=3){ #Corr must be positive and 3 out of 4 preliminary must be to tolerance
            
      #Correlation score
      cscr<-(corr.val) #Assumes only positive corr vals mad it through the check
      
      #Approx (relative) DTW score
      dscr<-(score.tolerances[2]-dtw.lb_keogh.val)/score.tolerances[2]
      if(dscr<0){ #If DTW dist exceeds the tolerance just set the score to 0 
        dscr<-0   #MAYBE CHANGE THIS TO DOUBLE PENALIZING THE SCORE INSTEAD??
      }
      
      #Euclidian (relative) score. Proxy for normalized DTW, which we can't do approximately
      dnscr<-(score.tolerances[3]-ssqd.val)/score.tolerances[3]
      if(dnscr<0){ #If  dist exceeds the tolerance just set the score to 0.
        dnscr<-0
      }
      
      mnscr<-mean(c(cscr,dscr,dnscr))
      
    } else {
      
      cscr<-0
      dscr<-0
      dnscr<-0
      mnscr<-0
      
    }
    mean.scores<-c(mean.scores,mnscr)    
  }
    
  scores.mat<-cbind(corr.vals,
                    dtw.lb_keogh.vals,
                    ssqd.vals,
                    line.widths,
                    mean.scores)
  
  return(scores.mat)
  
}