#library(dtw)
#library(bioDist)

#-----------------------------------------------------------------
#Segment a profile into "Lines" (Gives a possible answer to Mattia's retorical question: "What's a line?")
#Adapted from the function to barcode
#-----------------------------------------------------------------
find.line.idxs<-function(profil, tol, plotQ=FALSE, printQ=FALSE)
{

  #fit profile with splines and compute first derivative:
  #profile "function":
  psf<-splinefun(1:length(profil),profil)

  #First derivative:
  dpsf<-psf(1:length(profil),deriv=1)

  #Find zero crossings of derivative:
  deriv.zeros<-NULL
  deriv.vals<-NULL
  for(i in 1:(length(dpsf)-1))
  {
    if(dpsf[i]*dpsf[i+1]<0)
    {
      deriv.zeros<-c(deriv.zeros,i)
    }
  }

  #Drop out the ends if they had zero derivatives:
  if(1 %in% deriv.zeros)
  {
    deriv.zeros<-deriv.zeros[-which(deriv.zeros==1)]
  }
  if(length(profil) %in% deriv.zeros)
  {
    deriv.zeros<-deriv.zeros[-which(deriv.zeros==length(profil))]
  }

  #Pluck out any small extrema. They are probably noise from the
  #differencing procedure above, or due to a long flat set of
  #adjacent extrema points.
  extrema.height.diffs<-NULL
  for(i in 2:length(deriv.zeros))
  {
    #Subtract Right (current) extremum from extremum to the Left
    val<-profil[deriv.zeros[i]]-profil[deriv.zeros[i-1]]
    extrema.height.diffs<-c(extrema.height.diffs,val)
  }

  #Get rid of "noise" extrema:
  col1<-which(abs(extrema.height.diffs)>tol)

  #Keep the first (left most) extremum:
  col1<-col1+1 #Shift the indices by 1
  col1<-c(1,col1) #Tack on index for left most extremum

  #Extrema (df/dx~0) indices:
  deriv.zeros<-deriv.zeros[col1]

  #extrema.height.differences missing difference for first extremum.
  #To the left of first extremum may be a noise extremum (e.g.profile point 1).
  #To get a more definitive idea of whether it is a min or a max, subtract
  #it from the extremum to its right
  first.extremum.height.diff<-val<-profil[deriv.zeros[1]]-profil[deriv.zeros[2]]
  extrema.height.diffs<-c(first.extremum.height.diff,extrema.height.diffs)
  extrema.height.diffs<-extrema.height.diffs[col1]

  #Max=1,Min=-1. THERE SHOULD BE NO ZEROS!!!!!!!
  extrema.typ<-sign(extrema.height.diffs)
  if(0 %in% extrema.typ)
  {
    print("Undefined max/min!")
    print(cbind(col1,extrema.height.diffs,deriv.zeros,extrema.typ))
    return(0)
  }

  #Loop over all but the first and last extremum:
  if(printQ==TRUE){
    for(i in 2:(length(deriv.zeros)-1))
    {
      if(extrema.typ[i]!=extrema.typ[i+1])
      {
        print("PV or VP")
        print(paste("Left",deriv.zeros[i-1],"(",extrema.typ[i-1],")"))
        print(paste("Central",deriv.zeros[i],"(",extrema.typ[i],")"))
        print(paste("Right",deriv.zeros[i+1],"(",extrema.typ[i+1],")"))
        print("============================================================")
      }
      if(extrema.typ[i]==extrema.typ[i+1])
      {
        print("******PP or VV")
        print(paste("Left",deriv.zeros[i-1],"(",extrema.typ[i-1],")"))
        print(paste("Central",deriv.zeros[i],"(",extrema.typ[i],")"))
        print(paste("Right",deriv.zeros[i+1],"(",extrema.typ[i+1],")"))
        print("============================================================")
      }
    }
    #Now take care of the ends, ie left of first extremum and right of last extremum
    #Distance from beginning of profile to first extremum:
    start.to.first.span<-deriv.zeros[1]

    #Distance from last extremum to end of profile:
    last.to.end.span<-(length(profil)-deriv.zeros[length(deriv.zeros)])

    print(paste("Points from start to first extremum:",start.to.first.span))
    print(paste("Points from last extremum to end:",last.to.end.span))
  }

  #Tack the first and last profil idxs onto the "interesting" points
  aug.deriv.zeros<-c(1,deriv.zeros,length(profil))
  aug.extrema.typ<-c(0,extrema.typ,0) #Zeros for the endpoints

  if(plotQ==TRUE)
  {
    color.vec<-rep(NA,length(aug.extrema.typ))
    for(i in 1:length(color.vec)){
      if(aug.extrema.typ[i]==(-1)){
        color.vec[i]<-"blue"
      }
      if(aug.extrema.typ[i]==1){
        color.vec[i]<-"red"
      }
      if(aug.extrema.typ[i]==0){
        color.vec[i]<-"black"
      }
    }

    ymax<-max(profil)
    ymin<-min(profil)
    plot(1:length(profil),profil,typ="l",xlim=c(1,length(profil)),ylim=c(ymin,ymax))
    par(new=T)
    plot(aug.deriv.zeros,rep(0,length(aug.deriv.zeros)),col="green",xlim=c(1,length(profil)),ylim=c(ymin,ymax))
    par(new=T)
    plot(aug.deriv.zeros,profil[aug.deriv.zeros],col=color.vec,xlim=c(1,length(profil)),ylim=c(ymin,ymax))
  }

  line.idx.info.mat<-cbind(aug.deriv.zeros,aug.extrema.typ)

  return(line.idx.info.mat)

}

#-------------------------------------------------------------------------------------
#"Project" lines of a query profile onto reference profile (a chapter of the  dictionary)
#-------------------------------------------------------------------------------------
merge<-function(reference.signal, query.signal, reference.line.info.mat, score.tolerances, max.lag=10, rackQ=TRUE, renormalize.linesQ=FALSE, overlapped.linesQ=FALSE, plotQ=FALSE){

  #Align signals and determine shift to the line indices that will be needed:
  if(rackQ==TRUE){ #Rack the signals over each other with ccf to determine the "best" alignment
    aligned.signals<-align.mov.to.ref.pad.both(reference.signal,query.signal,lagmax=max.lag,padtyp="NAs",printQ=FALSE)
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
  corr.vals<-NULL
  #mi.vals<-NULL
  dtwDist.vals<-NULL
  dtwNDist.vals<-NULL
  #ssqd.vals<-NULL
  #sabd.vals<-NULL
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

    #Correlation:
    corr.val<-cor(reference.line,query.line)
    if(is.na(corr.val)){
      corr.val<-(-2) #This should trip the flag not to accept this line and not mess up the scoring system below
                    #NOTE though, it is NOT a valid correlation value and should not be treated as such
    }
    corr.vals<-c(corr.vals,corr.val)
    #DTW:
    dtw.info<-dtw(query.line,reference.line)
    dtwDist.val<-dtw.info$distance
    dtwDist.vals<-c(dtwDist.vals,dtwDist.val)
    dtwNDist.val<-dtw.info$normalizedDistance
    dtwNDist.vals<-c(dtwNDist.vals,dtwNDist.val)
    #Sum square deviations
    #ssqd.val<-sum((query.line-reference.line)^2)
    #ssqd.vals<-c(ssqd.val,ssqd.vals)
    ##Sum abs deviations
    #sabd.val<-sum(abs(query.line-reference.line))
    #sabd.vals<-c(sabd.vals,sabd.val)
    ##Mutual Info:
    #mi.val<-1-MIdist(rbind(reference.line,query.line))[[1]]
    #mi.vals<-c(mi.vals,mi.val)


    #Check scores:
    sc.ck.val<-sum(c(corr.val>=score.tolerances[1], dtwDist.val<=score.tolerances[2], dtwNDist.val<=score.tolerances[3]))
    #print(sc.ck.val)
    #print(sign(corr.val)==1)
    #print(sc.ck.val>=2)
    if(sign(corr.val)==1 & sc.ck.val>=2){ #Corr must be positive and 2 out of 3 distances must be to tolerance

      cscr<-(corr.val) #Assumes only positive corr vals mad it through the check

      dscr<-(score.tolerances[2]-dtwDist.val)/score.tolerances[2]
      if(dscr<0){ #If DTW dist exceeds the tolerance just set the score to 0
        dscr<-0   #MAYBE CHANGE THIS TO DOUBLE PENALIZING THE SCORE INSTEAD??
      }

      dnscr<-(score.tolerances[3]-dtwNDist.val)/score.tolerances[3]
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

  #scores.mat<-cbind(corr.vals,dtwDist.vals,dtwNDist.vals,ssqd.vals,sabd.vals)
  scores.mat<-cbind(corr.vals,dtwDist.vals,dtwNDist.vals,mean.scores,norm.scores)

  return(scores.mat)

}

#-----------------------------------------------------
#
#-----------------------------------------------------
plot.fv.on.reference<-function(reference.signal, query.signal, reference.line.info.mat, scores.mat, scores.typ=4, max.lag=10, rackQ=TRUE, renormalize.linesQ=FALSE, overlapped.linesQ=FALSE){

  #Align signals and determine shift to the line indices that will be needed:
  if(rackQ==TRUE){ #Rack the signals over each other with ccf to determine the "best" alignment
    aligned.signals<-align.mov.to.ref.pad.both(reference.signal,query.signal,lagmax=max.lag,padtyp="NAs",printQ=FALSE)
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

    #This is a jerry-rig not to kill the plot window at the last line:
    if(is.na(right.idx)){
      break()
    } else {
      mid.line.idx<-floor(mean(c(right.idx,left.idx))) #This idx will be used for plotting the summary plot below
      mid.line.idxs<-c(mid.line.idxs,mid.line.idx)
    }

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

  }

  #Plot the aligned signals:
  scores<-scores.mat[,scores.typ]
  ymax.tot<-max(c(padded.reference.signal, padded.query.signal),na.rm=T) #For the full profiles
  ymin.tot<-min(c(padded.reference.signal, padded.query.signal),na.rm=T)
  #ymax<-max(c(reference.line,query.line),na.rm=T) #For the lines
  #ymin<-min(c(reference.line,query.line),na.rm=T)

  plot(mid.line.idxs,scores,typ="h",lwd=4,col="blue",xlim=c(1,length(padded.reference.signal)),ylim=c(ymin.tot,ymax.tot),xlab="",ylab="")
  points(mid.line.idxs,scores,pch=16,col="blue")
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
