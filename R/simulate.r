#------------------------------------------------------------------------
#Simulate dwt output for a group of aligned signals.
#------------------------------------------------------------------------
simulate.group<-function(dmat,wbasis,split.level,num.sims) {
  
  num.prof<-nrow(dmat)
  num.pts<-ncol(dmat)
  max.level<-(floor((log(num.pts)/log(2))))
  
  #dwt the signals (uses waveslim) and re-arrange the coeficients into level matrices:
  transforms.info<-dwt.group(dmat, wbasis, max.level)
  na.idx.list<-transforms.info[[2]]
  grouped.level.mats<-group.level.coefs(transforms.info[[1]])
  
  #waveslim containers to hold dwt simulations
  coef.list.tmplate<-rep(list(NULL), max.level+1)
  sim.dwts.list<-rep(list(construct.dwt.obj(coef.list.tmplate, basis.choice=wbasis, boundary.typ="periodic")), num.sims)
  
  #Bulk simulate with one dist arcoss the whole level.
  #Simulate the level by drawing coefs across it (ie treat the level coefs as a random sample and sample from it)
  #This seems to be ok if the mags of the coefs are very small, and/or there are lots of them. Speeds up the 
  #simulation tremendously by doing this.
  if(split.level>0) {
    for(i in 1:split.level) {
      level.mat<-grouped.level.mats[[i]]
      
      #CLEAN LEVEL MAT HERE AND RE-STORE IN level.mat. Tones down edge coefs which can introduce spikes in the 
      #reconstructed signal.
      level.mat<-t(sapply(1:nrow(level.mat),function(x){clean.and.fill.level(level.mat[x,],typQ="mean",sig=3)}))
            
      #Initialize the matrix to hold simulated levels
      #sim.level.mat<-matrix(numeric(ncol(level.mat)*num.sims),nrow=num.sims,ncol=ncol(level.mat))
      
      for(j in 1:num.sims) {
        #Randomly select a row. This mixes in a little influence across signals.
        rand.row.num<-sample(1:nrow(level.mat),1)
        sim.level<-sample.coefs(level.mat[rand.row.num,],ncol(level.mat))
        #sim.level.mat[j,]<-sim.level
        #sim.dwts.list[[j]][[i]]<- sim.level.mat[j,]
        sim.dwts.list[[j]][[i]]<- sim.level
      }
    }
  }
  
  #Build a separate density fit for each shift. Slow for levels with lots of coefs (ie fine detail levels)
  if(split.level>0) {
    start.level<-split.level+1 #Start these right after the bulk level simulations
  } else {
    start.level <- 1           #If no levels were bulk simulated, start the simulated levels at level 1.
  }
  for(i in start.level:length(grouped.level.mats) ) {
    level.mat<-grouped.level.mats[[i]]
    sim.level.mat <- sapply(1:ncol(level.mat),function(x){simulate.column(level.mat[,x], num.sims)}) 
    for(j in 1:nrow(sim.level.mat)) {
      sim.dwts.list[[j]][[i]]<- sim.level.mat[j,]
    }
  }
  
  return(sim.dwts.list)
}

#---------------------------------------------------------------
#Wrapper to simulate dwts, transform them back to signals, 
#renormalize and chop off any tails due to padding
#---------------------------------------------------------------
simulate.signals<-function(dmat,wbasis,split.level,num.sims=1,clipQ=FALSE,num.clip, plotQ=FALSE) {
  
  #NOTE: Assumes no NAs at left of signal. All NAs are to the right.
  signals.length<-length(as.numeric(na.omit(dmat[1,])))

  #Zero indices for 0-padded signals. These can be to the left or the right of the signal.
  #IMPLEMENT LATER
  
  dwt.sims<-simulate.group(dmat,wbasis,split.level,num.sims)                 #sim dwts
  simulated.signals<-t(sapply(1:num.sims,function(x){ idwt(dwt.sims[[x]])})) #invert 
  simulated.signals<-simulated.signals[,1:signals.length]                    #Clipps off the tails from padding 
  
  #Clip the ends a bit
  if(clipQ==TRUE) {
    simulated.signals<-simulated.signals[,-c(1:num.clip[1],(signals.length-num.clip[2]):signals.length)]
  }
  
  #Re-normalize
  simulated.signals<-t(apply(simulated.signals,1,norm.profile))
  colnames(simulated.signals)<-NULL
  rownames(simulated.signals)<-NULL
  
  if(plotQ==TRUE) {
    plot(simulated.signals[1,],typ="l",col=1,ylim=c(0,1))
    if(num.sims>1) {
      for(i in 2:num.sims){
        par(new=T)
        plot(simulated.signals[i,],typ="l",col=i,ylim=c(0,1))
      }      
    }
  }
   
  return(simulated.signals)
    
}

#
#
#
refined.simulate.signals<-function(dmat,wbasis,split.level=0,num.sims=1,clipQ=FALSE,num.clip, 
                                   num.burn.in.iterations, num.burn.in.sims, num.keep.per.iter,
                                   reference.threshold, burn.in.threshold, keep.threshold, 
                                   plotQ=FALSE) {
  
  #First remove any NAs from input patterns
  dmat2<-t(sapply(1:nrow(dmat),function(x){as.numeric(na.omit(dmat[x,]))}))
  
  #Determine which input signals will be used in deciding which simulated signals to keep
  ref.signals<-NULL
  #corr.scrs<-NULL
  for(i in 1:nrow(dmat2)) {
    compto<-(1:nrow(dmat2))[-i]
    #print(compto)
    corr.vec<-sapply(compto,function(x){cor(dmat2[i,],dmat2[x,])})
    #corr.scrs<-rbind(corr.scrs,corr.vec)
    avg.cor<-100*mean(corr.vec)
    if(avg.cor>reference.threshold) {
      ref.signals<-c(ref.signals,i)
    }
  }
  print("Reference signals:")
  print(ref.signals)
  if(is.null(ref.signals)) {
    #print("BAD Correlation Scores????:")
    #print(corr.scrs)
    stop("BIG problem! All in-put signals correlate less than chosen threshold. Lower threshold or check signals!")
  }
  
  first.sims<-dmat2

  #burn-in simulations
  if(num.burn.in.iterations>0) {
    print("Burn-in")
    for(i in 1:num.burn.in.iterations) {
      first.sims<-simulate.signals(first.sims,wbasis,split.level,num.burn.in.sims,FALSE,c(0,0),FALSE)

      #Check which to keep
      compto<-1:num.burn.in.sims      
      cor.mat<-matrix(numeric(num.burn.in.sims*length(ref.signals)), nrow=num.burn.in.sims, ncol=length(ref.signals))
      for(j in 1:length(ref.signals)) {
        cors<-sapply(compto,function(x){cor(dmat2[ref.signals[j],],first.sims[x,])})
        cor.mat[,j]<-cors
      }
      keep.vec<-sapply(1:nrow(cor.mat),function(x){!(FALSE %in% ((100*cor.mat[x,])>=burn.in.threshold))})
      if(length(which(keep.vec==FALSE))) {
        first.sims<-first.sims[-which(keep.vec==FALSE),]
      }
      
      #Pick a few simulated signals at random for the next iteration
      first.sims<-first.sims[sample(1:nrow(first.sims),floor(0.1*nrow(first.sims)),replace=FALSE),]
      print(paste("Burn-in iter:",i,"done."))
    }
    burnin.sims<-first.sims
  }
  
  if(num.burn.in.iterations>0) {
    print(paste("Starting simumations with",nrow(dmat2),"real signals and",nrow(burnin.sims),"burn-in signals."))
    input.signals<-rbind(dmat2,burnin.sims)
  } else {
    print(paste("Starting simumations with",nrow(dmat2),"real signals."))
    input.signals<-dmat2
  }
  
  simulated.signals<-NULL
  count<-0
  iter<-1
  while(count<num.sims) {
    print(paste("Iteration:",iter))
    candidate.sims<-simulate.signals(input.signals,wbasis,split.level,num.sims,FALSE,c(0,0),FALSE)
    
    compto<-1:num.sims
    cor.mat<-matrix(numeric(num.sims*length(ref.signals)), nrow=num.sims, ncol=length(ref.signals))
    for(i in 1:length(ref.signals)) {
      cors<-sapply(compto,function(x){cor(dmat2[ref.signals[i],],candidate.sims[x,])})
      cor.mat[,i]<-cors
    }
    keep.vec<-sapply(1:nrow(cor.mat),function(x){!(FALSE %in% ((100*cor.mat[x,])>=keep.threshold))})
    if(length(which(keep.vec==FALSE))) {
      candidate.sims<-candidate.sims[-which(keep.vec==FALSE),]
    }
    
    simulated.signals<-rbind(simulated.signals,candidate.sims[1:num.keep.per.iter,])
    count<-count+num.keep.per.iter
    
    input.signals<-rbind(dmat2,simulated.signals) #Roger Xu suggestion. Thank him in paper.
    iter <- iter+1

  }
  
  #Clip the ends a bit
  if(clipQ==TRUE) {
    simulated.signals<-simulated.signals[,-c(1:num.clip[1],(signals.length-num.clip[2]):signals.length)]
  }
    
  if(plotQ==TRUE) {
    plot(simulated.signals[1,],typ="l",col=1,ylim=c(0,1))
    if(num.sims>1) {
      for(i in 2:num.sims){
        par(new=T)
        plot(simulated.signals[i,],typ="l",col=i,ylim=c(0,1))
      }      
    }
  }
  
  return(simulated.signals)
}

#------------------------------------------------------------------------
#Get dwts on a group of aligned signals
#NOTE: This will right-pad the signals if necessary (first with NAs then with 0s)
#------------------------------------------------------------------------
dwt.group<-function(dmat,wbasis,num.levels) {
  
  num.signals<-nrow(dmat)
  N<-ncol(dmat)
  #Is 2^num.levels an integer divisior of the signal lengths? Don't worry about NAs in dmat. They will be replaced later
  dmat2<-dmat
  #print(dim(dmat2))
  divQ <- (N/(2^num.levels) == floor(N/(2^num.levels)))
  if(divQ==FALSE) {
    ext.signal.leng<-(ceiling(N/(2^num.levels)) * 2^num.levels)
    pad.leng<-(ext.signal.leng - N)
    pad.mat<-matrix(rep(NA,num.signals*pad.leng), ncol=pad.leng, nrow=num.signals)
    dmat2<-cbind(dmat,pad.mat)
  }
  #print(dim(dmat2))
  
  #For each signal, make a note of which indices are NAs. Users can toss these later if desired
  na.idx.list<-rep(list(NULL),num.signals)
  for(i in 1:num.signals) {
    if(NA %in% dmat2[i,]) {
      signal.na.idxs<-which(is.na(dmat2[i,])==TRUE)
      na.idx.list[[i]]<-signal.na.idxs
    }
  }
  #print(na.idx.list)
  
  #Now replace NAs with 0s so dwt can be performed
  pad.typ<-list("zeros")
  dmat3<-replace.NAs(dmat2,padtype=pad.typ)
  dwts.list <- lapply(1:num.signals,function(x){dwt(dmat3[x,], n.levels=num.levels, wf=wbasis)})
  
  info.list<-list(dwts.list,na.idx.list)
  names(info.list)<-c("dwts","NA indices")
  
  return(info.list)
}

#--------------------------------------------------------------------------
#Rearrange level coeficents into matrices of each level across all signals
#--------------------------------------------------------------------------
group.level.coefs<-function(list.dwts) {
  
  num.signals<-length(list.dwts) #Below just using signal 1 as an example. All other signals should be the same info.
  max.level<-length(list.dwts[[1]]) #Note: Counts s too.
  level.lengths<-sapply(1:max.level, function(x){length(list.dwts[[1]][[x]]) })
  #print(max.level)
  
  #Initialize container to hold DWT coefs, organized by level
  coefs.mat.list<-lapply(1:max.level,function(x){matrix(rep(0,num.signals*level.lengths[x]),ncol=level.lengths[x],nrow=num.signals)})
  #print(sapply(1:max.level,function(x){dim(coefs.mat.list[[x]]) })  )
  
  #Collect the level coefs from the different signals together
  for(i in 1:num.signals) {
    signal.dwt.info<-list.dwts[[i]]
    for(j in 1:max.level) {
      level.coef.vector<-signal.dwt.info[[j]]
      coefs.mat.list[[j]][i,]<-level.coef.vector
    }
  }
  
  return(coefs.mat.list)
}

#------------------------------------------------------------------------
#Simulate a "column" (vector) of values, given the input values. Used to 
#simulate a column of dwt coefs at a level.
#------------------------------------------------------------------------
simulate.column<-function(cmat.col,samp.size) {
  
  #Test the column to see what to do with it
  zero.idxs<-which(cmat.col==0)
  if(length(zero.idxs)==length(cmat.col)) { #All input elements are zero. Just return a 0 vector
    sim.vec<-rep(0,samp.size)
  }
  if(length(zero.idxs)==(length(cmat.col)-1)) { #Only one non-zero input element. Sample from a fat tail dist.
    mt<-cmat.col[which(cmat.col!=0)]            #This may be a stupid idea
    sdt<-mt
    sim.vec<-rt(samp.size,4)
    sim.vec<-(sim.vec*sdt)+mt
  }
  if(length(zero.idxs)<=(length(cmat.col)-2)) { #Two or more non-zero input element.
   #
   if(length(zero.idxs)==0) {
     dfit<-density(cmat.col)                        #Just use default values for now  
   } else {
     dfit<-density(cmat.col[-zero.idxs])            #Just use default values for now
   }
   sim.vec<-sample(dfit$x,samp.size,replace=TRUE,prob=dfit$y) #Sample over the support returned by the fit function
  }
  
  return(sim.vec)
}

#--------------------------------------------------------------------
#Remove big spiked coefs from a level. These are probably edge coefs 
#This function is used to get a sample of "reasonable" coefs traversing
#a level or zero-ing out suspected spikes in prep for a traverse simulation. 
#It SHOULD NOT be used as a "denoising" tool
#--------------------------------------------------------------------
level.zero.out.outliers<-function(coef.vec,typQ="mean",sig) {
  
  zero.idxs<-which(coef.vec==0)
  non.zero.idxs<-((1:length(coef.vec))[-zero.idxs])
  cleaned.vec<-coef.vec[-zero.idxs]      #Pull out zeros to compute a better mean est. We put zeros in.
  if(typQ=="mean") {
    #Standardize the values to mean 0, sd 1.
    #Note, mean still does contain influence from outliers.
    std.cleaned.vec<-(cleaned.vec-mean(cleaned.vec))/sd(cleaned.vec)
  }
  if(typQ=="median") {
    #Try median-mad standardization instead of mean-sd bec there are probably big outliers, especially at the edges
    std.cleaned.vec<-(cleaned.vec-median(cleaned.vec))/mad(cleaned.vec)
  }
  
  too.pos.idxs<-which(std.cleaned.vec> sig )  
  too.neg.idxs<-which(std.cleaned.vec< -sig )  
  idxs.to.be.zeroed<-c(non.zero.idxs[too.neg.idxs],non.zero.idxs[too.pos.idxs])
  zeroed.coef.vec<-coef.vec
  zeroed.coef.vec[idxs.to.be.zeroed]<-0
  
  return(zeroed.coef.vec)
}

#--------------------------------------------------------------
#Sample based on a rough density over a fed in vector of coefs
#--------------------------------------------------------------
sample.coefs<-function(cvec,ssize) {
  dfit<-density(cvec)                                 #Just use default values for now
  samp<-sample(dfit$x,ssize,replace=TRUE,prob=dfit$y) #Sample over the support returned by the fit function
  return(samp)
}


#--------------------------------------------------------------------
#Clean and fill level with random coefs sampled from level's coef dist.
#--------------------------------------------------------------------
clean.and.fill.level<-function(coef.vec,typQ="mean",sig) {
  
  #Zero out designated outliers
  zeroed.coef.vec<-level.zero.out.outliers(coef.vec,typQ="mean",sig)
  
  #Pick out the zero indices. These will be replaced
  zero.idxs<-which(zeroed.coef.vec==0)
  
  #Use the non-zero coefs to build a density from which to sample
  dezeroed.coef.vec<-zeroed.coef.vec[-zero.idxs]
  samp.leng<-length(zero.idxs)
  replacement.coefs<-sample.coefs(dezeroed.coef.vec,samp.leng)
  
  #Replace the zero coefs with the sampled coefs
  filled.coef.vec<-zeroed.coef.vec
  for(i in 1:samp.leng) {
    filled.coef.vec[zero.idxs[i]] <- replacement.coefs[i]
  }
  
  return(filled.coef.vec)
  
}

#---------------------------------------------------------------
#Fill in "underhang" of signals with group grand mean chunks. 
#IE, where there are NAs in one signal but numbers
#in another, fill in the NAs with (possibly y-tranlated) 
#group grand mean values.
#---------------------------------------------------------------
initialize.signals<-function(dmat) {
  
  #Compute the group grand mean to fill in the gaps at the overhangs
  ggm<-colMeans(dmat,na.rm=T)
  ls.ggm<-which(is.na(ggm)==FALSE)[1]                          #left start ggm
  rs.ggm<-rev(1:ncol(dmat))[which(rev(is.na(ggm))==FALSE)[1]]  #right start ggm
  
  #loop over profiles
  #find points of each profile where they start/stop
  #fill in rest with y-translated ggm chunk
  dmat.new<-dmat #Initialize
  for(i in 1:nrow(dmat)) {
    signal<-dmat[i,]
    ls.sig<-which(is.na(signal)==FALSE)[1]                          #left start signal
    rs.sig<-rev(1:ncol(dmat))[which(rev(is.na(signal))==FALSE)[1]]  #right start signal
    
    #Check distance between ggm and signal at the start (left) of the signal
    ytrans.left<-(ggm[ls.sig]-signal[ls.sig])
    #print(ytrans.left)
    
    fill.left<-ggm[ls.ggm:(ls.sig-1)] #If no shift, this should stay the same
    
    #From the sign on the above, decide to shift the ggm chunk up or down
    if(sign(ytrans.left)==(-1)) { #ggm chunk is BELOW the signal at the attachment point. Shift the chunk UP.
      fill.left<-(ggm[ls.ggm:(ls.sig-1)]+abs(ytrans.left))
    }
    if(sign(ytrans.left)==(1)) { #ggm chunk is ABOVE the signal at the attachment point. Shift the chunk DOWN.
      fill.left<-(ggm[ls.ggm:(ls.sig-1)]-abs(ytrans.left))
    }
    #Fill the left gap:
    signal[ls.ggm:(ls.sig-1)]<-fill.left
 
    
    #Check distance between ggm and signal at the end (right) of the signal
    ytrans.right<-(ggm[rs.sig]-signal[rs.sig])
    
    fill.right<-ggm[rs.sig:rs.ggm] #If no shift, this should stay the same
    
    #From the sign on the above, decide to shift the ggm chunk up or down
    if(sign(ytrans.right)==(-1)) { #ggm chunk is BELOW the signal at the attachment point. Shift the chunk UP.
      fill.right<-(ggm[rs.sig:rs.ggm]+abs(ytrans.right))
    }
    if(sign(ytrans.right)==(1)) { #ggm chunk is ABOVE the signal at the attachment point. Shift the chunk DOWN.
      fill.right<-(ggm[rs.sig:rs.ggm]-abs(ytrans.right))
    }
    
    #Fill the right gap:
    signal[rs.sig:rs.ggm]<-fill.right
    
    #Replace the old row with the augmented row
    dmat.new[i,]<-signal
    
  }
  
  #Renormalize the signals:
  #****NOTE: This loads ALL the NAs (left and right of the profile) over to the right. That's ok, because at this point
  #they are all aligned AND the same length. IE, there are no NAs in one profile where another has numbers.
  #NOTE also that this process will destroy any between group alignment because NA padding on the left will be shifted 
  #to the right.
  dmat.new<-t(apply(dmat.new,1,norm.profile))
  colnames(dmat.new)<-NULL
  rownames(dmat.new)<-NULL
  
  return(dmat.new)
  
}