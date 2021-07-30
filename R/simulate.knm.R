#------------------------------------------------
#Peter and Martin's ideas to generate KNMs
#------------------------------------------------

#------------------------------------------------
#First feed in some groups of profiles. DWT each group 
#and collect the coefficients together. This preps
#the coefs for sampling. A group MUST
#have at least two observations at the moment.
#
#input.lbls is a vector of label names. It can have
#repeats and should be samples EXTERNAL to this 
#function
#
#Note that each different sampling of these coef
#matrices should give a new KNM group.
#------------------------------------------------
prep.dwt.groups.for.knm<-function(input.signals,input.lbls.selection,all.lbls,wbasis){
  
  #First determine witch group is the longest. That will define the max number of levels:
  maxJ<-(floor((log(max(sapply(1:nrow(input.signals),function(x){length(na.omit(input.signals[x,]))})))/log(2))))
  
  #Make a container (list) to hold the level matrices. It is maxJ+1 units long.
  #The first maxJ elements are for Detail levels. The last unit is for the Smooth level.
  #I think this should accomidate dwt's with differing number of levels, so long as it is only ~+/-1
  level.matrices<-rep(list(NULL),maxJ+1)
  names(level.matrices)<-c(paste("d",seq(1:maxJ),sep=""),paste("s",maxJ,sep=""))
  #print(level.matrices)

  #For each input group of signals, group dtw each group and collect level coefs together:
  for(i in 1:length(input.lbls.selection)){
    #Pick out a group:
    dmatg<-pick.out.groups(input.signals,all.lbls,c(input.lbls.selection[i]))[[1]]
    
    #Initialize it:
    dmatg.aug<-initialize.signals(dmatg) #Initialize (see simulate.r for details)
    
    #Remove any NAs from input patterns
    dmatg.aug2<-t(sapply(1:nrow(dmatg.aug),function(x){as.numeric(na.omit(dmatg.aug[x,]))}))
    
    #Max J for the group. May be less than maxJ but SHOULD NOT be greater:
    max.levg<-(floor((log(ncol(dmatg.aug2))/log(2))))
    if(max.levg>maxJ) {
      print(paste("Whoa! Level for group",input.lbls.selection[i],"exceeds maximum J computed for the inputsignals:",maxJ))
      print("Something's not right........")
      stop()
    }
    
    #Get the dwt coefs for the group of signals:
    group.dwt.container<-dwt.group(dmatg.aug2,wbasis=wbasis,num.levels=max.levg)[[1]]
    
    #Rearrange the coefs into matrices for each level:
    #Detail matrices are num.obs x 2^j
    #Smooth matrix is num.obs x 2^1 (or 2^0 if signal length is dyadic)
    level.mats.container<-group.level.coefs(group.dwt.container)
    
    #Store the level mats in the running level.matrices list:
    #They are stored side by side bec later we will randomly select columns.
    ##D
    for(j in 1:(max.levg)){
      colnames(level.mats.container[[j]])<-rep(input.lbls.selection[i],ncol(level.mats.container[[j]]))
      level.matrices[[j]] <- cbind(level.matrices[[j]], level.mats.container[[j]])
    }
    ##S
    level.matrices[[maxJ+1]] <- cbind(level.matrices[[maxJ+1]], level.mats.container[[max.levg+1]])
    print(paste("Group:",input.lbls.selection[i],"coefficients added."))
    #print(max.levg)
    #print(dim(dmatg.aug2))
    print("===========================================")    
  }
  
  return(level.matrices)
  
}


#------------------------------------------------
#Sample the column indices from the level matrices
#For each level, these will be the columns used to
#construct the KNM level matrices
#****Later add something to account for diadic length signal instead of the indicator variable????
#
#NOTE! THIS MAY SAMPLE CRAZY EDGE COEFS. WE NEED TO TEMPER THEM OR FIX THIS
#TO NOT SAMPLE EDGE COEFS!
#
#------------------------------------------------
sample.column.idxs<-function(lev.mats, dyadicQ=FALSE){
  
  col.widths<-sapply(1:length(lev.mats), function(x){ncol(lev.mats[[x]])})
  Jmax<-(length(lev.mats)-1)
  
  #Compute the number of dwt coefs needed for each level
  #CAUTION, below is a heuristic that may cause the code to BARF!
  if(dyadicQ==FALSE) {
    #For this option we went as high as we could in J
    lev.lengs<-sapply(Jmax:1,function(x){2^x})
  }
  if(dyadicQ==TRUE) {
    #For this option, we want the signals forming the group to be dyadic length. 
    lev.lengs<-sapply(Jmax:0,function(x){2^x})
  }  
  lev.lengs<-c(lev.lengs,lev.lengs[length(lev.lengs)]) #Added for the s "level"
  #print(lev.lengs)
  
  #Peter Zoon's idea:
  sampled.idxs<-rep(list(NULL), length(lev.mats))
  for(i in 1:length(col.widths)){
    lev.sampled.col.idxs<-sample(1:col.widths[i], lev.lengs[i])
    sampled.idxs[[i]]<-lev.sampled.col.idxs
  }
  
  return(sampled.idxs)
  
}

#------------------------------------------------
#Sample from columns of each level matrix output
#by prep.dwt.groups.for.knm EXTERNAL to this function.
#THEN send in sampled column indices. This way we
#can decide if we want to keep using the sampled
#columns
#------------------------------------------------
simulate.knm.group<-function(lev.mats,column.idx.list,wbasis,split.level,num.sims) {
  
  max.level<-(length(lev.mats)-1) 
  grouped.level.mats<-rep(list(NULL),max.level+1)
  for(i in 1:length(lev.mats)) {
    tmp<-lev.mats[[i]][,column.idx.list[[i]]]
    grouped.level.mats[[i]]<-tmp
  }
  #print(t(sapply(1:length(grouped.level.mats),function(x){dim(grouped.level.mats[[x]])})))
  
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
    
    #CLEAN LEVEL MAT HERE AND RE-STORE IN level.mat. Tones down edge coefs which can introduce spikes in the 
    #reconstructed signal.
    #level.mat<-t(sapply(1:nrow(level.mat),function(x){clean.and.fill.level(level.mat[x,],typQ="mean",sig=3)}))
    
    sim.level.mat <- sapply(1:ncol(level.mat),function(x){simulate.column(level.mat[,x], num.sims)}) 
    for(j in 1:nrow(sim.level.mat)) {
      sim.dwts.list[[j]][[i]]<- sim.level.mat[j,]
    }
  }
  
  return(sim.dwts.list)
}


#------------------------------------------------
#Clean up and back-transform the simulated dwt
#into a rough set of simulated signals.
#Renormalize and chop off any tails due to padding
#---------------------------------------------------------------
simulate.knm.signals<-function(dwt.level.mats, sampled.idxs, wbasis, split.level,num.sims=1,clipQ=FALSE,num.clip, plotQ=FALSE) {
    
  #dwt.sims<-simulate.group(dmat,wbasis,split.level,num.sims)                #sim dwts
  dwt.sims<-simulate.knm.group(dwt.level.mats, sampled.idxs, wbasis, split.level, num.sims)
  simulated.signals<-t(sapply(1:num.sims,function(x){ idwt(dwt.sims[[x]])})) #invert 
  #simulated.signals<-simulated.signals[,1:signals.length]                    #Clipps off the tails from padding 
  
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

#------------------------------------------------
#
#------------------------------------------------
refined.simulate.knm.signals<-function(dwt.level.mats, sampled.idxs,wbasis,split.level=0,num.sims=1,                                    
                                   num.references, reference.lengths, reference.threshold,
                                   num.burn.in.iterations, num.burn.in.sims, burn.in.threshold,
                                   num.keep.per.iter, keep.threshold, 
                                   plotQ=FALSE) {
  
  #Generate some reference signals using the rough knm simulator:
  ref.idxs<-sample(1:num.sims,num.references,replace=FALSE)
  refs<-simulate.knm.signals(dwt.level.mats, sampled.idxs, wbasis, split.level,num.sims,clipQ=FALSE,num.clip, plotQ=FALSE)[ref.idxs,1:reference.lengths]
  
  print("Initiallizing Signals Generated>")
  
  knm.sims<-refined.simulate.signals(refs,wbasis,split.level,num.sims,clipQ=FALSE,num.clip=NULL, 
                                     num.burn.in.iterations, num.burn.in.sims, num.keep.per.iter,
                                     reference.threshold, burn.in.threshold, keep.threshold,plotQ)
  
  return(knm.sims)

}