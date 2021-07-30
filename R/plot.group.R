#------------------------------------------------------
#Plot all the profiles in the group, reference first
#------------------------------------------------------
plot.group<-function(dmat,ref=NULL,ret) {
  num.profiles<-dim(dmat)[1]
  num.pts<-dim(dmat)[2]
  
  prof.lengs<-sapply(1:num.profiles,function(x){length(na.omit(dmat[x,]))})
  if(is.null(ref)){
    ref.idx<-which(prof.lengs==max(prof.lengs))[1] #Grab the first if many have the same (max) length
  } else {
    ref.idx<-ref
  }
  
  prof.idxs<-(1:num.profiles)[-ref.idx]
  
  #Plot the reference first:
  par(new=FALSE)
  plot(1:num.pts,dmat[ref.idx,],col=ref.idx,typ="l")
  par(new=T)
  print(paste("Plotted: ", ref.idx," (reference)",sep=""))
  Sys.sleep(ret)
  
  #Plot the other profiles:
  for(i in 1:length(prof.idxs)) {
    plot(1:num.pts,dmat[prof.idxs[i],],col=prof.idxs[i],typ="l")
    par(new=T)
    #print(prof.idxs[i])
    print(paste("Plotted: ", prof.idxs[i],sep=""))
    Sys.sleep(ret)
    
  }
  par(new=FALSE)
    
}


#------------------------------------------------------
#Plot all the profiles in the group, reference first.
#Scale the y axis to be the same for all plots
#------------------------------------------------------
plot.group2<-function(dmat,ref=NULL,ret,ylims) {
  num.profiles<-dim(dmat)[1]
  num.pts<-dim(dmat)[2]
  
  prof.lengs<-sapply(1:num.profiles,function(x){length(na.omit(dmat[x,]))})
  if(is.null(ref)){
    ref.idx<-which(prof.lengs==max(prof.lengs))[1] #Grab the first if many have the same (max) length
  } else {
    ref.idx<-ref
  }
  
  prof.idxs<-(1:num.profiles)[-ref.idx]
  
  #Plot the reference first:
  par(new=FALSE)
  plot(1:num.pts,dmat[ref.idx,],col=ref.idx,typ="l",ylim=ylims)
  par(new=T)
  print(paste("Plotted: ", ref.idx," (reference)",sep=""))
  Sys.sleep(ret)
  
  #Plot the other profiles:
  for(i in 1:length(prof.idxs)) {
    plot(1:num.pts,dmat[prof.idxs[i],],col=prof.idxs[i],typ="l",ylim=ylims)
    par(new=T)
    #print(prof.idxs[i])
    print(paste("Plotted: ", prof.idxs[i],sep=""))
    Sys.sleep(ret)
    
  }
  par(new=FALSE)
  
}


#----------------------------------------------------------
#Plot the refenence profile vs all the rest, sequentially
#----------------------------------------------------------
plot.ref.vs.others<-function(dmat,ref=NULL,ret) {
  num.profiles<-dim(dmat)[1]
  num.pts<-dim(dmat)[2]
  
  prof.lengs<-sapply(1:num.profiles,function(x){length(na.omit(dmat[x,]))})
  if(is.null(ref)){
    ref.idx<-which(prof.lengs==max(prof.lengs))[1] #Grab the first if many have the same (max) length 
  } else {
    ref.idx<-ref
  }
  
  prof.idxs<-(1:num.profiles)[-ref.idx]
  
  par(new=FALSE)  
  #Plot the other profiles:
  for(i in 1:length(prof.idxs)) {
    plot(1:num.pts,dmat[ref.idx,],col=ref.idx,typ="l")
    par(new=T)
    plot(1:num.pts,dmat[prof.idxs[i],],col=prof.idxs[i],typ="l", main=paste("Ref=",ref.idx,"vs. profile",prof.idxs[i]))
    par(new=F)
    #print(prof.idxs[i])
    Sys.sleep(ret)
    
  }
  par(new=FALSE)
  
}