#--------------------
#MRA filter profile
#--------------------
mra.filter.profile<-function(orig.profile, mra.obj, keep.levels, plotQ=FALSE) {
  
  filtered.signal<-numeric(length(orig.profile))
  for(i in 1:length(keep.levels)) {
    filtered.signal = filtered.signal + mra.obj[[keep.levels[i]]]    
  }
  #par( mfrow = c( 2, 1 ))
  if(plotQ==TRUE) {
    ymin<-min(orig.profile)
    ymax<-max(orig.profile)
    plot(orig.profile,typ="l",col="black",ylab="",xlab="",ylim=c(ymin,ymax))
    par(new=TRUE)
    plot(filtered.signal,typ="l",col="blue",ylab="",xlab="",ylim=c(ymin,ymax),lwd=3)
    
    print(paste("Mean original profile: ", mean(orig.profile)))
    print(paste("Mean filtered profile: ", mean(filtered.signal)))
    
  }
  
  return(filtered.signal)
    
}

mra.filter.routine<-function(profile.vec, basis, req.levl.depth, keep.levs, scaling.typ="None") {
  
  #Filter profile
  prof<-na.omit(profile.vec)
  prof.decomp<-wavMODWT(prof,  wavelet=basis, n.levels=req.levl.depth, keep.series=TRUE)
  prof.mra<-wavMRD(prof.decomp)
  filt.prof<-mra.filter.profile(prof, prof.mra, keep.levs, plotQ=FALSE)
  
  #Scale filtered profile if requested
  if(scaling.typ=="l2") {
    filt.prof<-l2.norm.profile(filt.prof)
  }
  if(scaling.typ=="zero-one") {
    filt.prof<-norm.profile(filt.prof)
  }
  if(scaling.typ=="auto scale") {
    filt.prof<-auto.scale.profile(filt.prof,typ="robust")
  }
  
  #NA pad filterd profile to the same length as the profile sent in
  na.pad<-(length(profile.vec)-length(filt.prof))
  filt.prof<-c(filt.prof,rep(NA,na.pad))
  
  return(filt.prof)
  
}