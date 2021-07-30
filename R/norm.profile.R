#################################################################
#Normalize z-heights of a profile
#IE-sets z to scale between 0 and 1
#CAREFUL. ASSUMES NAs are on the right hand side only!
#################################################################
norm.profile<-function(profile)
{

numNAs<-length(profile)-length(na.omit(profile))

nprofl<-na.omit(profile)
minp<-min(nprofl)
maxp<-max(nprofl)

nprofl<-apply(as.array(nprofl),1,function(x){(x-minp)/(maxp-minp)})
	
nprofl<-c(nprofl,rep(NA,numNAs))
return(nprofl)
 
}


#################################################################
#Normalize z-heights of a profile
#IE-sets z to scale between 0 and 1
#This version can handle NAs anywhere
#################################################################
norm.profile2<-function(profile)
{
  
  #numNAs<-length(profile)-length(na.omit(profile))
  na.idx.vec<-is.na(profile)
  
  nprofl<-na.omit(profile)
  minp<-min(nprofl,na.rm=T)
  maxp<-max(nprofl,na.rm=T)
  
  nprofl<-apply(as.array(nprofl),1,function(x){(x-minp)/(maxp-minp)})
  
  nprofl.corrected<-rep(NA,length(profile))
  count<-1
  #Put the normalized profile value back in the slot it started out in:
  for(i in 1:length(profile)){
    if(na.idx.vec[i]==FALSE){
      nprofl.corrected[i] <- nprofl[count]
      count <- count+1
    }
  }
    
  return(nprofl.corrected)
  
}


#################################################################
#Un-0/1 Normalize z-heights of a profile
#Handy for using with the simulator output.
#ASSUMES NO NAs!
#################################################################
unnorm.profile<-function(nprofl,maxp,minp)
{
  profl<-apply(as.array(nprofl),1,function(x){x*(maxp-minp)+minp})
  return(profl)
}


#-------------------------------------------------------------------
#Norm all the profiles making up a surface
#-------------------------------------------------------------------
norm.surface<-function(profiles)
{

num.profiles<-dim(profiles)[1]
normsurf<-t(apply(profiles,1,norm.profile))

return(normsurf)
	
}
