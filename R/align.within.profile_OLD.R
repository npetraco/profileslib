#################################################################
#Align profiles in a data frame within each group
#Return within group aligned profiles as a list
#Necessary because each group of profiles will have different lengths
#################################################################
align.within.profile.old<-function(dmat,lbls,lagmax,printQ=TRUE)
{
 num.groups<-nlevels(lbls)
 dat.aligned<-NULL
 
 for(i in 1:num.groups)
  {
   group.idxs<-which(lbls==i)
   dat.group<-dmat[group.idxs,]
   dat.group.aligned<-align.profile.group(dat.group,lagmax,printQ)
   print(paste("Done with group: ",i,sep=""))
   dat.aligned<-c(dat.aligned,list(dat.group.aligned))
  }
 
return(dat.aligned) 
}