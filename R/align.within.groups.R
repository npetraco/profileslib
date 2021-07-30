#################################################################
#Align profiles in a data frame within each group
#Return within group aligned profiles as a list
#Necessary because each group of profiles will have different lengths
#################################################################
align.within.groups<-function(dmat,lbls,ref.choice,lagmax,padtyp,printQ=TRUE,writeQ=FALSE,write.dir)
{
 num.groups<-nlevels(lbls)
 dat.aligned<-NULL
 
 for(i in 1:num.groups)
  {
   group.idxs<-which(lbls==i)
   dat.group<-dmat[group.idxs,]
   dat.group.aligned<-align.within.a.group(dat.group, ref.choice, lagmax, padtyp, printQ=FALSE)
   print(paste("Done with group: ",i,sep=""))
   #print(dim(dat.group))
   print(dim(dat.group.aligned))
   dat.aligned<-c(dat.aligned,list(dat.group.aligned))
  }

#NEED TO ADD WRITE FUNCTION BACK IN 
 
return(dat.aligned) 
}
