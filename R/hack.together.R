#################################################################
#Roughly put together a data frame by hacking off the ends of 
#groups so that they are all the same length 
#################################################################
hack.together<-function(dmat.list)
{
 minleng<-min(sapply(dmat.list,function(x){dim(x)[2]}))
 num.groups<-length(dmat.list)

 dmat<-NULL
 for(i in 1:num.groups)
  {
   group.dmat<-dmat.list[[i]][,1:minleng]
   dmat<-rbind(dmat,group.dmat)
  }

return(dmat)
 
}