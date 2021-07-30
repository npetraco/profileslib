#################################################################
#"Best" Align BETWEEN groups of profile vectors.
#"Best" is cross-correlation between grand mean profiles for 
#each group.
#
#Each group is then shifted en-mass by the parameters found
#by aligning the group grand means.
#
#In the code below:
#ref is the stationary (longest) group grand mean profile. 
#mov is one of the other (shorter) group grand mean profiles 
#that is sfhited
#
#################################################################
align.between.profile<-function(within.aligned.dmat.list,lbls,lagmax,printQ=FALSE)
{

#Compute a grand mean for each aligned group of profiles
grp.means<-lapply(within.aligned.dmat.list,colMeans)

#Mean vectors are all generally different in length
#Even them out by tacking on NA vectors
grp.means.lengs<-unlist(lapply(grp.means,length))
max.leng<-max(grp.means.lengs)
leng.na.vecs<-(max.leng-grp.means.lengs)
grp.means.mat<-NULL
for(i in 1:length(leng.na.vecs))
 {
  na.vec<-rep(NA,leng.na.vecs[i])
  aug.grp.mean.vec<-c(grp.means[[i]],na.vec)
  grp.means.mat<-rbind(grp.means.mat,aug.grp.mean.vec)
 }

#Borrowing from align.group function.
#Find "best" alignments between group grand means.
#Use these shifts later to shift the group BLOCKS the same amount as the group grand mean was shifted.

#Determine an "anchor" group which will serve as the refenence with which to align the
#other groups 
refp.idx<-which(grp.means.lengs==max(grp.means.lengs))[1] #pick the first if there are many "longest"
if(printQ==TRUE)
 {
  print("***************************************")
  print(paste("Longest group grand mean profile is #",refp.idx,sep=""))
  print("It is the reference")
  print("***************************************")
 }

grp.names<-levels(lbls)
reordered.grp.names<-c(grp.names[-refp.idx],grp.names[refp.idx])

#Count the number of replicates per group.
#Drop the number of replicates in the reference group.
num.reps<-sapply(levels(lbls),function(x){sum(lbls==x)})
num.reps.reduced<-num.reps[-refp.idx] #Needed to replicate the shift vectors for the moving groups


#Determine between group alignment by registering group grand means:
refp<-grp.means.mat[refp.idx,]
sftmat<-NULL
for(i in 1:nrow(grp.means.mat))
 {
  if(i!=refp.idx)
   {
    croscor<-ccf(refp,grp.means.mat[i,],lag.max=lagmax,plot=FALSE,na.action=na.omit) #Using Cross-Correlation for alignment
    lags<-croscor$lag
    max.corr.idx<-which(croscor$acf[,,1]==max(croscor$acf[,,1]))
    maxlag<-lags[max.corr.idx]
    if(printQ==TRUE)
     {
      print(paste("Group grand mean profiles #",i," vs. #",refp.idx,"(ref.)",sep=""))     
      print(paste("Index of max corr: ",max.corr.idx,sep=""))
      print(paste("Max corr at lag: ",maxlag,sep=""))
      print(paste("**********Max Corr: ",max(croscor$acf[,,1]),sep=""))
     }

    lag.na.vec<-rep(NA,abs(maxlag))

    if(printQ==TRUE)
     {
      print(paste("Length of Reference: ",length(refp),sep=""))
      print(paste("Length of Profile: ",length(grp.means.mat[i,]),sep=""))
      print("XX NOW REMOVE NAs XX")
     }
         
    refp<-na.omit(refp)
    movp<-na.omit(grp.means.mat[i,])

    if(printQ==TRUE)
     {
      print(paste("Length of Reference NOW: ",length(refp),sep=""))
      print(paste("Length of Profile NOW: ",length(movp),sep=""))
     }
    
    dif<-length(refp)-length(movp)

    if(printQ==TRUE)
     {
      print(paste("Length difference: ",dif,sep=""))
      print("================================")
     }
     
    #shifts: prepend ref, append ref, prepend mov, append mov
    sftvec<-c(0,0,0,0)        
    if(sign(maxlag)==1) #Pushing mov profile forward
     {
      sftvec[3]<-maxlag #prepend lag to mov profile
      if(sign(dif-maxlag)==1)
       {
       	sftvec[4]<-(dif-maxlag) #if mov profile does not go past ref, append dif-lag to mov profile 
       }
      if(sign(dif-maxlag)==-1)
       {
       	sftvec[2]<-(maxlag-dif) #if mov profile goes past ref, append lag-dif to ref profile
       }
     }
     
    if(sign(maxlag)==-1) #Lagging mov profile backward
     {
      sftvec[1]<-abs(maxlag) #prepend lag to ref
      sftvec[4]<-(abs(maxlag)+dif)	#append lag+dif to mov
     }

    sftvec<-c(length(movp),sign(maxlag),sign(dif-maxlag),sftvec)        
    sftmat<-rbind(sftmat,sftvec)
    
   }#end if i!=refp.idx	 	
 }#end for

if(printQ==TRUE)
 {
  colnames(sftmat)<-c("len mov|","sgn lag|","sgn dif-lag|","ref pre|","ref app|","mov pre|","mov app")
  print(sftmat)
 }

#Build a prepend and append block of NAs for each of the shifting groups
all.grp.names<-as.numeric(levels(lbls))
mov.grp.names<-all.grp.names[-refp.idx]

mp<-max(sftmat[,4]) #Max prepend length
ma<-max(sftmat[,5]) #Max append length
refp.leng<-length(refp)

#Drop the reference group from the data to be aligned. It will be tacked on
#to the results at the end of the process. The data will then be rearranged
#into its original order.
within.aligned.dmat.dropped.ref.list<-within.aligned.dmat.list[-refp.idx]
#print(refp.idx)
#print(length(within.aligned.dmat.list))
#print(length(within.aligned.dmat.dropped.ref.list))

#Make an empty list to hold the aligned blocks:
num.grps<-length(levels(lbls))
shifted.profiles.list<-rep(list(NULL),num.grps)

#Align the moving groups according to the parameters in sftmat:
for(i in 1:nrow(sftmat))
 {
  rows.na<-as.numeric(num.reps.reduced[i])
  pre.cols.na<-sftmat[i,6]
  app.cols.na<-sftmat[i,7]
  
  shifted.grp.block<-within.aligned.dmat.dropped.ref.list[[i]]
  #print(i)
  #print(dim(shifted.grp.block))
  
  #Case 1:
  #Lag backward, or Lag forward and get overhang
  if(sftmat[i,2]<0 || (sftmat[i,2]>0 && sftmat[i,3]>0))
   {
   	#print("*CASE1")
   	#Prepend 1. Prepend chunk needed to account for shift wrt the reference.
    if(pre.cols.na>0)
     {
      nas<-rep(NA,rows.na*pre.cols.na)
      grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
      #print(dim(grp.pre.block))
      #print(dim(shifted.grp.block))
      shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
     }
    #Append 1. Append chunck for same thing as above.
    if(app.cols.na>0)
     {
      nas<-rep(NA,rows.na*app.cols.na)
      grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
      #print(rows.na)
      #print(dim(grp.app.block))
      #print(dim(shifted.grp.block))      
      shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
     }
    #Prepend 2. Extra prepend chunk need to account for shifting of the other blocks.
    ex.na.cols<-(refp.leng + mp - ncol(shifted.grp.block))
    if(ex.na.cols>0)
     {
      nas<-rep(NA,rows.na*ex.na.cols)
      ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
      shifted.grp.block<-cbind(ex.grp.pre.block,shifted.grp.block)
     }
    #Append 2. Extra append chunck for same thing as above.
    if(ma>0)
     {
      nas<-rep(NA,rows.na*ma)
      ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ma)
      shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
     }
   	shifted.profiles.list[[i]]<-shifted.grp.block
   }
   
   #Case 2:
   #Lag forward but get no overhang
   if(sftmat[i,2]>0 && sftmat[i,3]<0)
    {
     #print("*CASE2")
     #Prepend 1
     if(pre.cols.na>0)
      {
       nas<-rep(NA,rows.na*pre.cols.na)
       grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
       shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
      }
     #Append 1
     if(app.cols.na>0)
      {
       nas<-rep(NA,rows.na*app.cols.na)
       grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
       shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
      }
     #Prepend 2
     if(mp>0)
      {
       nas<-rep(NA,rows.na*mp)
       ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=mp)
       shifted.grp.block<-cbind(shifted.grp.block,ex.grp.pre.block)
      }
     #Append 2
     ex.na.cols<-(refp.leng + ma - ncol(shifted.grp.block))
     if(ex.na.cols>0)
      {
       nas<-rep(NA,rows.na*ex.na.cols)
       ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
       shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
      }
     shifted.profiles.list[[i]]<-shifted.grp.block
    }
    
   #Case 3:
   #No shifting with respect to the reference required: 
   if(sftmat[i,2]==0) 
    {
     #print(paste("Group:",i))	
     #print("*CASE3")
     #Prepend 1
     if(pre.cols.na>0)
      {
       nas<-rep(NA,rows.na*pre.cols.na)
       grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
       shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
      }
     #Append 1
     if(app.cols.na>0)
      {
       nas<-rep(NA,rows.na*app.cols.na)
       grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
       shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
      }
     #Prepend 2
     if(mp>0)
      {
       #print(mp)
       nas<-rep(NA,rows.na*mp)
       ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=mp)
       #print(dim(ex.grp.pre.block))
       #print(dim(shifted.grp.block))
       shifted.grp.block<-cbind(shifted.grp.block,ex.grp.pre.block)
       #print(dim(shifted.grp.block))
      }
     #Append 2                     Tacked on lagmax to sum. Fixed problem with screwdrivers but DONT KNOW IF THIS IS ALWAYS OK!!!!!!! 
     ex.na.cols<-(refp.leng + ma + lagmax - ncol(shifted.grp.block))
     #print(ma)
     #print(refp.leng)
     #print(ncol(shifted.grp.block))
     #print(ex.na.cols)
     if(ex.na.cols>0)
      {
       nas<-rep(NA,rows.na*ex.na.cols)
       ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
       shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
      }
     shifted.profiles.list[[i]]<-shifted.grp.block      
    }
  
 }

#Finally pad the reference block with NA mats:
refp.block<-within.aligned.dmat.list[[refp.idx]]
refp.block.before<-within.aligned.dmat.list[[refp.idx]]

if(mp>0)
 {
  prepend.na.block<-t(replicate(nrow(refp.block),rep(NA,mp)))
  refp.block<-cbind(prepend.na.block,refp.block) 
 }
if(ma>0)
 {
  append.na.block<-t(replicate(nrow(refp.block),rep(NA,ma)))
  refp.block<-cbind(refp.block,append.na.block) 
 }

shifted.profiles.list[[num.grps]]<-refp.block

#Rearrange the group blocks into their original order:
reordered.grp.idxs<-as.numeric(reordered.grp.names)
#print(reordered.grp.idxs)
shifted.profiles<-NULL
for(i in 1:length(shifted.profiles.list))
 {
  perm.idx<-which(reordered.grp.idxs==i)
  #print(paste("Group:",perm.idx, dim(shifted.profiles.list[[perm.idx]]) ))
  #print(dim(shifted.profiles.list[[perm.idx]]))
  shifted.profiles<-rbind(shifted.profiles,shifted.profiles.list[[perm.idx]])
 }

na.removed.shifted.profiles<-t(na.omit(data.frame(t(shifted.profiles))))

return(list(shifted.profiles,na.removed.shifted.profiles))
	
}


#################################################################
#Align an unknown of arbitrary length to a specific group of
#knowns.
#Useful for repeating a model fit for different choices of group
#to align unknown to.
#################################################################
align.unknown.to.group<-function(unk.obs,specific.grp,X.within.groups.aligned.list,lagmax,printQ=FALSE,plotQ=FALSE)
{

#grp.unk.aligned.to<-specific.grp
print("==============================================================================")
print(paste("Aligning UNKNOWN to REFERENCE GROUP #", specific.grp))

#**aligning unknown on this group:
ref.grp<-X.within.groups.aligned.list[[specific.grp]]

#align unknown to reference group's mean:
ref.grp.mean<-colMeans(ref.grp)

#align unknown here:
unk.obs.mov.list<-align.mov.to.ref(ref.grp.mean,unk.obs,lagmax,padtyp=list("NAs"),printQ)
unk.obs.mov<-unk.obs.mov.list[[1]]

#Chop reference group to be the same length as the aligned unknown:
ref.grp.chop.left<-unk.obs.mov.list[[2]][3]
ref.grp.chop.right<-unk.obs.mov.list[[2]][5]
print(paste("CHOP reference group LEFT:",ref.grp.chop.left,"CHOP reference group RIGHT:", ref.grp.chop.right))

if(ref.grp.chop.right>0)
 {
  idxs.drop.right<-sort((dim(ref.grp)[2]):((dim(ref.grp)[2])-(ref.grp.chop.right-1)))
  ref.grp.chop<-ref.grp[,-idxs.drop.right]         #chop ref on right
  unk.obs.mov<-unk.obs.mov[-idxs.drop.right]       #remove padding
 }
if(ref.grp.chop.left>0)
 {
  idxs.drop.left<-(1:ref.grp.chop.left)
  ref.grp.chop<-ref.grp[,-idxs.drop.left]         #chop ref on left
  unk.obs.mov<-unk.obs.mov[-idxs.drop.left]       #remove padding
 }
if((ref.grp.chop.left==0) & (ref.grp.chop.right==0))
 {
  ref.grp.chop<-ref.grp
 }

#Show alignment if requested:
if(plotQ==TRUE)
 {
  for(j in 1:dim(ref.grp.chop)[1])
   {
    plot(1:length(ref.grp.chop[j,]),ref.grp.chop[j,],typ="l",col=j,ylim=c(0,1))
    par(new=TRUE)
   }  
  plot(1:length(unk.obs.mov),unk.obs.mov,col="blue",ylim=c(0,1))
 }

#print(length(X.within.groups.aligned.list))
X.within.groups.aligned.list.MOD<-X.within.groups.aligned.list
#print(length(X.within.groups.aligned.list.MOD))
#X.within.groups.aligned.list.MOD[[specific.grp]]<-NULL
X.within.groups.aligned.list.MOD[[specific.grp]]<-ref.grp.chop
#print(length(X.within.groups.aligned.list.MOD))

return(list(unk.obs.mov, X.within.groups.aligned.list.MOD))

}

#####################################
#Another attempt at this
#####################################
# align.between.groups.with.unknown.in.a.group2<-function(grp.idx,within.aligned.dmat.list,lbls,lagmax,printQ=FALSE)
# {
# 
# within.aligned.dmat.list.tmp<-within.aligned.dmat.list
# grp.idx.mat<-within.aligned.dmat.list.tmp[[grp.idx]]
# grp.idx.mat<-grp.idx.mat[-1,]
# within.aligned.dmat.list.tmp[[grp.idx]]<-grp.idx.mat
#   
# #Compute a grand mean for each aligned group of profiles
# grp.means<-lapply(within.aligned.dmat.list.tmp,colMeans)
# 
# #Mean vectors are all generally different in length
# #Even them out by tacking on NA vectors
# grp.means.lengs<-unlist(lapply(grp.means,length))
# max.leng<-max(grp.means.lengs)
# leng.na.vecs<-(max.leng-grp.means.lengs)
# grp.means.mat<-NULL
# for(i in 1:length(leng.na.vecs))
#  {
#   na.vec<-rep(NA,leng.na.vecs[i])
#   aug.grp.mean.vec<-c(grp.means[[i]],na.vec)
#   grp.means.mat<-rbind(grp.means.mat,aug.grp.mean.vec)
#  }
# 
# #Borrowing from align.group function.
# #Find "best" alignments between group grand means.
# #Use these shifts later to shift the group BLOCKS the same amount as the group grand mean was shifted.
# 
# #Determine an "anchor" group which will serve as the refenence with which to align the
# #other groups 
# refp.idx<-which(grp.means.lengs==max(grp.means.lengs))[1] #pick the first if there are many "longest"
# if(printQ==TRUE)
#  {
#   print("***************************************")
#   print(paste("Longest group grand mean profile is #",refp.idx,sep=""))
#   print("It is the reference")
#   print("***************************************")
#  }
# 
# grp.names<-levels(lbls)
# reordered.grp.names<-c(grp.names[-refp.idx],grp.names[refp.idx])
# 
# #Count the number of replicates per group.
# #Drop the number of replicates in the reference group.
# num.reps<-sapply(levels(lbls),function(x){sum(lbls==x)})
# num.reps.reduced<-num.reps[-refp.idx] #Needed to replicate the shift vectors for the moving groups
# 
# 
# #Determine between group alignment by registering group grand means:
# refp<-grp.means.mat[refp.idx,]
# sftmat<-NULL
# for(i in 1:nrow(grp.means.mat))
#  {
#   if(i!=refp.idx)
#    {
#     croscor<-ccf(refp,grp.means.mat[i,],lag.max=lagmax,plot=FALSE,na.action=na.omit) #Using Cross-Correlation for alignment
#     lags<-croscor$lag
#     max.corr.idx<-which(croscor$acf[,,1]==max(croscor$acf[,,1]))
#     maxlag<-lags[max.corr.idx]
#     if(printQ==TRUE)
#      {
#       print(paste("Group grand mean profiles #",i," vs. #",refp.idx,"(ref.)",sep=""))     
#       print(paste("Index of max corr: ",max.corr.idx,sep=""))
#       print(paste("Max corr at lag: ",maxlag,sep=""))
#       print(paste("**********Max Corr: ",max(croscor$acf[,,1]),sep=""))
#      }
# 
#     lag.na.vec<-rep(NA,abs(maxlag))
# 
#     if(printQ==TRUE)
#      {
#       print(paste("Length of Reference: ",length(refp),sep=""))
#       print(paste("Length of Profile: ",length(grp.means.mat[i,]),sep=""))
#       print("XX NOW REMOVE NAs XX")
#      }
#          
#     refp<-na.omit(refp)
#     movp<-na.omit(grp.means.mat[i,])
# 
#     if(printQ==TRUE)
#      {
#       print(paste("Length of Reference NOW: ",length(refp),sep=""))
#       print(paste("Length of Profile NOW: ",length(movp),sep=""))
#      }
#     
#     dif<-length(refp)-length(movp)
# 
#     if(printQ==TRUE)
#      {
#       print(paste("Length difference: ",dif,sep=""))
#       print("================================")
#      }
#      
#     #shifts: prepend ref, append ref, prepend mov, append mov
#     sftvec<-c(0,0,0,0)        
#     if(sign(maxlag)==1) #Pushing mov profile forward
#      {
#       sftvec[3]<-maxlag #prepend lag to mov profile
#       if(sign(dif-maxlag)==1)
#        {
#          sftvec[4]<-(dif-maxlag) #if mov profile does not go past ref, append dif-lag to mov profile 
#        }
#       if(sign(dif-maxlag)==-1)
#        {
#        	sftvec[2]<-(maxlag-dif) #if mov profile goes past ref, append lag-dif to ref profile
#        }
#      }
#      
#     if(sign(maxlag)==-1) #Lagging mov profile backward
#      {
#       sftvec[1]<-abs(maxlag) #prepend lag to ref
#       sftvec[4]<-(abs(maxlag)+dif)	#append lag+dif to mov
#      }
# 
#     sftvec<-c(length(movp),sign(maxlag),sign(dif-maxlag),sftvec)        
#     sftmat<-rbind(sftmat,sftvec)
#     
#    }#end if i!=refp.idx	 	
#  }#end for
# 
# if(printQ==TRUE)
#  {
#   colnames(sftmat)<-c("len mov|","sgn lag|","sgn dif-lag|","ref pre|","ref app|","mov pre|","mov app")
#   print(sftmat)
#  }
# 
# #Build a prepend and append block of NAs for each of the shifting groups
# all.grp.names<-as.numeric(levels(lbls))
# mov.grp.names<-all.grp.names[-refp.idx]
# 
# mp<-max(sftmat[,4]) #Max prepend length
# ma<-max(sftmat[,5]) #Max append length
# refp.leng<-length(refp)
# 
# #Drop the reference group from the data to be aligned. It will be tacked on
# #to the results at the end of the process. The data will then be rearranged
# #into its original order.
# within.aligned.dmat.dropped.ref.list<-within.aligned.dmat.list[-refp.idx]
# #print(refp.idx)
# #print(length(within.aligned.dmat.list))
# #print(length(within.aligned.dmat.dropped.ref.list))
# 
# #Make an empty list to hold the aligned blocks:
# num.grps<-length(levels(lbls))
# shifted.profiles.list<-rep(list(NULL),num.grps)
# 
# #Align the moving groups according to the parameters in sftmat:
# for(i in 1:nrow(sftmat))
#  {
#   rows.na<-as.numeric(num.reps.reduced[i])
#   pre.cols.na<-sftmat[i,6]
#   app.cols.na<-sftmat[i,7]
#   
#   shifted.grp.block<-within.aligned.dmat.dropped.ref.list[[i]]
#   #print(i)
#   #print(dim(shifted.grp.block))
#   
#   #Case 1:
#   #Lag backward, or Lag forward and get overhang
#   if(sftmat[i,2]<0 || (sftmat[i,2]>0 && sftmat[i,3]>0))
#    {
#    	#print("*CASE1")
#    	#Prepend 1. Prepend chunk needed to account for shift wrt the reference.
#     if(pre.cols.na>0)
#      {
#       nas<-rep(NA,rows.na*pre.cols.na)
#       grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
#       #print(dim(grp.pre.block))
#       #print(dim(shifted.grp.block))
#       shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
#      }
#     #Append 1. Append chunck for same thing as above.
#     if(app.cols.na>0)
#      {
#       nas<-rep(NA,rows.na*app.cols.na)
#       grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
#       #print(rows.na)
#       #print(dim(grp.app.block))
#       #print(dim(shifted.grp.block))      
#       shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
#      }
#     #Prepend 2. Extra prepend chunk need to account for shifting of the other blocks.
#     ex.na.cols<-(refp.leng + mp - ncol(shifted.grp.block))
#     if(ex.na.cols>0)
#      {
#       nas<-rep(NA,rows.na*ex.na.cols)
#       ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
#       shifted.grp.block<-cbind(ex.grp.pre.block,shifted.grp.block)
#      }
#     #Append 2. Extra append chunck for same thing as above.
#     if(ma>0)
#      {
#       nas<-rep(NA,rows.na*ma)
#       ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ma)
#       shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
#      }
#    	shifted.profiles.list[[i]]<-shifted.grp.block
#    }
#    
#    #Case 2:
#    #Lag forward but get no overhang
#    if(sftmat[i,2]>0 && sftmat[i,3]<0)
#     {
#      #print("*CASE2")
#      #Prepend 1
#      if(pre.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*pre.cols.na)
#        grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
#        shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
#       }
#      #Append 1
#      if(app.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*app.cols.na)
#        grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
#        shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
#       }
#      #Prepend 2
#      if(mp>0)
#       {
#        nas<-rep(NA,rows.na*mp)
#        ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=mp)
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.pre.block)
#       }
#      #Append 2
#      ex.na.cols<-(refp.leng + ma - ncol(shifted.grp.block))
#      if(ex.na.cols>0)
#       {
#        nas<-rep(NA,rows.na*ex.na.cols)
#        ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
#       }
#      shifted.profiles.list[[i]]<-shifted.grp.block
#     }
#     
#    #Case 3:
#    #No shifting with respect to the reference required: 
#    if(sftmat[i,2]==0) 
#     {
#      #print(paste("Group:",i))	
#      #print("*CASE3")
#      #Prepend 1
#      if(pre.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*pre.cols.na)
#        grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
#        shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
#       }
#      #Append 1
#      if(app.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*app.cols.na)
#        grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
#        shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
#       }
#      #Prepend 2
#      if(mp>0)
#       {
#        #print(mp)
#        nas<-rep(NA,rows.na*mp)
#        ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=mp)
#        #print(dim(ex.grp.pre.block))
#        #print(dim(shifted.grp.block))
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.pre.block)
#        #print(dim(shifted.grp.block))
#       }
#      #Append 2                     Tacked on lagmax to sum. Fixed problem with screwdrivers but DONT KNOW IF THIS IS ALWAYS OK!!!!!!! 
#      ex.na.cols<-(refp.leng + ma + lagmax - ncol(shifted.grp.block))
#      #print(ma)
#      #print(refp.leng)
#      #print(ncol(shifted.grp.block))
#      #print(ex.na.cols)
#      if(ex.na.cols>0)
#       {
#        nas<-rep(NA,rows.na*ex.na.cols)
#        ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
#       }
#      shifted.profiles.list[[i]]<-shifted.grp.block      
#     }
#   
#  }
# 
# #Finally pad the reference block with NA mats:
# refp.block<-within.aligned.dmat.list[[refp.idx]]
# refp.block.before<-within.aligned.dmat.list[[refp.idx]]
# 
# if(mp>0)
#  {
#   prepend.na.block<-t(replicate(nrow(refp.block),rep(NA,mp)))
#   refp.block<-cbind(prepend.na.block,refp.block) 
#  }
# if(ma>0)
#  {
#   append.na.block<-t(replicate(nrow(refp.block),rep(NA,ma)))
#   refp.block<-cbind(refp.block,append.na.block) 
#  }
# 
# shifted.profiles.list[[num.grps]]<-refp.block
# 
# #Rearrange the group blocks into their original order:
# reordered.grp.idxs<-as.numeric(reordered.grp.names)
# #print(reordered.grp.idxs)
# shifted.profiles<-NULL
# for(i in 1:length(shifted.profiles.list))
#  {
#   perm.idx<-which(reordered.grp.idxs==i)
#   #print(paste("Group:",perm.idx, dim(shifted.profiles.list[[perm.idx]]) ))
#   #print(dim(shifted.profiles.list[[perm.idx]]))
#   shifted.profiles<-rbind(shifted.profiles,shifted.profiles.list[[perm.idx]])
#  }
# 
# na.removed.shifted.profiles<-t(na.omit(data.frame(t(shifted.profiles))))
# 
# return(list(shifted.profiles,na.removed.shifted.profiles))
# 	
# }
# 
# 
# #################################################
# #Just what it says
# #################################################
# align.between.groups.with.unknown.in.a.group<-function(unk.obs.aligned, grp.idx, within.aligned.dmat.list.no.unk, lbls.no.unk, lagmax, printQ=FALSE)
# {
# 
# #Copy for modification later:
# within.aligned.dmat.list<-within.aligned.dmat.list.no.unk
# lbls<-lbls.no.unk
#   
# #Compute a grand mean for each aligned group of profiles
# grp.means<-lapply(within.aligned.dmat.list,colMeans)
# 
# #Mean vectors are all generally different in length
# #Even them out by tacking on NA vectors
# grp.means.lengs<-unlist(lapply(grp.means,length))
# max.leng<-max(grp.means.lengs)
# leng.na.vecs<-(max.leng-grp.means.lengs)
# grp.means.mat<-NULL
# for(i in 1:length(leng.na.vecs))
#  {
#   na.vec<-rep(NA,leng.na.vecs[i])
#   aug.grp.mean.vec<-c(grp.means[[i]],na.vec)
#   grp.means.mat<-rbind(grp.means.mat,aug.grp.mean.vec)
#  }
# 
# #Borrowing from align.group function.
# #Find "best" alignments between group grand means.
# #Use these shifts later to shift the group BLOCKS the same amount as the group grand mean was shifted.
# 
# #Determine an "anchor" group which will serve as the refenence with which to align the
# #other groups 
# refp.idx<-which(grp.means.lengs==max(grp.means.lengs))[1] #pick the first if there are many "longest"
# if(printQ==TRUE)
#  {
#   print("***************************************")
#   print(paste("Longest group grand mean profile is #",refp.idx,sep=""))
#   print("It is the reference")
#   print("***************************************")
#  }
# 
# #First determine where the group the unknown is aligned to is in the data (this should be a contigous ordered list).
# #The first element in this vector will be the index of the unknown in the between-group aligned data
# unk.idx<-which(lbls==grp.idx)[1]
# #Temporarily put the unknown in the specified group. It ALWAYS goes FIRST in this group:
# #within.aligned.dmat.list[[grp.idx]]<-rbind(unk.obs,within.aligned.dmat.list[[grp.idx]])
# #print(dim(within.aligned.dmat.list[[grp.idx]]))
# #print(length(unk.obs.aligned))
# within.aligned.dmat.list[[grp.idx]]<-rbind(unk.obs.aligned, within.aligned.dmat.list[[grp.idx]])
# #print(dim(within.aligned.dmat.list[[grp.idx]]))
# #print(length(lbls))
# 
# num.samps.vec<-count.group.replicates(lbls)
# num.samps.vec[grp.idx]<-(num.samps.vec[grp.idx]+1)
# print(num.samps.vec)
# lbls<-generate.label.vec(num.samps.vec)
# #lbls<-factor(c(lbls[1:(unk.idx-1)],grp.idx,lbls[unk.idx:length(lbls)]))
# 
# #print(length(lbls))
# 
# grp.names<-levels(lbls)
# reordered.grp.names<-c(grp.names[-refp.idx],grp.names[refp.idx])
# 
# #Count the number of replicates per group.
# #Drop the number of replicates in the reference group.
# num.reps<-sapply(levels(lbls),function(x){sum(lbls==x)})
# num.reps.reduced<-num.reps[-refp.idx] #Needed to replicate the shift vectors for the moving groups
# 
# 
# #Determine between group alignment by registering group grand means:
# refp<-grp.means.mat[refp.idx,]
# sftmat<-NULL
# for(i in 1:nrow(grp.means.mat))
#  {
#   if(i!=refp.idx)
#    {
#     croscor<-ccf(refp,grp.means.mat[i,],lag.max=lagmax,plot=FALSE,na.action=na.omit) #Using Cross-Correlation for alignment
#     lags<-croscor$lag
#     max.corr.idx<-which(croscor$acf[,,1]==max(croscor$acf[,,1]))
#     maxlag<-lags[max.corr.idx]
#     if(printQ==TRUE)
#      {
#       print(paste("Group grand mean profiles #",i," vs. #",refp.idx,"(ref.)",sep=""))     
#       print(paste("Index of max corr: ",max.corr.idx,sep=""))
#       print(paste("Max corr at lag: ",maxlag,sep=""))
#       print(paste("**********Max Corr: ",max(croscor$acf[,,1]),sep=""))
#      }
# 
#     lag.na.vec<-rep(NA,abs(maxlag))
# 
#     if(printQ==TRUE)
#      {
#       print(paste("Length of Reference: ",length(refp),sep=""))
#       print(paste("Length of Profile: ",length(grp.means.mat[i,]),sep=""))
#       print("XX NOW REMOVE NAs XX")
#      }
#          
#     refp<-na.omit(refp)
#     movp<-na.omit(grp.means.mat[i,])
# 
#     if(printQ==TRUE)
#      {
#       print(paste("Length of Reference NOW: ",length(refp),sep=""))
#       print(paste("Length of Profile NOW: ",length(movp),sep=""))
#      }
#     
#     dif<-length(refp)-length(movp)
# 
#     if(printQ==TRUE)
#      {
#       print(paste("Length difference: ",dif,sep=""))
#       print("================================")
#      }
#      
#     #shifts: prepend ref, append ref, prepend mov, append mov
#     sftvec<-c(0,0,0,0)        
#     if(sign(maxlag)==1) #Pushing mov profile forward
#      {
#       sftvec[3]<-maxlag #prepend lag to mov profile
#       if(sign(dif-maxlag)==1)
#        {
#          sftvec[4]<-(dif-maxlag) #if mov profile does not go past ref, append dif-lag to mov profile 
#        }
#       if(sign(dif-maxlag)==-1)
#        {
#        	sftvec[2]<-(maxlag-dif) #if mov profile goes past ref, append lag-dif to ref profile
#        }
#      }
#      
#     if(sign(maxlag)==-1) #Lagging mov profile backward
#      {
#       sftvec[1]<-abs(maxlag) #prepend lag to ref
#       sftvec[4]<-(abs(maxlag)+dif)	#append lag+dif to mov
#      }
# 
#     sftvec<-c(length(movp),sign(maxlag),sign(dif-maxlag),sftvec)        
#     sftmat<-rbind(sftmat,sftvec)
#     
#    }#end if i!=refp.idx	 	
#  }#end for
# 
# if(printQ==TRUE)
#  {
#   colnames(sftmat)<-c("len mov|","sgn lag|","sgn dif-lag|","ref pre|","ref app|","mov pre|","mov app")
#   print(sftmat)
#  }
# 
# #Build a prepend and append block of NAs for each of the shifting groups
# all.grp.names<-as.numeric(levels(lbls))
# mov.grp.names<-all.grp.names[-refp.idx]
# 
# mp<-max(sftmat[,4]) #Max prepend length
# ma<-max(sftmat[,5]) #Max append length
# refp.leng<-length(refp)
# 
# #Drop the reference group from the data to be aligned. It will be tacked on
# #to the results at the end of the process. The data will then be rearranged
# #into its original order.
# within.aligned.dmat.dropped.ref.list<-within.aligned.dmat.list[-refp.idx]
# #print(refp.idx)
# #print(length(within.aligned.dmat.list))
# #print(length(within.aligned.dmat.dropped.ref.list))
# 
# #Make an empty list to hold the aligned blocks:
# num.grps<-length(levels(lbls))
# shifted.profiles.list<-rep(list(NULL),num.grps)
# 
# #Align the moving groups according to the parameters in sftmat:
# for(i in 1:nrow(sftmat))
#  {
#   rows.na<-as.numeric(num.reps.reduced[i])
#   pre.cols.na<-sftmat[i,6]
#   app.cols.na<-sftmat[i,7]
#   
#   shifted.grp.block<-within.aligned.dmat.dropped.ref.list[[i]]
#   #print(i)
#   #print(dim(shifted.grp.block))
#   
#   #Case 1:
#   #Lag backward, or Lag forward and get overhang
#   if(sftmat[i,2]<0 || (sftmat[i,2]>0 && sftmat[i,3]>0))
#    {
#    	#print("*CASE1")
#    	#Prepend 1. Prepend chunk needed to account for shift wrt the reference.
#     if(pre.cols.na>0)
#      {
#       nas<-rep(NA,rows.na*pre.cols.na)
#       grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
#       #print(dim(grp.pre.block))
#       #print(dim(shifted.grp.block))
#       shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
#      }
#     #Append 1. Append chunck for same thing as above.
#     if(app.cols.na>0)
#      {
#       nas<-rep(NA,rows.na*app.cols.na)
#       grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
#       #print(rows.na)
#       #print(dim(grp.app.block))
#       #print(dim(shifted.grp.block))      
#       shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
#      }
#     #Prepend 2. Extra prepend chunk need to account for shifting of the other blocks.
#     ex.na.cols<-(refp.leng + mp - ncol(shifted.grp.block))
#     if(ex.na.cols>0)
#      {
#       nas<-rep(NA,rows.na*ex.na.cols)
#       ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
#       shifted.grp.block<-cbind(ex.grp.pre.block,shifted.grp.block)
#      }
#     #Append 2. Extra append chunck for same thing as above.
#     if(ma>0)
#      {
#       nas<-rep(NA,rows.na*ma)
#       ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ma)
#       shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
#      }
#    	shifted.profiles.list[[i]]<-shifted.grp.block
#    }
#    
#    #Case 2:
#    #Lag forward but get no overhang
#    if(sftmat[i,2]>0 && sftmat[i,3]<0)
#     {
#      #print("*CASE2")
#      #Prepend 1
#      if(pre.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*pre.cols.na)
#        grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
#        shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
#       }
#      #Append 1
#      if(app.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*app.cols.na)
#        grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
#        shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
#       }
#      #Prepend 2
#      if(mp>0)
#       {
#        nas<-rep(NA,rows.na*mp)
#        ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=mp)
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.pre.block)
#       }
#      #Append 2
#      ex.na.cols<-(refp.leng + ma - ncol(shifted.grp.block))
#      if(ex.na.cols>0)
#       {
#        nas<-rep(NA,rows.na*ex.na.cols)
#        ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
#       }
#      shifted.profiles.list[[i]]<-shifted.grp.block
#     }
#     
#    #Case 3:
#    #No shifting with respect to the reference required: 
#    if(sftmat[i,2]==0) 
#     {
#      #print(paste("Group:",i))	
#      #print("*CASE3")
#      #Prepend 1
#      if(pre.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*pre.cols.na)
#        grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
#        shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
#       }
#      #Append 1
#      if(app.cols.na>0)
#       {
#        nas<-rep(NA,rows.na*app.cols.na)
#        grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
#        shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
#       }
#      #Prepend 2
#      if(mp>0)
#       {
#        #print(mp)
#        nas<-rep(NA,rows.na*mp)
#        ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=mp)
#        #print(dim(ex.grp.pre.block))
#        #print(dim(shifted.grp.block))
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.pre.block)
#        #print(dim(shifted.grp.block))
#       }
#      #Append 2                     Tacked on lagmax to sum. Fixed problem with screwdrivers but DONT KNOW IF THIS IS ALWAYS OK!!!!!!! 
#      ex.na.cols<-(refp.leng + ma + lagmax - ncol(shifted.grp.block))
#      #print(ma)
#      #print(refp.leng)
#      #print(ncol(shifted.grp.block))
#      #print(ex.na.cols)
#      if(ex.na.cols>0)
#       {
#        nas<-rep(NA,rows.na*ex.na.cols)
#        ex.grp.app.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
#        shifted.grp.block<-cbind(shifted.grp.block,ex.grp.app.block)
#       }
#      shifted.profiles.list[[i]]<-shifted.grp.block      
#     }
#   
#  }
# 
# #Finally pad the reference block with NA mats:
# refp.block<-within.aligned.dmat.list[[refp.idx]]
# refp.block.before<-within.aligned.dmat.list[[refp.idx]]
# 
# if(mp>0)
#  {
#   prepend.na.block<-t(replicate(nrow(refp.block),rep(NA,mp)))
#   refp.block<-cbind(prepend.na.block,refp.block) 
#  }
# if(ma>0)
#  {
#   append.na.block<-t(replicate(nrow(refp.block),rep(NA,ma)))
#   refp.block<-cbind(refp.block,append.na.block) 
#  }
# 
# shifted.profiles.list[[num.grps]]<-refp.block
# 
# #Rearrange the group blocks into their original order:
# reordered.grp.idxs<-as.numeric(reordered.grp.names)
# #print(reordered.grp.idxs)
# shifted.profiles<-NULL
# for(i in 1:length(shifted.profiles.list))
#  {
#   perm.idx<-which(reordered.grp.idxs==i)
#   #print(paste("Group:",perm.idx, dim(shifted.profiles.list[[perm.idx]]) ))
#   #print(dim(shifted.profiles.list[[perm.idx]]))
#   shifted.profiles<-rbind(shifted.profiles,shifted.profiles.list[[perm.idx]])
#  }
# 
# na.removed.shifted.profiles<-t(na.omit(data.frame(t(shifted.profiles))))
# 
# #Now separate the unknown from the aligned data set:
# #print(dim(na.removed.shifted.profiles))
# aligned.unk<-na.removed.shifted.profiles[unk.idx,]
# na.removed.shifted.profiles<-na.removed.shifted.profiles[-unk.idx,]
# #print(dim(na.removed.shifted.profiles))
# 
# return(list(shifted.profiles,na.removed.shifted.profiles,aligned.unk))
# 	
# }


#
#
# align.between.groups.with.unknown.in.a.group<-function(unk.obs.aligned, grp.idx, X.aligned.list.dat, lbls, lagmax, anchor.grp.idx, printQ=FALSE)
# {
# 
# grp.mns<-lapply(X.aligned.list.dat,colMeans)
# grp.means.lengs<-unlist(lapply(grp.mns,length))
# 
# grp.names<-levels(lbls)
# print(grp.names)
# 
# if(anchor.grp.idx=="longest")
#  {
#   refp.idx<-which(grp.means.lengs==max(grp.means.lengs))[1] #pick the first if there are many "longest"
#   if(printQ==TRUE)
#    {
#     print("***************************************")
#     print(paste("Longest group grand mean profile is #",refp.idx," (Group: ",grp.names[refp.idx],")",sep=""))
#     print("It is the reference")
#     print("***************************************")
#    } 
#  }
# if(class(anchor.grp.idx)=="numeric")
#  {
#   refp.idx<-anchor.grp.idx
#   if(printQ==TRUE)
#    {
#     print("***************************************")
#     print(paste("Aligning All Groups to CHOSEN group grand mean#",refp.idx," (Group: ",grp.names[refp.idx],")",sep=""))
#     print("***************************************")
#    }   
#  }
# refp<-grp.mns[[refp.idx]]
# #print(refp.idx)
# 
# X.aligned.list.dat.mod<-rep(list(NULL),length(X.aligned.list.dat))
# sftmat<-NULL
# for(i in 1:length(grp.mns))
#  {
#   if(i!=refp.idx)
#    {
#     mov.list<-align.mov.to.ref(refp,grp.mns[[i]],lagmax,padtyp=list("NAs"),printQ)
#     sftinfo<-mov.list[[2]]
#     #sftinfo format c(ref shorter/longer -1/1, shift right/left 1/-1, num.pad.left, num.chop.left, num.pad.right, num.chop.right)
#     sftmat<-rbind(sftmat,sftinfo)
#     mov.chop.left<-sftinfo[4]
#     mov.chop.right<-sftinfo[6]
#     mov<-X.aligned.list.dat[[i]]
#     
#     if(mov.chop.right>0)
#      {
#       idxs.drop.right<-sort((dim(mov)[2]):((dim(mov)[2])-(mov.chop.right-1)))
#       mov.chop<-mov[,-idxs.drop.right]         #chop moving group on right
#       #unk.obs.mov<-unk.obs.mov[-idxs.drop.right]       #remove padding
#      }
#     if(mov.chop.left>0)
#      {
#       idxs.drop.left<-(1:mov.chop.left)
#       mov.chop<-mov[,-idxs.drop.left]         #chop moving group on left
#       #unk.obs.mov<-unk.obs.mov[-idxs.drop.left]       #remove padding
#      }
#     if((mov.chop.left==0) & (mov.chop.right==0))
#      {
#       mov.chop<-mov.grp
#      }
#     print(length(refp))
#     print(dim(mov))
#     print(dim(mov.chop))
#    }
#  }
# colnames(sftmat)<-c("Ref shorter/longer (-1/1)", "Shift right/left (1/-1)", "Pad Left", "Chop Left", "Pad Right", "Chop Right")
# print(sftmat)
# 
# 
# }