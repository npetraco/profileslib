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
align.between.profile.padrandom<-function(within.aligned.dmat.list,lbls,lagmax,printQ=FALSE)
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
      nas<-runif(rows.na*pre.cols.na,min=0,max=0.0001)
      grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
      #print(dim(grp.pre.block))
      #print(dim(shifted.grp.block))
      shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
     }
    #Append 1. Append chunck for same thing as above.
    if(app.cols.na>0)
     {
      nas<-runif(rows.na*app.cols.na,min=0,max=0.0001)
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
      nas<-runif(rows.na*ex.na.cols,min=0,max=0.0001)
      ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=ex.na.cols) 
      shifted.grp.block<-cbind(ex.grp.pre.block,shifted.grp.block)
     }
    #Append 2. Extra append chunck for same thing as above.
    if(ma>0)
     {
      nas<-runif(rows.na*ma,min=0,max=0.0001)
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
       nas<-runif(rows.na*pre.cols.na,min=0,max=0.0001)
       grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
       shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
      }
     #Append 1
     if(app.cols.na>0)
      {
       nas<-runif(rows.na*app.cols.na,min=0,max=0.0001)
       grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
       shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
      }
     #Prepend 2
     if(mp>0)
      {
       nas<-runif(rows.na*mp,min=0,max=0.0001)
       ex.grp.pre.block<-matrix(nas,nrow=rows.na,ncol=mp)
       shifted.grp.block<-cbind(shifted.grp.block,ex.grp.pre.block)
      }
     #Append 2
     ex.na.cols<-(refp.leng + ma - ncol(shifted.grp.block))
     if(ex.na.cols>0)
      {
       nas<-runif(rows.na*ex.na.cols,min=0,max=0.0001)
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
       nas<-runif(rows.na*pre.cols.na,min=0,max=0.0001)
       grp.pre.block<-matrix(nas,nrow=rows.na,ncol=pre.cols.na)
       shifted.grp.block<-cbind(grp.pre.block,shifted.grp.block)
      }
     #Append 1
     if(app.cols.na>0)
      {
       nas<-runif(rows.na*app.cols.na,min=0,max=0.0001)
       grp.app.block<-matrix(nas,nrow=rows.na,ncol=app.cols.na)
       shifted.grp.block<-cbind(shifted.grp.block,grp.app.block)
      }
     #Prepend 2
     if(mp>0)
      {
       #print(mp)
       nas<-runif(rows.na*mp,min=0,max=0.0001)
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
       nas<-runif(rows.na*ex.na.cols,min=0,max=0.0001)
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
  prepend.na.block<-t(replicate(nrow(refp.block),runif(mp,min=0,max=0.0001)))
  refp.block<-cbind(prepend.na.block,refp.block) 
 }
if(ma>0)
 {
  append.na.block<-t(replicate(nrow(refp.block),runif(ma,min=0,max=0.0001)))
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

#na.removed.shifted.profiles<-t(na.omit(data.frame(t(shifted.profiles))))
na.removed.shifted.profiles<-shifted.profiles

return(list(shifted.profiles,na.removed.shifted.profiles))
	
}