#################################################################
#"Best" Align a set of profile vectors.
#"Best is cross-correlation"
#In the code below:
#ref is the stationary (longest) profile of the group. 
#mov is one of the other (shorter) profiles that is sfhited
#################################################################
align.profile.group.padrandom<-function(profile.group,lagmax,printQ=FALSE)
{

nump<-dim(profile.group)[1]  
plns<-as.numeric(lapply(apply(profile.group,1,na.omit),length))

refp.idx<-which(plns==max(plns))[1] #pick the first if there are many "longest"
if(printQ==TRUE)
 {
  print("***********************************")
  print(paste("Longest profile is #",refp.idx,sep=""))
  print("It is the reference for this group")
  print("***********************************")
 }

refp<-profile.group[refp.idx,]
sftmat<-NULL
for(i in 1:nump)
 {
  if(i!=refp.idx)
   {
    croscor<-ccf(refp,profile.group[i,],lag.max=lagmax,plot=FALSE,na.action=na.omit) #Using Cross-Correlation for alignment
    lags<-croscor$lag
    max.corr.idx<-which(croscor$acf[,,1]==max(croscor$acf[,,1]))
    maxlag<-lags[max.corr.idx]
    if(printQ==TRUE)
     {
      print(paste("Profiles #",i," vs. #",refp.idx,"(ref.)",sep=""))     
      print(paste("Index of max corr: ",max.corr.idx,sep=""))
      print(paste("Max corr at lag: ",maxlag,sep=""))
      print(paste("Max Corr: ",max(croscor$acf[,,1]),sep=""))
      #print("================================")
     }

    lag.na.vec<-rep(NA,abs(maxlag))

    if(printQ==TRUE)
     {
      print(paste("Length of Reference: ",length(refp),sep=""))
      print(paste("Length of Profile: ",length(profile.group[i,]),sep=""))
      print("XX NOW REMOVE NAs XX")
     }
         
    refp<-na.omit(refp)
    movp<-na.omit(profile.group[i,])

    if(printQ==TRUE)
     {
      print(paste("Length of Reference NOW: ",length(refp),sep=""))
      print(paste("Length of Profile NOW: ",length(movp),sep=""))
     }
    
    dif<-length(refp)-length(movp)
#    print("+_+_+_+_+_+_+_+_+_DIAGNOSTICS_+_+_+_+_+_+_+")
#    print(i)
#    print(dif)
#    print(maxlag)
#    print(dif-maxlag)
#    print(sign(dif-maxlag))
#    print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

    if(printQ==TRUE)
     {
      print(paste("Length difference: ",dif,sep=""))
      print("================================")
     }
     
    #shifts: prepend ref, append ref, prepend mov, append mov
    sftvec<-c(0,0,0,0)
    sgn.dif.lag<-sign(dif-maxlag)#***********ADDED HERE!!!!!!!!
    if(sgn.dif.lag==0)
     {
      sgn.dif.lag<-1
     }

    if(sign(maxlag)==1) #Pushing mov profile forward
     {      
      sftvec[3]<-maxlag #prepend lag to mov profile
      #if(sign(dif-maxlag)==1)
      if(sgn.dif.lag==1) #CHANGED!!!!!!!!!!
       {
       	sftvec[4]<-(dif-maxlag) #if mov profile does not go past ref, append dif-lag to mov profile 
       }
      #if(sign(dif-maxlag)==-1)
      if(sgn.dif.lag==-1) #CHANGED!!!!!!!!!
       {
       	sftvec[2]<-(maxlag-dif) #if mov profile goes past ref, append lag-dif to ref profile
       }
     }
     
    if(sign(maxlag)==-1) #Lagging mov profile backward
     {
      sftvec[1]<-abs(maxlag) #prepend lag to ref
      sftvec[4]<-(abs(maxlag)+dif)	#append lag+dif to mov
     }

    #sftvec<-c(length(movp),sign(maxlag),sign(dif-maxlag),sftvec)
    sftvec<-c(length(movp),sign(maxlag),sgn.dif.lag,sftvec) #CHANGED!!!!!!!!!    
    sftmat<-rbind(sftmat,sftvec)
    
   }#end if i!=refp.idx	

 	
 }#end for

if(printQ==TRUE)
 {
  colnames(sftmat)<-c("len mov|","sgn lag|","sgn dif-lag|","ref pre|","ref app|","mov pre|","mov app")
  print(sftmat)
 }
 
if(printQ=="SPECIAL")
 {
  colnames(sftmat)<-c("len mov|","sgn lag|","sgn dif-lag|","ref pre|","ref app|","mov pre|","mov app")
  print(sftmat)
 }


#Now do shifting by prepending and appending appropriate amounts:
mp<-max(sftmat[,4]) #Max prepend length
ma<-max(sftmat[,5]) #Max append length
#print(mp)
#print(ma)

#Prep the shifting profiles by dropping the NAs
movplist<-apply(profile.group,1,na.omit)
#IF apply spits out a matrix instead of a list make it into a list.
#This can happen if all the profiles are the same length (I think) => get matrix instead of a list 
if(class(movplist)=="matrix")
 {
  #print("PROBLEM?!")
  movplist<-lapply(1:ncol(movplist),function(x){movplist[,x]})
 }
#print(movplist)
movplist<-movplist[-refp.idx] #This drops the reference profile from the list! NOTE: [] used instead of [[]] !! 
#print(length(movplist))
#print(refp.idx)
#print(as.numeric(lapply(movplist,length)))

#NEED THIS YET??
#refp<-c(rep(NA,mp),refp) #
#print(length(refp)) #At this point NAs should have been removed from reference profile
#print(refp)

#mod.profiles<-NULL #Old way
#Make a list for the modified profiles in permuted order. This list will be rearranged back to the original order
#after processing is finished.
mod.profiles.list<-rep(list(NULL),nump)
for(i in 1:dim(sftmat)[1])
 {
  if(sftmat[i,2]<0 || (sftmat[i,2]>0 && sftmat[i,3]>0))
   {
   	#print("*CASE1")
   	modmovp<-c(runif(sftmat[i,6],min=0,max=0.0001),movplist[[i]],runif(sftmat[i,7],min=0,max=0.0001)) #prepend and append appropriate amounts to mov profile
   	modmovp<-c(runif((length(refp)+mp-length(modmovp)),min=0,max=0.0001),modmovp,runif(ma,min=0,max=0.0001))
   	#print(length(refp)-length(modmovp))
   	#print(length(modmovp))
#   	mod.profiles<-rbind(mod.profiles,modmovp) #Old way
   	mod.profiles.list[[i]]<-modmovp
   }
  
  if(sftmat[i,2]>0 && sftmat[i,3]<0) 
   {
    #print("*CASE2")
    modmovp<-c(runif(sftmat[i,6],min=0,max=0.0001),movplist[[i]],runif(sftmat[i,7],min=0,max=0.0001)) #prepend and append appropriate amounts to mov profile
    modmovp<-c(runif(mp,min=0,max=0.0001),modmovp,runif((length(refp)+ma-length(modmovp)),min=0,max=0.0001))
#    mod.profiles<-rbind(mod.profiles,modmovp) #Old way
    mod.profiles.list[[i]]<-modmovp
   }

  if(sftmat[i,2]==0) 
   {
   	#print("*CASE3)
    #print("HEY IN HERE! NO SHIFT")
    modmovp<-c(runif(sftmat[i,6],min=0,max=0.0001),movplist[[i]],runif(sftmat[i,7],min=0,max=0.0001)) #prepend and append appropriate amounts to mov profile
    modmovp<-c(runif(mp,min=0,max=0.0001),modmovp,runif((length(refp)+ma-length(modmovp)),min=0,max=0.0001))
#    mod.profiles<-rbind(mod.profiles,modmovp) #Old way
    mod.profiles.list[[i]]<-modmovp
   }
      
 }

refp<-c(runif(mp,min=0,max=0.0001),refp,runif(ma,min=0,max=0.0001)) #finally pad reference 
#mod.profiles<-rbind(refp,mod.profiles) #Old way
mod.profiles.list[[nump]]<-refp

#Now rearrange the profiles into their original order:
original.order.idxs<-1:nump
permuted.order.idxs<-c(original.order.idxs[-refp.idx],refp.idx)
#print(original.order.idxs)
#print(permuted.order.idxs)
mod.profiles<-NULL
#print(nump)
for(i in 1:nump)
 {
  #print(i)
  perm.idx<-which(permuted.order.idxs==i)
  #print(paste(i,perm.idx))
  #print(dim(mod.profiles))
  #print(length(mod.profiles.list[[perm.idx]]))
  mod.profiles<-rbind(mod.profiles,mod.profiles.list[[perm.idx]])
  #print(dim(mod.profiles))
  #print("**")
 }

#print(mod.profiles)
#print(nump)
#print(permuted.order.idxs)
#print(length(mod.profiles.list))
print(dim(mod.profiles))
print(dim(t(na.omit(t(mod.profiles)))))


#mod.profiles<-t(na.omit(data.frame(t(mod.profiles)))) #DONT CUT THESE. SHOULDNT HAVE ANY NAs!
rownames(mod.profiles)<-NULL
colnames(mod.profiles)<-NULL
return(mod.profiles)
	
}