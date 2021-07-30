align.between.groups<-function(dat.wga.list, lbls, lagmax, ref.choice, paddingtype, printQ=FALSE)
{ 

grp.means.list<-lapply(1:length(dat.wga.list),function(x){colMeans(dat.wga.list[[x]],na.rm=TRUE)})

pro.lengs<-sapply(1:length(grp.means.list),function(x){length(grp.means.list[[x]])})
leng.diffs<-max(pro.lengs)-pro.lengs

#Pad length diffs with NAs and see if code for align within a group can be used
padded.grp.means<-NULL
for(i in 1:length(grp.means.list))
 {
  tmp.pad<-rep(NA,leng.diffs[i])
  padded.grp.mean<-c(grp.means.list[[i]],tmp.pad)
  padded.grp.means<-rbind(padded.grp.means,padded.grp.mean)
 }

dat.group<-padded.grp.means #rename to use align.within.a.group code
anchor.profile.idx<-ref.choice

#If the longest profile is requested as the reference:  
if(anchor.profile.idx=="longest")
 {
  profile.lengths<-sapply(1:nrow(dat.group),function(x){length(na.omit(dat.group[x,]))})
  refp.idx<-which(profile.lengths==max(profile.lengths))[1] #pick the first if there are many "longest"

  if(printQ==TRUE)
   {
    print("***********************************")
    print(paste("Longest profile is #",refp.idx,sep=""))
    print("It is the reference for this group")
    print("***********************************")
   } 
 }
#If the best (highest average) CCF max is requested as the reference:
max.ccfs<-NULL
if(anchor.profile.idx=="best ccf")
 {
  print("CAUTION! This will be slow.")
  for(i in 1:nrow(dat.group))
   {
    #avg.max.ccf<-0 
    tmp<-0          #initialize for prifile i
    for(j in 1:nrow(dat.group))
     {
      if(i!=j)
       {
        croscor<-ccf(dat.group[i,], dat.group[j,], lag.max=lagmax,plot=FALSE,na.action=na.omit)
        tmp<-(tmp+max(croscor$acf[,,1]))
       }
     }
    avg.max.ccf<-tmp/(nrow(dat.group)-1)
    print(avg.max.ccf)
    max.ccfs<-c(max.ccfs,avg.max.ccf)
   }

  refp.idx<-which(max.ccfs==max(max.ccfs))[1]
  if(printQ==TRUE)
   {
    print("***************************************")
    print(paste("Aligning All Profiles to Best Avg CCF Max reference profile# ",refp.idx,sep=""))
    print(paste("Higest Average Max CCF: ",max(max.ccfs),sep=""))
    print("***************************************")
   }   
 }
#You pick the reference:
if(class(anchor.profile.idx)=="numeric")
 {
  refp.idx<-anchor.profile.idx
  if(printQ==TRUE)
   {
    print("***************************************")
    print(paste("Aligning All Profiles to CHOSEN reference profile#",refp.idx,sep=""))
    print("***************************************")
   }   
 }

#Grab the reference profile for the group:
refp<-dat.group[refp.idx,]

padded.mov.list<-rep(list(NULL),nrow(dat.group)) #container to hold the shifting profiles
smat<-NULL
for(i in 1:nrow(dat.group))
 {
  if(i !=refp.idx)
   {
    movp<-dat.group[i,]
    info.list<-align.mov.to.ref.pad.both(refp,movp,lagmax,padtyp=paddingtype,printQ)
    pdref<-info.list[[1]]
    pdmov<-info.list[[2]]
    #padded.mov.list[[i]]<-pdmov
    svec<-info.list[[3]]
    svec<-c(i,svec, length(pdref), length(pdmov))
    smat<-rbind(smat,svec)
    #print("=========")
    #print(length(pdref))
    #print(length(pdmov))
    #print("=========")
   }
 }
colnames(smat)<-c("Mov idx |","Ref shorter/longer (-1/1) |", "Shift right/left (1/-1) |", "Pad Mov Left |", "Pad Ref Left |", "Pad Mov Right |", "Pad Ref Right |", "Length PdMov |", "Length PdRef")
#rownames(smat)<-NULL
if(printQ==TRUE)
 {
  print(smat)  
 }

mov.max.left<-max(smat[,4])
mov.max.right<-max(smat[,6])
ref.max.left<-max(smat[,5])
ref.max.right<-max(smat[,7])

max.left<-max(c(mov.max.left, ref.max.left))
max.right<-max(c(mov.max.right, ref.max.right))

#Do the padding on the moving profiles:
for(i in 1:nrow(smat))
 {
  ex.pad.mov.left<-NULL  #Making sure to clear out
  ex.pad.mov.right<-NULL #Making sure to clear out
  movp<-NULL             #Making sure to clear out
  ex.pad.mov.left<-(ref.max.left-smat[i,5]+smat[i,4])
  ex.pad.mov.right<-(ref.max.right-smat[i,7]+smat[i,6])
  movp.idx<-smat[i,1]
  
  movp<-dat.wga.list[[movp.idx]]
  if(ex.pad.mov.left>0)
   {
    #Make the left padding block for the moving group
    left.pad.mat<-matrix(rep(padding(paddingtype,ex.pad.mov.left),nrow(movp)), nrow=nrow(movp),ncol=ex.pad.mov.left)
    movp<-cbind(left.pad.mat,movp)
   }
  if(ex.pad.mov.right>0)
   {
    #Make the right padding block for the moving group
    right.pad.mat<-matrix(rep(padding(paddingtype,ex.pad.mov.right),nrow(movp)), nrow=nrow(movp),ncol=ex.pad.mov.right)
    movp<-cbind(movp,right.pad.mat)
   }
  padded.mov.list[[movp.idx]]<-movp
 }

#Now pad the reference:
refp<-dat.wga.list[[refp.idx]]
if(ref.max.left>0)
 {
  left.ref.pad.mat<-matrix(rep(padding(paddingtype,ref.max.left),nrow(refp)), nrow=nrow(refp),ncol=ref.max.left)
  refp<-cbind(left.ref.pad.mat,refp)
 }
if(ref.max.right>0)
 {
  right.ref.pad.mat<-matrix(rep(padding(paddingtype,ref.max.right),nrow(refp)), nrow=nrow(refp),ncol=ref.max.right)
  refp<-cbind(refp,right.ref.pad.mat)
 }
padded.mov.list[[refp.idx]]<-refp

#Tack all the matricex into one big matrix:
dat.within.group.aligned<-NULL
for(i in 1:length(padded.mov.list))
 {
  dat.within.group.aligned<-rbind(dat.within.group.aligned,padded.mov.list[[i]])
 }
#print(dim(dat.within.group.aligned))
 
return(dat.within.group.aligned)

}
