align.within.a.group<-function(dat.group, anchor.profile.idx, lagmax, padtyp, printQ=FALSE)
{                              
  
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
    info.list<-align.mov.to.ref.pad.both(refp,movp,lagmax,padtyp,printQ)
    pdref<-info.list[[1]]
    pdmov<-info.list[[2]]
    padded.mov.list[[i]]<-pdmov
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
#print(sapply(1:nrow(dat.group), function(x){length(padded.mov.list[[x]])}))
for(i in 1:nrow(smat))
 {
  ex.pad.mov.left<-NULL  #Making sure to clear out
  ex.pad.mov.right<-NULL #Making sure to clear out
  movp<-NULL             #Making sure to clear out
  ex.pad.mov.left<-(ref.max.left-smat[i,5])
  ex.pad.mov.right<-(ref.max.right-smat[i,7])
  movp.idx<-smat[i,1]
  movp<-c(padding(padtyp,ex.pad.mov.left), padded.mov.list[[movp.idx]], padding(padtyp,ex.pad.mov.right))
  #print(length(movp))
  padded.mov.list[[movp.idx]]<-movp 
 }
#Now pad the reference:
refp<-na.omit(dat.group[refp.idx,])
refp<-c(padding(padtyp,ref.max.left), refp, padding(padtyp,ref.max.right))
padded.mov.list[[refp.idx]]<-refp
dat.within.group.aligned<-t(sapply(1:nrow(dat.group), function(x){padded.mov.list[[x]]}))

return(dat.within.group.aligned)

}