align.between.groups.with.unknown.in.a.group<-function(unk.obs.aligned, grp.idx, within.aligned.list, lbls, maxlag, anchor.grp.idx, printQ=FALSE)
{

grp.mns<-lapply(within.aligned.list,colMeans)
grp.means.lengs<-unlist(lapply(grp.mns,length))

grp.names<-levels(lbls)
#print(grp.names)

if(anchor.grp.idx=="longest")
 {
  refp.idx<-which(grp.means.lengs==max(grp.means.lengs))[1] #pick the first if there are many "longest"
  if(printQ==TRUE)
   {
    print("***************************************")
    print(paste("Longest group grand mean profile is #",refp.idx," (Group: ",grp.names[refp.idx],")",sep=""))
    print("It is the reference")
    print("***************************************")
   } 
 }
if(class(anchor.grp.idx)=="numeric")
 {
  refp.idx<-anchor.grp.idx
  if(printQ==TRUE)
   {
    print("***************************************")
    print(paste("Aligning All Groups to CHOSEN group grand mean#",refp.idx," (Group: ",grp.names[refp.idx],")",sep=""))
    print("***************************************")
   }   
 }
refp<-grp.mns[[refp.idx]]

#Note: the reference group should stay the same.
within.between.aligned.list<-within.aligned.list
#print(sapply(1:length(within.between.aligned.list),function(x){dim(within.between.aligned.list[[x]])}))

#Put the unknown onto its assigned group, to be moved:
within.between.aligned.list[[grp.idx]]<-rbind(unk.obs.aligned,within.between.aligned.list[[grp.idx]])
#print(sapply(1:length(within.between.aligned.list),function(x){dim(within.between.aligned.list[[x]])}))

sftmat<-NULL
for(i in 1:length(grp.mns))
 {
  if(i!=refp.idx)
   {
    mov.list<-align.mov.to.ref(refp,grp.mns[[i]],maxlag,padtyp=list("NAs"),printQ)
    sftinfo<-mov.list[[2]]
    #sftinfo format c(ref shorter/longer -1/1, shift right/left 1/-1, num.pad.left, num.chop.left, num.pad.right, num.chop.right)
    sftmat<-rbind(sftmat,sftinfo)
    mov.chop.left<-sftinfo[4]
    mov.chop.right<-sftinfo[6]
    mov.pad.left<-sftinfo[3]
    mov.pad.right<-sftinfo[5]
    
    mov<-within.between.aligned.list[[i]]
    #print(length(refp))
    #print(dim(mov))
    
    if(mov.chop.left>0)
     {
      idxs.drop.left<-(1:mov.chop.left)
      mov<-mov[,-idxs.drop.left]         #chop moving group on left
     }
    if(mov.chop.right>0)
     {
      idxs.drop.right<-sort((dim(mov)[2]):((dim(mov)[2])-(mov.chop.right-1)))
      mov<-mov[,-idxs.drop.right]         #chop moving group on right
     }
    if(mov.pad.right>0)
     {
      rws<-nrow(mov)
      cls<-mov.pad.right
      right.pad<-matrix(padding("NAs",rws*cls),nrow=rws,ncol=cls)
      mov<-cbind(mov,right.pad)
     }
    if(mov.pad.left>0)
     {
      rws<-nrow(mov)
      cls<-mov.pad.left
      left.pad<-matrix(padding("NAs",rws*cls),nrow=rws,ncol=cls)
      mov<-cbind(left.pad,mov)
     }
    #print(dim(mov))
    within.between.aligned.list[[i]]<-mov
   }
 }
#colnames(sftmat)<-c("Ref shorter/longer (-1/1)", "Shift right/left (1/-1)", "Pad Left", "Chop Left", "Pad Right", "Chop Right")
#print(sftmat)
#print(sapply(1:length(within.between.aligned.list),function(x){dim(within.between.aligned.list[[x]])}))
grp.idx.mat<-within.between.aligned.list[[grp.idx]]
#print(dim(within.between.aligned.list[[grp.idx]]))
unk.obs.wb.aligned<-grp.idx.mat[1,]
#print(dim(grp.idx.mat))
grp.idx.mat<-grp.idx.mat[-1,]
within.between.aligned.list[[grp.idx]]<-grp.idx.mat
#print(dim(within.between.aligned.list[[grp.idx]]))

within.between.aligned.mat<-NULL
for(i in 1:length(within.between.aligned.list))
 {
   within.between.aligned.mat<-rbind(within.between.aligned.mat,within.between.aligned.list[[i]])
 }
#print(dim(within.between.aligned.mat))
within.between.aligned.mat<-rbind(unk.obs.wb.aligned,within.between.aligned.mat)
#print(dim(within.between.aligned.mat))
within.between.aligned.mat<-t(na.omit(data.frame(t(within.between.aligned.mat))))
unk.obs.wb.aligned<-within.between.aligned.mat[1,]
within.between.aligned.mat<-within.between.aligned.mat[-1,]
#print(dim(within.between.aligned.mat))

return(list(unk.obs.wb.aligned,within.between.aligned.mat))

}
