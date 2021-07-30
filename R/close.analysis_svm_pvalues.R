################################
#LOOP OVER "CLOSE" GROUPS:
################################
close.analysis_svm.pvalues<-function(unk.obs.vec, close.groups, within.aligned.close.list, var.tol, unk.maxlag.within, maxlag.between)
{
  
#im.mat<-NULL
pvalue.vec<-NULL
for(i in 1:length(close.grps))
 {
  #Now we have to check the classification results for the unknown aligned against EACH of the "close" groups
  specific.grp<-close.grps[i]
  print("*********************************************************")
  print(paste("Iteration: ",i," Alignment with group: ", close.groups[i]))
  print("*********************************************************")
  specific.grp.idx<-which(close.grps==specific.grp)
  tmp<-align.unknown.to.group(unk.obs.vec, specific.grp.idx, within.aligned.close.list, lagmax=unk.maxlag.within, printQ=FALSE, plotQ=FALSE)
  unk.obs.mov<-tmp[[1]] #The unknown, aligned to one of the "close" groups
  within.aligned.close.list.tmp<-tmp[[2]] #The within group aligned data

  #Align between groups
  within.between.aligned.list.tmp<-align.between.groups.with.unknown.in.a.group(unk.obs.mov, specific.grp.idx, within.aligned.close.list.tmp, lbl.close, maxlag.between, "longest", printQ=FALSE)
  Xdat<-within.between.aligned.list.tmp[[2]]                   #Within and Between group aligned data
  wb.aligned.unk.obs.vec<-within.between.aligned.list.tmp[[1]] #Within and Between group aligned unknown, wrt close group i

  num.pts<-dim(Xdat)[2] #number of points in the profiles used, once alighed
  
  Xdat<-scale(Xdat,center=TRUE,scale=FALSE)[,]
  pca.model<-prcomp(Xdat,scale=FALSE)
  #im.mat<-cbind(im.mat,summary(pca.model)$importance[3,1:50])
  Mmax<-which(summary(pca.model)$importance[3,]>=var.tol)[1]
  print("")
  print(paste(100*var.tol,"% variance occurs at dimension: ", Mmax, sep=""))
  print("Begin HOO-CV model dimension determination.")

  #Find optimal discrimination dimension with HOO-CV
  err.vec<-NULL
  ind.mat<-NULL
  for(j in 2:Mmax)
   {
    Z<-predict(pca.model)[,1:j]

    ind.vec<-NULL
    for(k in 1:nrow(Z))
     {
      Z.heldout<-t(as.matrix(Z[k,]))
      lbl.heldout<-lbl.close[k]
  
      Z.kept<-Z[-k,]
      lbl.kept<-lbl.close[-k]
      svm.model<-svm(Z.kept,lbl.kept,scale=FALSE,type="C-classification",kernel="linear",cost=0.1,fitted=TRUE,probability=TRUE)
      pred<-predict(svm.model,Z.heldout)
      #prob.vec<-attr(pred, "probabilities")[,]
      ind.vec<-c(ind.vec,pred==lbl.heldout)

     } #end for k
    
    ind.mat<-cbind(ind.mat,ind.vec)  
    ccp<-(sum(ind.vec)/nrow(Z) )
    err<-(1-ccp)*100
    print(paste(j,err))
    err.vec<-c(err.vec,err)

   } #end for j
  cv.err.mat<-cbind(2:Mmax,err.vec)
  colnames(cv.err.mat)<-c("Dimension", "HOO-CV error")
  print(cv.err.mat)
  print("")
  Mmin<-(which(err.vec==min(err.vec))+1)[1]
  print(paste("Minimal error dimension is: ", Mmin,"D. Minimum HOO-CV error is: ",min(err.vec),sep=""))
  #plot(2:Mmax,err.vec,typ="l")

  #Predict unknown with chosen dimension:
  Mpred<-Mmin
  Z<-predict(pca.model)[,1:Mpred]                                   #Grab PCA scores
  
  #Project unknown into mimimal dimension space:
  Apc<-pca.model$rotation[,1:Mpred]
  Zunk.vec<-wb.aligned.unk.obs.vec%*%Apc

  #Put the unknown into the "bag" with lbl.unk = specific.grp (ie the unk labeled as the group it is aligned on)
  Zaug<-rbind(Z,Zunk.vec)
  lblaug<-as.factor(c(as.character(lbl.close), as.character(specific.grp)))
  
  #Determine nonconformity measures of each item in the bag, cf. Algorithmic Learning in a Random World and Tutorial on Conformal Prediction
  noncf.vec<-NULL
  for(kk in 1:nrow(Zaug))
   {
    Zaug.heldout<-t(as.matrix(Zaug[kk,]))
    lblaug.heldout<-lblaug[kk]
  
    Zaug.kept<-Zaug[-kk,]
    lblaug.kept<-lblaug[-kk]
    svm.model<-svm(Zaug.kept,lblaug.kept,scale=FALSE,type="C-classification",kernel="linear",cost=0.1,fitted=TRUE,probability=TRUE)
    #pred<-predict(svm.model,Zaug.heldout)
    grp.platt.scores<-attr(predict(svm.model, Zaug.heldout, probability=TRUE), "probabilities")[,]
    lblaug.heldout.idx<-which(names(grp.platt.scores)==as.character(lblaug.heldout))
    platt.pred.score<-grp.platt.scores[lblaug.heldout.idx]
    noncf.vec<-c(noncf.vec,platt.pred.score)
   } #end for kk
  pvalue<-sum(noncf.vec>=noncf.vec[length(noncf.vec)])/length(noncf.vec)
  print(paste("Assignment group",specific.grp,"p-value:",pvalue))
  pvalue.vec<-c(pvalue.vec,pvalue)
  
 } #end for i
print(cbind(as.character(close.groups),as.character(pvalue.vec)))
print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

}