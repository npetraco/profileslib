#library(rgl)

#-------------------------------------------------------------------
#Compute univariate Fisher ratios between two samples
#Reference: Sahota and Morgan, Anal.Chem.64,2383,1992
#-------------------------------------------------------------------
fisher.ratios<-function(dmat,lbls,grp1,grp2,plotQ=TRUE,plot.typ="l")
{

grp1.idxs<-which(lbls==grp1)
grp2.idxs<-which(lbls==grp2)

means.g1<-colMeans(dmat[grp1.idxs,])
means.g2<-colMeans(dmat[grp2.idxs,])

var.g1<-apply(dmat[grp1.idxs,],2,var)
var.g2<-apply(dmat[grp2.idxs,],2,var)
#var.g1<-var(dmat[grp1.idxs,])
#var.g2<-var(dmat[grp2.idxs,])

ratios<-((means.g1-means.g2)^2)/(var.g1+var.g2)

if(plotQ==TRUE)
 {
  plot(1:length(ratios),ratios,type=plot.typ)
 }

return(ratios)

}


#-------------------------------------------------------------------
#3D rotatable bar plot. Good for visualizing loadings across
#many variables.
#-------------------------------------------------------------------
# loading.bar3D<-function(loading.mat,bar.width)
# {
#
# dr<-dim(loading.mat)[1]
# dc<-dim(loading.mat)[2]
#
# bp.mat<-NULL
# for(i in 1:dr)
#  {
#   for(j in 1:dc)
#    {
#     bp.vec<-c(i,j,loading.mat[i,j])
#     bp.mat<-rbind(bp.mat,bp.vec)
#    }
#  }
# bp.mat<-abs(bp.mat)
#
# #X-axis is loadings. Y-axis is variables.
#
# yLabels <- c(1:nrow(loading.mat))
# xLabels <- rownames(loading.mat)
#
# plot3d(bp.mat[,1],bp.mat[,2],bp.mat[,3],axes=T,type="h",lwd=bar.width,xlab="Variables",ylab="L.V. #",zlab="|Loading|")
# axis3d(edge='x-', at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
# axis3d(edge='y-',at=1:length(yLabels), labels=yLabels, cex.axis=0.7)
#
# }
#
