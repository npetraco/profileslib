pca.fv<-function(fpath, keep.dim, chop.at.row) {
  ptr<-file(fpath, "rb")
  header.info<-read.digital.surf.header(ptr)
  surface<-read.digital.surf.profiles(ptr,header.info)
  close(ptr)
  print("Surface width and height:")
  print(dim(surface))
  
  #Look at it's column-wise variance:
  print("Running PCA")
  pca.model<-prcomp(surface,center=TRUE,scale=TRUE)
  keep.dim.cumvar<-summary(pca.model)$importance[3,keep.dim]
  
  dim100<-which(summary(pca.model)$importance[3,]==1.0)[1] #Dimension where 100% variance first occurs
  print(paste("100% Variance starts at dim:",dim100))
  print(paste("Keeping only up to dim: ", keep.dim, " = ", 100*keep.dim.cumvar, "% variance", sep=""))
  Mdim<-keep.dim                                         
  Zscores<-predict(pca.model)[,1:Mdim]
  
  Zscores<-Zscores[1:chop.at.row,]
  fv<-as.numeric(Zscores) #This should stack columns into a feature vector
  
  return(fv)
  
}