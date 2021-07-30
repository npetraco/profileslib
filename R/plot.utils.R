#---------------------------------------------------
#2D plot (image) of a Digital Surf formated surface
#---------------------------------------------------
plot.digital.surf.file<-function(digital.surf.file.info, x.decimation.factor, y.decimation.factor,contour.percent=100,titleQ=FALSE) {

  head.info<-digital.surf.file.info[[1]]
  surf.mat<-digital.surf.file.info[[2]]
  
  dec.col.idxs<-seq(from=1,to=as.numeric(head.info["num.pts.line"]),by=floor(as.numeric(head.info["num.pts.line"])/x.decimation.factor))
  dec.row.idxs<-seq(from=1,to=as.numeric(head.info["num.lines"]),by=floor(as.numeric(head.info["num.lines"])/y.decimation.factor))
  
  decimated.surf.mat<-surf.mat[dec.row.idxs,dec.col.idxs]
  
  num.pts.per.line<-as.numeric(head.info["num.pts.line"])
  xinc<-as.numeric(head.info["x.inc"])*1000 #microns
  xaxis<-seq(from=0,to=xinc*(num.pts.per.line-1),by=xinc)
  dec.xaxis<-xaxis[dec.col.idxs]
  
  num.lines<-as.numeric(head.info["num.lines"])
  yinc<-as.numeric(head.info["y.inc"])*1000 #microns
  yaxis<-seq(from=0,to=yinc*(num.lines-1),by=yinc)
  dec.yaxis<-yaxis[dec.row.idxs]
  
  #image(t(decimated.surf.mat))
  image(dec.xaxis, dec.yaxis, t(decimated.surf.mat), col = terrain.colors(24), xlab="x (microns)", ylab="y (microns)")
  if(contour.percent<100) {
    contour(dec.xaxis, dec.yaxis, t(decimated.surf.mat), 
            levels = seq(from=floor(min(decimated.surf.mat)), to=floor(max(decimated.surf.mat)), 
            by=floor((contour.percent/100)*(max(decimated.surf.mat)-min(decimated.surf.mat)))), 
            add = TRUE, col = "peru")
  }
  if(titleQ==TRUE) {
    title(main=as.character(head.info["name.obj"]))    
  }
}


#---------------------------------------------------
#3D plot  of a Digital Surf formated surface
#---------------------------------------------------
plot3D.digital.surf.file<-function(digital.surf.file.info, x.decimation.factor, y.decimation.factor) {
  
  head.info<-digital.surf.file.info[[1]]
  surf.mat<-digital.surf.file.info[[2]]
  
  dec.col.idxs<-seq(from=1,to=as.numeric(head.info["num.pts.line"]),by=floor(as.numeric(head.info["num.pts.line"])/x.decimation.factor))
  dec.row.idxs<-seq(from=1,to=as.numeric(head.info["num.lines"]),by=floor(as.numeric(head.info["num.lines"])/y.decimation.factor))
  
  decimated.surf.mat<-surf.mat[dec.row.idxs,dec.col.idxs]
  
  num.pts.per.line<-as.numeric(head.info["num.pts.line"])
  xinc<-as.numeric(head.info["x.inc"])*1000 #microns
  xaxis<-seq(from=0,to=xinc*(num.pts.per.line-1),by=xinc)
  dec.xaxis<-xaxis[dec.col.idxs]
  
  num.lines<-as.numeric(head.info["num.lines"])
  yinc<-as.numeric(head.info["y.inc"])*1000 #microns
  yaxis<-seq(from=0,to=yinc*(num.lines-1),by=yinc)
  dec.yaxis<-yaxis[dec.row.idxs]
  
  coords<-cbind(expand.grid(X=dec.xaxis, Y=dec.yaxis), as.numeric(t(decimated.surf.mat)))
  #print(dim(coords))
  #print(coords[,3])
  #open3d()                                          #For RGL
  #plot3d(coords[,1],coords[,2],coords[,3],type="s",radius=2,col="red",xlab="x",ylab="y",zlab="z")

  x<-coords[,1]
  y<-coords[,2]
  z<-coords[,3]
  
  zlim<-range(z)
  zlen <- zlim[2] - zlim[1] + 1
  #colorlut <- terrain.colors(zlen) # height color lookup table
  #colorlut <- rainbow(zlen) # height color lookup table
  #col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
  
  s=interp(x,y,z)
  open3d()
  plot3d(x,y,z,type="s",radius=0.01,xlab="x",ylab="y",zlab="z")
  surface3d(s$x,s$y,s$z,color="red")
  
}