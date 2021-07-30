splt<-function(x,vec,colr,yrng,xrng) {
  if(x>=0){
    nvec<-c(rep(NA,x),vec)
    n<-length(nvec)
    plot(1:n,nvec,typ="l",col=colr,ylim=yrng,xlim=xrng)  
  }
  if(x<0) {
    nvec<-c(vec,rep(NA,-x))
    n<-length(nvec)
    plot(1:n,nvec,typ="l",col=colr,ylim=yrng,xlim=xrng)
  }
}
stkp<-function(pm1,pm2,pm3,pm4,pm5,vec1,vec2,vec3,vec4,vec5,yrange,xrange) {
  splt(pm1,vec1,"red",yrange,xrange)
  par(new=T)
  splt(pm2,vec2,"orange",yrange,xrange)
  par(new=T)
  splt(pm3,vec3,"violet",yrange,xrange)
  par(new=T)
  splt(pm4,vec4,"green",yrange,xrange)
  par(new=T)
  splt(pm5,vec5,"blue",yrange,xrange)
}
