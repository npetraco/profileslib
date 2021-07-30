#-------------------------------------------------------------------
#Arrange the x,y,z, coordinates of the surface
#-------------------------------------------------------------------
xyz.surface.coords<-function(zmat,header)
{
	
num.profiles<-dim(zmat)[1]
num.pts<-dim(zmat)[2]
#print(num.profiles)
#print(num.pts)

#Both these are in units of microns:
x.len<-(header$num.pts.line-1)*(header$x.inc)*1000
y.len<-(header$num.lines-1)*(header$y.inc)*1000

xgrid<-seq(0, x.len, 1000*header$x.inc)
ygrid<-seq(0, y.len, 1000*header$y.inc)
#print(length(xgrid))
#print(length(ygrid))
numpoints<-length(xgrid)*length(ygrid)
print(paste("Grid is",length(xgrid),"by",length(ygrid),"(",numpoints,"points",")"))


coordlist<-lapply(1:length(xgrid),function(i)
                       {
                        xtmp<-rep(xgrid[i],length(ygrid))
                        cltmp<-cbind(xtmp,ygrid)
                       })
coordlist<-do.call(rbind,coordlist)
#print(dim(coordlist)[1])
#print(length(unlist(zmat)))
#print(class(unlist(zmat)))
#print(dim(as.vector(zmat)))
coordlist<-cbind(coordlist,as.vector(zmat))
#print(dim(coordlist))

#This seems to be too many points. "kilo-mate" the points:
pt.select<-seq(0,dim(coordlist)[1],1000)
print("Skipping every 1000 points for speed")
#print(length(pt.select))


return(coordlist[pt.select,])

}