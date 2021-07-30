#-----------------------------------------------------------------
#Try to turn the profile into something that looks like a bar
#code
#
#-----------------------------------------------------------------
continuous.to.rectangular.profile<-function(profil,tol,plotQ=FALSE)
{

#fit profile with splines and compute first derivative:
#profile "function":
psf<-splinefun(1:length(profil),profil)

#First derivative:
dpsf<-psf(1:length(profil),deriv=1)

#Find zero crossings of derivative:
deriv.zeros<-NULL
deriv.vals<-NULL
for(i in 1:(length(dpsf)-1))
 {
  if(dpsf[i]*dpsf[i+1]<0)
   {   	
    deriv.zeros<-c(deriv.zeros,i)
   }
 }

#Drop out the ends if they had zero derivatives:
if(1 %in% deriv.zeros)
 {
  deriv.zeros<-deriv.zeros[-which(deriv.zeros==1)]
 }
if(length(profil) %in% deriv.zeros)
 {
  deriv.zeros<-deriv.zeros[-which(deriv.zeros==length(profil))]
 }

#Pluck out any small extrema. They are probably noise from the 
#differencing procedure above, or due to a long flat set of 
#adjacent extrema points.
extrema.height.diffs<-NULL 
for(i in 2:length(deriv.zeros))
 {
  #Subtract Right (current) extremum from extremum to the Left
  val<-profil[deriv.zeros[i]]-profil[deriv.zeros[i-1]]
  extrema.height.diffs<-c(extrema.height.diffs,val)
 }
 
#Get rid of "noise" extrema: 
col1<-which(abs(extrema.height.diffs)>tol)

#Keep the first (left most) extremum:
col1<-col1+1 #Shift the indices by 1
col1<-c(1,col1) #Tack on index for left most extremum

#Extrema (df/dx~0) indices:
deriv.zeros<-deriv.zeros[col1]

#extrema.height.differences missing difference for first extremum.
#To the left of first extremum may be a noise extremum (e.g.profile point 1).
#To get a more definitive idea of whether it is a min or a max, subtract
#it from the extremum to its right
first.extremum.height.diff<-val<-profil[deriv.zeros[1]]-profil[deriv.zeros[2]]
extrema.height.diffs<-c(first.extremum.height.diff,extrema.height.diffs)
extrema.height.diffs<-extrema.height.diffs[col1]

#Max=1,Min=-1. THERE SHOULD BE NO ZEROS!!!!!!!
extrema.typ<-sign(extrema.height.diffs)
if(0 %in% extrema.typ)
 {
  print("Undefined max/min!")
  print(cbind(col1,extrema.height.diffs,deriv.zeros,extrema.typ))
  return(0)	
 } 

#Loop over all but the last extremum:
barcode<-NULL
for(i in 1:(length(deriv.zeros)-1))
 {
  if(extrema.typ[i]!=extrema.typ[i+1])
   {
    peak.valley.span<-(deriv.zeros[i+1]-deriv.zeros[i])
    #print(paste("Half PV or VP span:",peak.valley.span/2))
    motif<-c(rep(extrema.typ[i],floor(peak.valley.span/2)),rep(extrema.typ[i+1],ceiling(peak.valley.span/2)))
   }
  if(extrema.typ[i]==extrema.typ[i+1])
   {
    pp.or.vv.span<-(deriv.zeros[i+1]-deriv.zeros[i])
    #print(paste("*PP or VV span:",pp.or.vv.span))
    motif<-rep(extrema.typ[i],pp.or.vv.span)
   }
  barcode<-c(barcode,motif)
 }

#Now take care of the ends, ie left of first extremum and right of last extremum
#Distance from beginning of profile to first extremum:
start.to.first.span<-deriv.zeros[1]

#Distance from last extremum to end of profile:
last.to.end.span<-(length(profil)-deriv.zeros[length(deriv.zeros)])

#print(start.to.first.span)
#print(last.to.end.span)

start.motif<-c(rep(0,floor(start.to.first.span/2)),rep(extrema.typ[1],ceiling(start.to.first.span/2)) )
end.motif<-c(rep(extrema.typ[length(extrema.typ)],floor(last.to.end.span/2)), rep(0,ceiling(last.to.end.span/2)))
barcode<-c(start.motif,barcode,end.motif)

#print(length(profil))
#print(length(barcode))
#plot(1:length(barcode),barcode,type="s")

if(plotQ==TRUE)
 {
  ymax<-max(profil)
  ymin<-min(profil)
  plot(1:length(profil),profil,typ="l",xlim=c(1,length(profil)),ylim=c(ymin,ymax))
  par(new=T)
  plot(deriv.zeros,rep(0,length(deriv.zeros)),col="green",xlim=c(1,length(profil)),ylim=c(ymin,ymax))
  par(new=T)
  plot(deriv.zeros,profil[deriv.zeros],col="red",xlim=c(1,length(profil)),ylim=c(ymin,ymax))
  par(new=T)
  plot(1:length(barcode),barcode,type="s")	
 }

return(barcode)	
	
}

#-----------------------------------------------------------------
#Turn entire surface into a "barcode surface" or a "mode barcode"
#
#output can be a surface representation of a profile
#Control with output.typ="profile" (default) or 
#output.typ="surface".
#"profile" type is constructed by taking the MODE of each column
#of the barcode "surface"
#-----------------------------------------------------------------
surface.barcode.representation<-function(surface.mat,tolerance,output.typ="profile")
{

num.profiles<-dim(surface.mat)[1]

barcode.surface<-t(sapply(1:num.profiles,function(x){continuous.to.rectangular.profile(surface.mat[x,],tolerance,plotQ=FALSE)}))

if(output.typ=="surface")
 {
  return(barcode.surface)	
 }

if(output.typ=="profile")
 {
  #This computes the MODE down each column of the barcode surface:
  barcode.profile<-sapply(1:ncol(barcode.surface),function(x){as.numeric(names(sort(-table(barcode.surface[,x])))[1])       })
  return(barcode.profile)
 }
 
#as.numeric(names(sort(-table(surf.dat2[,800])))[1])
	
}


#----------------------------------------------------------------
#Turn measurments of position and width of striation 
#lines into barcode.
#
#Works well with calibrated ImageJ line measurements
#----------------------------------------------------------------
measured.barcode.representation<-function(file.path,pattern.width,bin.width)
{

raw.dat<-read.csv(file.path,header=TRUE)
#print(raw.dat)

#Get the measurment data from imageJ output. These are the
#distances of the left and right edge of each line, 
#from the left edge of the striation pattern. 
dat<-raw.dat[,7]
#print(dat)

#Compute width and position of each "line"
pwmat<-NULL
loop.seq<-seq(from=1,to=length(dat),by=2)
#print(loop.seq)
for(i in loop.seq)
 {
  #print(dat[i+1]-dat[i])
  #print(i)
  pwmat<-rbind(pwmat,c(dat[i],(dat[i+1]-dat[i])))
 }
#print(pwmat)

#Compute the number of bins required for the input length
#of the striation pattern and input resolution between 
#points on the pattern (i.e. the sampling interval, 2*Nyquist freq??) 
nbins<-ceiling(pattern.width/bin.width)

#Zero out a barcode feature vector
barcode<-rep(0,nbins)
#print(barcode)

lr.bin.mat<-NULL
for(i in 1:nrow(pwmat))
 {
  left.bin.num<-floor(pwmat[i,1]/bin.width)
  #This is for a groove starting at length=0. Needed otherwise alg looks for bin 0!
  if(left.bin.num==0)
   {
    left.bin.num<-1
   }

  right.bin.num<-floor((pwmat[i,1] + pwmat[i,2])/bin.width)
  #This is for a groove ending at length=0. Needed otherwise alg looks for bin 0!
  if(right.bin.num==0)
   {
    right.bin.num<-1
   }

   lr.bin.mat<-rbind(lr.bin.mat,c(left.bin.num,right.bin.num))

   groove.range<-seq(from=left.bin.num,to=right.bin.num,by=1) #NOTE: the groove range itself may make a good feature vector. Perhaps for use with dtw or needleman-wunsuch
   #print(groove.range)
   tmp<-barcode
   tmp[groove.range]<-1
   barcode<-tmp
   #print(barcode) 
 }
#print(cbind(pwmat,lr.bin.mat))
#plot(barcode,type="l")

return(barcode)

}