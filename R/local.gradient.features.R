#-----------------------------------------------------------------
#Get first derivative critical values from a splined signal
#
#Much taken from continuous.to.rectangular.profile finction
#
#-----------------------------------------------------------------
profile.crit.vals<-function(profil,tol)
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
  
  return(cbind(deriv.zeros,extrema.typ))
  
}


#----------
#Grab a number of gradient valued to the left and right of the critical points
#----------
get.gradients<-function(profil,crit.val.info,gradvecs.leng) {
  
  deriv.zeros<-crit.val.info[,1]
  extrema.typ<-crit.val.info[,2]
  
  #fit profile with splines and compute first derivative:
  #profile "function":
  psf<-splinefun(1:length(profil),profil)
  
  #First derivative:
  dpsf<-psf(1:length(profil),deriv=1)
  
  
  #print(deriv.zeros)
  #right.grads<-matrix(rep(0,length(deriv.zeros)*gradvecs.leng), nrow=length(deriv.zeros),ncol=gradvecs.leng)
  #left.grads<-matrix(rep(0,length(deriv.zeros)*gradvecs.leng), nrow=length(deriv.zeros),ncol=gradvecs.leng)
  #print(dim(right.grads))
  right.peak.grads<-NULL
  left.peak.grads<-NULL
  right.valley.grads<-NULL
  left.valley.grads<-NULL
  for(i in 1:length(deriv.zeros)) {
    didx<-deriv.zeros[i]
    if(extrema.typ[i]==1) {
      #A peak
      peak.grad.right<-rep(0,gradvecs.leng)
      peak.grad.left<-rep(0,gradvecs.leng)
      for(j in 1:gradvecs.leng){
        peak.grad.right[j]<-(dpsf[didx+j])
        peak.grad.left[j]<-(dpsf[didx-j])
      }
      right.peak.grads<-rbind(right.peak.grads,peak.grad.right)
      left.peak.grads<-rbind(left.peak.grads,peak.grad.left)
    }
    
    if(extrema.typ[i]==-1) {
      #A valley
      valley.grad.right<-rep(0,gradvecs.leng)
      valley.grad.left<-rep(0,gradvecs.leng)
      for(j in 1:gradvecs.leng){
        valley.grad.right[j]<-(dpsf[didx+j])
        valley.grad.left[j]<-(dpsf[didx-j])
      }
      right.valley.grads<-rbind(right.valley.grads,valley.grad.right)
      left.valley.grads<-rbind(left.valley.grads,valley.grad.left)
    }
    
  }
#   print(right.peak.grads)
#   print(left.peak.grads)
#   print(right.valley.grads)
#   print(left.valley.grads)
#     print(c(min(right.peak.grads),max(right.peak.grads)))
#     print(c(min(left.peak.grads),max(left.peak.grads)))
#     print(c(min(right.valley.grads),max(right.valley.grads)))
#     print(c(min(left.valley.grads),max(left.valley.grads)))
  rownames(left.peak.grads)<-NULL
  rownames(right.peak.grads)<-NULL
  rownames(left.valley.grads)<-NULL
  rownames(right.valley.grads)<-NULL
  
  grad.chuncks.info<-list(left.peak.grads,right.peak.grads,left.valley.grads,right.valley.grads)
  names(grad.chuncks.info)<-c("Peak left grads","Peak right grads","Valley left grads", "Valley right grads")
  
  return(grad.chuncks.info)
}

#--------------------
#Histogram features
#--------------------
# hist.features(grad.list.info) {
#   
#   peak.mat<-cbind(grad.list.info[[1]],grad.list.info[[2]])
#   valley.mat<-cbind(grad.list.info[[4]],grad.list.info[[4]])
# 
#   peakdensityfeatures<-NULL
#   for(i in 1:nrow(peak.mat)) {
#     rawfeature<-na.omit(peak.mat[i,])
#     normalizedfeature<-(mean(rawfeature)-rawfeature)/sd(rawfeature)
#     peakdensityfeature<-hist(normalizedfeature,breaks=seq(-3,3,0.5),plot=F)
#     peakdensityfeatures<-rbind(peakdensityfeatures,peakdensityfeature)
#   }
#   
#   valleydensityfeatures<-NULL
#   for(i in 1:nrow(valley.mat)) {
#     rawfeature<-na.omit(valley.mat[i,])
#     normalizedfeature<-(mean(rawfeature)-rawfeature)/sd(rawfeature)
#     valleydensityfeature<-hist(normalizedfeature,breaks=seq(-3,3,0.5),plot=F)
#     valleydensityfeatures<-rbind(valleydensityfeatures,valleydensityfeature)
#   }
#   
#   features<-list(peakdensityfeatures,valleydensityfeature)
#   
#   
#   
# }