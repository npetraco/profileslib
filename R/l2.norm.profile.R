#------------------------------------------------------
#L^2 norm a profile
#|z|_2 = sqrt(sum(z_i^2)) = z dot z
#z_n = z/|z|
#------------------------------------------------------
l2.norm.profile<-function(profile){
  
  #pick out the numbers, ignore any NAs
  number.idxs<-which(!is.na(profile)==TRUE)
  actual.profile<-profile[number.idxs]
  
  normed.actual.profile<-actual.profile*(1/sqrt(actual.profile%*%actual.profile))
  
  normed.profile<-rep(NA,length(profile))
  for(i in 1:length(number.idxs)) {
    normed.profile[number.idxs[i]] <- normed.actual.profile[i]
  }
#  sapply(1:length(number.idxs), function(x){normed.profile[number.idxs[x]] <- normed.actual.profile[x]})
  
  return(normed.profile)
}