#-----------------------------------------------------------------
#Find the lag (shift) that best aligns two profiles according to the
#ccf
#
#First profile is the reference profile
#Second profile is the shifted profile
#
#Function returns shift and an error code
#Error code is 0 if max correlation is < tolerance, 1 if all's well.
#-----------------------------------------------------------------
maxcorr.lag.two.profiles<-function(profile.ref,profile.mov,maximum.lag,tolerance)
{

ccflist<-ccf(profile.ref,profile.mov,lag.max=maximum.lag,na.action=na.omit,plot=FALSE,type="correlation")

#Maximum correlation of ccf within +/-maximum lag window:
maxcor<-max(ccflist$acf)

#Lag at which maximum of ccf occurs:
maxlag<-ccflist$lag[which(ccflist$acf==maxcor)]

#If maximum of ccf is leff than some tolerance indicate a problem
errcod<-1 #1 means all ok
if(maxcor<tolerance) errcod<-0
      
return(c(errcod,maxcor,maxlag))
	
}


#-----------------------------------------------------------------
#Find lags that give best ccf alignment wrt chosen reference
#-----------------------------------------------------------------
maxcorr.lag.all.profiles<-function(all.profiles,profile.ref,maximum.lag,tolerance)
{

num.profiles<-dim(all.profiles)[1]

#Drop the reference profile index from the list of profiles:
shifting.profiles<-1:numprofiles
shifting.profiles<-shifting.profiles[-profile.ref]

shiftmat<-t(sapply(shifting.profiles,function(x){align.two.profiles(all.profiles[profile.ref,],all.profiles[x,],maximum.lag,tolerance)}))

#put the refenence into the lag list at lag = 0
shiftmat<-insert(shiftmat,ats=profile.ref,values=0)

return(shiftmat)
	
}
