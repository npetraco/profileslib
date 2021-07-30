#-----------------------------------------------------------------
#Align any number of profiles making up a retangular surface
#
#ALL PROFILES ARE OF THE SAME LENGTH.
#All profiles are sent in as a rectangular array (n-profiles by p-zheights)
#NO NAs ALLOWED!
#-----------------------------------------------------------------
align.profiles.within.surface<-function(profile.group,refp.idx,lagmax,tolerance,printQ=FALSE)
{
  
if(NA %in% as.matrix(profile.group))
 {
  print("No NAs allowed")
  print("Either NAs within profile (PROBLEM!) or NAs on edges")
  print("If NAs on edges try align.within.profile algorithm instead")
  return(0)
 }
 
#Number of profiles:
nump<-dim(profile.group)[1]   
#print(paste("Number of profiles: ",nump))

if(printQ==TRUE)
 {
  print("**********************************************")
  print(paste("Reference profile for this surface is #",refp.idx,sep=""))
  print(paste("Surface is ",nump," profiles by ",dim(profile.group)[2]," points",sep=""))
  print("**********************************************")
 }

refp<-profile.group[refp.idx,]
sftmat<-NULL
for(i in 1:nump)
 {
  if(i!=refp.idx)
   {
    #The shifting profile:   
    movp<-profile.group[i,]

    #Compute the required shift:
    lag.info<-maxcorr.lag.two.profiles(refp,movp,lagmax,tolerance)
    maxlag<-lag.info[3]
    
    if(lag.info[1]==0)
     {
      print(paste("CAUTION. Max.Corr. between refenence and profile ",i," is less than tolerance",sep=""))
      print(lag.info)
     }
    
    if(printQ==TRUE)
     {
      print(paste("Profiles #",i," vs. #",refp.idx,"(ref.)",sep=""))     
      print(paste("Maximum correlation at lag: ",maxlag,sep=""))
      print(paste("Max correlation value: ",lag.info[2],sep=""))
     }
    
    #Make up an NA vector for the shifting profile so it can be aligned with the reference: 
    #lag.na.vec<-rep(NA,abs(maxlag))

    #shifts: prepend ref, append ref, prepend mov, append mov
    sftvec<-c(0,0,0,0)        
    if(sign(maxlag)==1) #Pushing mov profile forward
     {
      sftvec[2]<-maxlag #mov profile goes past ref => append lag to ref profile
      sftvec[3]<-maxlag #mov profile goes past ref => prepend lag to mov profile
     }
    if(sign(maxlag)==-1) #Pushing mov profile backward
     {
      sftvec[1]<-abs(maxlag)    #prepend lag to ref
      sftvec[4]<-abs(maxlag)    #append lag to mov
     }
     
    #Append shift sign to the shift vec and put info in shift matrix 
    sftvec<-c(sign(maxlag),sftvec)    
    sftmat<-rbind(sftmat,sftvec)
    
   }
 }
#print("DIM SHIFT MAT:")
#print(dim(sftmat))
#print(count)
#print("")

if(printQ==TRUE)
 {
  colnames(sftmat)<-c("sgn lag|","ref pre|","ref app|","mov pre|","mov app")
  print(sftmat)
  print("Now aligning:....")
 }

#Now do shifting by prepending and appending appropriate amounts:
mp<-max(sftmat[,2]) #Max prepend length
ma<-max(sftmat[,3]) #Max append length

all.movp<-profile.group[-refp.idx,] #This drops the reference profile from the stack
#print("DIM all.movp:")
#print(dim(all.movp))

##Problem if sftmat is only one row. Then it is recognized as a numeric with only one dimension (length)
##Its most convenient to have two dimensions so change its class to a matrix:
if(class(all.movp)=="numeric")
 {
  print("Numeric!")
  all.movp<-t(as.matrix(all.movp)) #Change its class into a row matrix.
 }

#print(dim(sftmat)[1])
mod.profiles<-NULL
for(i in 1:dim(sftmat)[1])
 {
  if(sftmat[i,1]<0) #****If moving profile shifts left:
   {
   	#prepend and append appropriate amounts to mov profile wrt reference:
   	modmovp<-c(rep(NA,sftmat[i,4]),all.movp[i,],rep(NA,sftmat[i,5]))         
   	
   	#account for other moving profiles that may have shifted more:
   	modmovp<-c(rep(NA,length(refp)+mp-length(modmovp)),modmovp,rep(NA,ma))
   	 
   	mod.profiles<-rbind(mod.profiles,modmovp)
   }
  if(sftmat[i,1]>=0) #****If moving profile shifts right or doesn't need to shift:
   {
    #prepend and append appropriate amounts to mov profile wrt reference
    modmovp<-c(rep(NA,sftmat[i,4]),all.movp[i,],rep(NA,sftmat[i,5]))

    #account for other moving profiles that may have shifted more:
    modmovp<-c(rep(NA,mp),modmovp,rep(NA,(length(refp)+ma-length(modmovp))))

    mod.profiles<-rbind(mod.profiles,modmovp)
   }

 }

#Pad the reference:
refp<-c(rep(NA,mp),refp,rep(NA,ma))

#print(dim(mod.profiles))

#Reinsert reference profile into the surface
#reference was the first profile:
if(refp.idx==1)
 {
  #print("IN 1")	
  mod.profiles<-rbind(refp,mod.profiles)	
 }
#reference was the last profile: 
if(refp.idx==dim(profile.group)[1])
 {
  #print("IN 2")	
  mod.profiles<-rbind(mod.profiles,refp)	
 }
#reference was somewhere in between the top and bottom:
if((refp.idx!=1)&(refp.idx!=dim(profile.group)[1]))
 {
  #print("IN 3")
  mod.profiles<-rbind(mod.profiles[1:refp.idx-1,],refp,mod.profiles[refp.idx:dim(mod.profiles)[1],])
 }

#print(dim(mod.profiles))

#Even out the profiles:
#Recipe: Transpose, chop out rows with NAs, transpose back
mod.profiles<-t(na.omit(data.frame(t(mod.profiles))))

#print(dim(mod.profiles))

#Fix up the row and column names
rownames(mod.profiles)<-NULL
colnames(mod.profiles)<-NULL


return(mod.profiles)
 
}
