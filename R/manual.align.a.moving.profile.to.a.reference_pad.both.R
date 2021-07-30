#################################################################
#MANUALLY Align arbitrary mov to reference, padding both
#################################################################
manual.align.mov.to.ref.pad.both<-function(ref,mov,specified.lag,padtyp,printQ=FALSE)
{

refp<-na.omit(ref)
movp<-na.omit(mov)
#croscor<-ccf(refp,movp,lag.max=lagmax,plot=FALSE,na.action=na.omit) #Using Cross-Correlation for alignment
#refp<-ref
#movp<-mov
  
#lags<-croscor$lag
#max.corr.idx<-which(croscor$acf[,,1]==max(croscor$acf[,,1]))
#maxlag<-lags[max.corr.idx]
maxlag<-specified.lag
dif<-length(refp)-length(movp)

if(printQ==TRUE)
 {
  print(paste("Length of Reference: ",length(refp),sep=""))
  print(paste("  Length of Profile: ",length(movp),sep=""))
  print(paste("  Length difference: ",dif,sep=""))
  #print(paste("  Index of max corr: ",max.corr.idx,sep=""))
  print(paste("  SPECIFIED lag: ",maxlag,sep=""))
  #print(paste("           Max Corr: ",max(croscor$acf[,,1]),sep=""))
 }

#*********sftvec format c(ref shorter/longer -1/1, shift right/left 1/-1, num.pad.left, num.ref.pad.left, num.pad.right, num.ref.pad.right)
    
if(dif<=0) #Case: ref shorter than mov
 {
  #Subcase 1.1
  #shift mov forward (get overhang on right, need padding for mov on left and padding for ref on right)
  if(maxlag>=0) 
   {
    #initialize just in case:
    num.pad.left<-0
    num.ref.pad.right<-0
    
    #determine params:
    num.pad.left<-maxlag
    num.ref.pad.right<-(abs(dif) + maxlag)
    if(printQ==TRUE)
     print(paste("Ref: SHORTER, 1.1. Shift mov ---->. PAD mov LEFT:",num.pad.left,"PAD ref RIGHT:", num.ref.pad.right))
    sftvec<-c(-1,1,num.pad.left,0,0,num.ref.pad.right)
      
    #carry out action on profiles
    sftd.movp<-c(padding(padtyp,num.pad.left), movp)      #pad mov left    #PPPPPPmmmmmmmmmmmmmmmm
    sftd.refp<-c(refp,padding(padtyp,num.ref.pad.right))  #pad ref right   #rrrrrrrrrrrPPPPPPPPPPP
   }
  
  #Subcase 1.2a
  # shift move backward, but get overhang on both right and left
  if(maxlag<0 & ( abs(maxlag) <= abs(dif) )  ) 
   {
    #initialize just in case:
    num.ref.pad.left<-0
    num.ref.pad.right<-0
    
    #determine params:
    num.ref.pad.left<-abs(maxlag)
    num.ref.pad.right<-( abs(dif)-abs(maxlag) )
    if(printQ==TRUE)
     print(paste("Ref: SHORTER 1.2a. Shift mov <----. PAD ref LEFT:",num.ref.pad.left,"PAD ref RIGHT:", num.ref.pad.right))
    sftvec<-c(-1,-1,0, num.ref.pad.left,0, num.ref.pad.right)
    
    #carry out action on ref:                                          
    sftd.movp<-movp                                                                                                    #mmmmmmmmmmmmmmmmmmmmmm
    sftd.refp<-c(padding(padtyp,num.ref.pad.left),refp,padding(padtyp,num.ref.pad.right))  #pad ref left and right    #PPPrrrrrrrrPPPPPPPPPPP
    
   }
  
  #Subcase 1.2b
  # shift move backward, get overhang on left but need padding for mov on right, ref on left
  if(maxlag<0 & ( abs(maxlag) > abs(dif) ) )
   {
    #initialize just in case:
    num.ref.pad.left<-0
    num.pad.right<-0
    
    #determine params:
    num.ref.pad.left<-abs(maxlag)
    num.pad.right<- ( abs(maxlag) - abs(dif) )
    if(printQ==TRUE)
     print(paste("Ref: SHORTER. 1.2b. Shift mov <----. PAD ref LEFT:",num.ref.pad.left,"PAD mov RIGHT:", num.pad.right))
    sftvec<-c(-1,-1,0, num.ref.pad.left, num.pad.right,0)
    
    #carry out action on profiles:                                                
    sftd.refp<-c(padding(padtyp,num.ref.pad.left), refp)  #pad ref left   PPPPPrrrrrrrrrrr  
    sftd.movp<-c(movp, padding(padtyp,num.pad.right) )    #pad mov right  mmmmmmmmmmmmmmPP  
   }
 }

if(dif>0) #Case: ref longer than mov
 {
  #Subcase 2.1
  #shift mov backward (get overhang on left, need padding for mov on right and ref on left)
  if(maxlag<=0) 
   {
    #initialize just in case:
    num.ref.pad.left<-0
    num.pad.right<-0
    
    #determine params:
    num.ref.pad.left<-abs(maxlag)
    num.pad.right<-(dif + abs(maxlag))
    if(printQ==TRUE)
     print(paste("Ref: LONGER 2.1. Shift mov <----. PAD ref LEFT:",num.ref.pad.left,"PAD mov RIGHT:", num.pad.right))
    sftvec<-c(1,-1, 0, num.ref.pad.left, num.pad.right, 0)
    
    #carry out action on mov:                                             
    sftd.refp<-c(padding(padtyp,num.ref.pad.left), refp)    #pad ref left   PPPPPPPPrrrrrrrrrrrrrrrrr          
    sftd.movp<-c(movp, padding(padtyp,num.pad.right) ) #pad mov right  mmmmmmmmmmmmmmmmPPPPPPPPP
   }

  #Subcase 2.2a
  #shift mov forward, need padding on both left and right
  if(maxlag>0 & ( maxlag<=dif ) )
   {
    #initialize just incase:
    num.pad.left<-0
    num.pad.right<-0
    
    #determine params:
    num.pad.left<-maxlag
    num.pad.right<-(dif-maxlag)
    if(printQ==TRUE)
     print(paste("Ref: LONGER. 2.2a. Shift mov ---->. PAD mov LEFT:",num.pad.left,"PAD mov RIGHT:", num.pad.right))
    sftvec<-c(1, 1, num.pad.left, 0, num.pad.right, 0)
      
    #carry out action on mov:                                                                        rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
    sftd.refp<-refp
    sftd.movp<-c(padding(padtyp,num.pad.left), movp, padding(padtyp,num.pad.right) ) #pad left/right PPPPPPPPPPmmmmmmmmmmmmmmmmPPPPPPP
   }
   
  #Subcase 2.2b
  #shift mov forward, need padding on left of mov and padd overhang over mov on right or ref
  if(maxlag>0 & ( maxlag>dif ) )
   {
    #initialize just incase:
    num.pad.left<-0
    num.ref.pad.right<-0
    
    #determine params:
    num.pad.left<-maxlag
    num.ref.pad.right<-(maxlag - dif)
    if(printQ==TRUE)
     print(paste("Ref: LONGER. 2.2b. Shift mov ---->. PAD mov LEFT:",num.pad.left,"PAD ref RIGHT:", num.ref.pad.right))
    sftvec<-c(1, 1, num.pad.left, 0, 0, num.ref.pad.right)
    
    #carry out action on profiles:
    sftd.movp<-c(padding(padtyp,num.pad.left), movp)      #pad mov left   PPPPPPPPPPPmmmmmmmmmmmmm
    sftd.refp<-c(refp, padding(padtyp,num.ref.pad.right)) #pad ref right  rrrrrrrrrrrrrrrrrrrrrPPP
   }
 }

if(length(sftd.refp)!=length(sftd.movp))
 {
  print("********WARNING, WARNING, WARNING!!!!! Reference and shifted profiles are different lengths! They are not supposed to be!******")
  print("CHECK FOR AN ERROR IN THE CODE OR WITH THE PROFILES!")
 }

return(list(sftd.refp, sftd.movp, sftvec))

}