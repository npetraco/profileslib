#---------------------------------------------------------------------
#Generate a structuring element for morphological filtering
#Can be a ball or a line
#---------------------------------------------------------------------
structuring.element<-function(type,leng,header.info)
{

xinc<-header.info$x.inc

#Ball:	
if(type=="ball")
 {
  print(paste("The ball radius is: ",leng))	
 	
  #The length is a ball radius:
  ball.radius<-leng
  
  #x-axis over diameter of ball:
  xr<-seq(from=-ball.radius,to=ball.radius,by=xinc)
  #print(xr^2)
  
  #Ball structuring element:
  elmt<-sqrt((ball.radius^2)-(xr^2))
  
  #Length of ball:
  num.pts<-floor(((length(elmt))-1)/2)
 }

#Line:
if(type=="line")
 {
  print(paste("The line length is: ",leng))	

  #The length is a line element length:
  line.len<-leng
  
  #x-axis over length of line:
  xr<-seq(from=-line.len,to=line.len,by=xinc)
  
  #Line structuring element:
  elmt<-rep(1,length(xr))
  
  #Length of line:
  num.pts<-floor(line.len/xinc)
  
 } 

struc.el.info<-list(elmt,num.pts)
names(struc.el.info)<-c("elmt","num.pts")

return(struc.el.info)
	
}


#---------------------------------------------------------------------
#Dilation operator
#---------------------------------------------------------------------
dilation<-function(profile,struc.elem,header.info)
{

#This is n in Muralikrishnan and Raja:
prof.len<-length(profile)

#make a ball:
#For now, assume units are mm
xinc<-header.info$x.inc
#print(xinc)

#x-axis over diameter of ball:
#xr<-seq(from=-ball.radius,to=ball.radius,by=xinc)
#print(xr^2)
#B<-sqrt((ball.radius^2)-(xr^2))
#print(B)
B<-struc.elem$elmt

#Num pts in one radius of ball:
#m<-floor((length(B))-1)/2
#print(m)
m<-struc.elem$num.pts

#Pad the profile on the ends to take ball radius into account:
#z1<-c(numeric(m),profile,numeric(m))
#print(z1)

#Process interior (between ends) region of profile:
c1mid<-sapply((m+1):(prof.len-m),
     function(j)
      {
       #print(paste(length(profile[(j-m):(j+m)])," ",length(B)))
       max(profile[(j-m):(j+m)]+B)
      }   )
      
#plot(1:length(profile),profile,typ="l")      
#par(new=T)      
#plot(1:length(c1),c1,typ="l",col="green")      

#Process the first m points:
c1left<-sapply(1:m,
         function(j)
          {
           #print(paste(length(profile[1:(j+m)])," ",length(B[(m-j+2):(2*m+1)])))
           max(profile[1:(j+m)] + B[(m-j+2):(2*m+1)])
          }   )

#plot(1:length(c1left),c1left,typ="l",col="green")      

#Process the last m points:
c1right<-sapply((prof.len-m+1):prof.len,
          function(j)
           {
           	#print(paste(length(profile[(j-m):prof.len])," ",length(B[1:(m+(prof.len-j)+1)])))
           	max(profile[(j-m):prof.len]+B[1:(m+(prof.len-j)+1)])
           }   )

c1<-c(c1left,c1mid,c1right)

#plot(1:length(profile),profile,typ="l")
#par(new=T) 
#plot(1:length(c1),c1,typ="l",col="green")      
#print(length(profile))
#print(length(c1))

return(c1)

}


#---------------------------------------------------------------------
#Erosion operator
#---------------------------------------------------------------------
erosion<-function(profile,struc.elem,header.info)
{

c2<-((-1)*dilation(-1*profile,struc.elem,header.info))

#plot(1:length(profile),profile,typ="l")
#par(new=T) 
#plot(1:length(c2),c2,typ="l",col="green")      
#print(length(profile))
#print(length(c2))

return(c2)

}


#---------------------------------------------------------------------
#Closing operator
#---------------------------------------------------------------------
closing<-function(profile,struc.elem,header.info)
{

c3<- erosion(dilation(profile,struc.elem,header.info),struc.elem,header.info)

#plot(1:length(profile),profile,typ="l")
#par(new=T) 
#plot(1:length(c3),c3,typ="l",col="green")      
#print(length(profile))
#print(length(c3))

return(c3)

}


#---------------------------------------------------------------------
#Opening operator
#---------------------------------------------------------------------
opening<-function(profile,struc.elem,header.info)
{

c4<- dilation(erosion(profile,struc.elem,header.info),struc.elem,header.info)

#plot(1:length(profile),profile,typ="l")
#par(new=T) 
#plot(1:length(c4),c4,typ="l",col="green")      
#print(length(profile))
#print(length(c4))

#return(c4)

}


#---------------------------------------------------------------------
#Multi-scale decomposition
#---------------------------------------------------------------------
multi.scale.decomposition<-function(profil,levs,elem.type,num.x.incs,header.info)
{

#First make sure computations will work with the input parameters:
#Structuring  element length:
elem.leng<-num.x.incs*header.info$x.inc

#Max allowed struct. elem. leng must be less than half the profile length:
max.elem.leng<-(length(profil)*header.info$x.inc)/2
requested.max.elem.leng<-elem.leng*10^(levs-1)
if(requested.max.elem.leng>max.elem.leng)
 {
  print("******Requested structuring element maximum length exceeds allowed for this algorithm")
  print("Readjust num.x.incs and levs")
  print(paste("Request max leng:",requested.max.elem.leng))
  print(paste("Allowed max leng:",max.elem.leng))
  print("NULL!!!!!!!!!!!")
  return(NULL)
 }

decomp.pros<-NULL
diff.pros<-NULL
approx.prev<-profil
#elem.prev<-structuring.element(elem.type,elem.leng,header.info)
for(i in 1:levs)
 {
  #Scale: uses an element 10x bigger than previous:
  elem.prev<-structuring.element(elem.type,(elem.leng*(10^(i-1))),header.info)

  #Compute the scale approx according to struc. elem size:
  approx<-opening(closing(approx.prev,elem.prev,header.info),elem.prev,header.info)
  
  #Comprte the difference between the sucessive scale approximations:
  diff<-approx.prev-approx
  
  #Store the results:
  decomp.pros<-rbind(decomp.pros,approx)
  diff.pros<-rbind(diff.pros,diff)
  print(paste("Approx Done: ",i))
  
  #Start the next scale analysis with the previous approximation:
  approx.prev<-approx

 	
 }

total.decomp<-list(decomp.pros,diff.pros)
names(total.decomp)<-c("approx.profils","diff.profils")

return(total.decomp)
}