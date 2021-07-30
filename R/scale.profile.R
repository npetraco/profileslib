#################################################################
#Subtrace mean of profile from profile and divide by profiles
#standard deviation. Sometimes this is called z-scaling or
#standardizing
#################################################################
auto.scale.profile<-function(profile, typ="standard")
{
  
  z.std<-NULL
  
  if(typ=="standard") {
    zb<-mean(profile)
    zsd<-sd(profile)
    z.std<-(profile-zb)
    z.std<-z.std/zsd
  }
  if(typ=="robust") {
    zm<-median(profile)
    zmad<-mad(profile)
    z.std<-(profile-zm)
    z.std<-z.std/zmad
  }
  
  return(z.std)
}