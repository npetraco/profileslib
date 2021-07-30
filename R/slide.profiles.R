#library(manipulate)

#-----------------------------------------------
#Push signal left or right with NA padding
#-----------------------------------------------
rotateLR<-function(profiledata,shift){
  if(shift>=0) {
    shiftedpattern<-c(rep(NA,shift),profiledata)
    #shiftedpattern=ArrayPad[profiledata,{shift,0}];
  } else {
    shiftedpattern<-c(profiledata,rep(NA,abs(shift)))
  }

  return(shiftedpattern);

}

#-----------------------------------------------
#Plot balanced shifted signals
#shift>0 shifts right
#shift<0 shifts left
#-----------------------------------------------
slide.over<-function(ref,sld,shift){

  ymax<-max(na.omit(c(ref,sld)))
  ymin<-min(na.omit(c(ref,sld)))
  #print(c(ymin,ymax))

  plot(rotateLR(ref,-shift),col="blue",ylab="",typ="l",ylim=c(ymin,ymax))
  par(new=TRUE)
  plot(rotateLR(sld,shift),col="red",typ="l",ylab="",main="Ref=Blue, Sliding=Red",ylim=c(ymin,ymax))

}

#-----------------------------------------------
#Interactively slide two signals over each other
#-----------------------------------------------
islide<-function(ref,sld,max.shift){
  manipulate(
    slide.over(ref,sld+sep.shift,shift),
    shift = slider(-max.shift,max.shift,initial=0),
    sep.shift = slider(-1,1,initial=0,step=0.1,ticks=FALSE)
    )
}

#-----------------------------------------------
#Feed in a matrix of signals and select which to
#slide
#-----------------------------------------------
islide2<-function(sdat,max.shift,max.y.sep){
  manipulate(
    slide.over(sdat[sig1,],sdat[sig2,]+sep.shift,shift),
    sig1 = slider(1,nrow(sdat),label="Ref (blue)"),
    sig2 = slider(1,nrow(sdat),label="Slide (red)"),
    shift = slider(-max.shift,max.shift,initial=0),
    sep.shift = slider(-max.y.sep,max.y.sep,initial=0,step=0.1,ticks=FALSE)
  )
}

#-----------------------------------------------
#Interactively slide two signals over each other.
#Signals are differnt lengths however
#-----------------------------------------------
islide3<-function(ref,sld,max.shift){
  if(length(ref)>length(sld)) {
    ref2<-ref
    sld2<-c(sld,rep(NA,length(ref)-length(sld)))
    #print(length(ref2))
    #print(length(sld2))
  }
  if(length(ref)<length(sld)) {
    sld2<-sld
    ref2<-c(ref,rep(NA,length(sld)-length(ref)))
    #print(length(ref2))
    #print(length(sld2))
  }
  if(length(ref)==length(sld)) {
    ref2<-ref
    sld2<-sld
    #print(length(ref2))
    #print(length(sld2))
  }

  manipulate(
    slide.over(ref2,sld2+sep.shift,shift),
    shift = slider(-max.shift,max.shift,initial=0),
    sep.shift = slider(-1,1,initial=0,step=0.1,ticks=FALSE)
  )
}
