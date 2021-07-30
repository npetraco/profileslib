#--------------------------------------------
#Pad signal to length divisible by 2^Jmax
#Aleviates having to compute Jmax
#--------------------------------------------
pad.signal.for.dwt<-function(signal,pad.typ) {
  
  N<-length(signal)
  Jmax<-floor(logb(N,2))
  ext.signal.leng<-(ceiling(N/(2^Jmax)) * 2^Jmax)
  pad.leng<-(ext.signal.leng - N)
  signal.ext<-c(signal,padding(pad.typ,pad.leng))
  
  return(signal.ext)
  
}


#--------------------------------------------
#Pad signal to length divisible by 2^J where J is the
#max level specified by the user.
#--------------------------------------------
pad.signal.for.dwtJ<-function(signal,J,pad.typ) {
  
  N<-length(signal)
  ext.signal.leng<-(ceiling(N/(2^J)) * 2^J)
  pad.leng<-(ext.signal.leng - N)
  #print(pad.leng)
  signal.ext<-c(signal,padding(pad.typ,pad.leng))
  
  return(signal.ext)
  
}


#--------------------------------------------
#Pad or chop 2D signal to length divisible by 
#2^J where J is the
#max level specified by the user.
#--------------------------------------------
adjust.signal.for.dwt.2d.J<-function(signal,J,pad.typ) {
  
  
  #The extended signal will need this many rows and columns:
  ssize<-2^J
  
  #I have this many rows and cols currently:
  dims<-dim(signal)
  
  #I need this many more/less cols:
  pad.leng.cols <- (ssize - dims[2])
  
  #I need this many more/less rows:
  pad.leng.rows <- (ssize - dims[1])
  
  print(paste("Curr rows: ",dims[1]))
  print(paste("Curr cols: ",dims[2]))
  
  print(paste("Pad rows: ",pad.leng.rows))
  print(paste("Pad cols: ",pad.leng.cols))
  
  print(paste("Ext rows: ",dims[1] + pad.leng.rows))
  print(paste("Ext cols: ",dims[2] + pad.leng.cols))
  
  #Don't do anything:
  if(pad.leng.rows == 0 && pad.leng.cols == 0) {
    signal.ext<-signal
    return(signal.ext)
  }
  #Dont do anything to cols, pad rows:
  if(pad.leng.rows > 0 && pad.leng.cols == 0) {
    pad.mat<-array( padding(pad.typ,pad.leng.rows*dims[2]), c(pad.leng.rows,dims[2]))
    signal.ext<-rbind(signal,pad.mat)             
    return(signal.ext)                            
  }
  #Dont do anything to cols, chop rows:
  if(pad.leng.rows < 0 && pad.leng.cols == 0) {
    signal.ext <- signal[1:(dim[1]+pad.leng.rows), ]  
    return(signal.ext)
  }
  #Pad cols, don't do anything to rows:
  if(pad.leng.rows == 0 && pad.leng.cols > 0) {
    pad.mat <- array( padding(pad.typ, dims[1]*pad.leng.cols), c(dims[1],pad.leng.cols))
    signal.ext <- cbind(signal,pad.mat)
    return(signal.ext)
  }
  #Chop cols, don't do anything to rows:
  if(pad.leng.rows == 0 && pad.leng.cols < 0) {
    signal.ext <- signal[,1:(dims[2]+pad.leng.cols)]
    return(signal.ext)
  }
  #Pad cols, Pad rows
  if(pad.leng.rows > 0 && pad.leng.cols > 0) {
    #Pad the columns first:
    #print("I'm here")
    pad.mat.cols <- array( padding(pad.typ, dims[1]*pad.leng.cols), c(dims[1],pad.leng.cols)) 
    #print("Done 1")
    signal.ext <- cbind(signal,pad.mat.cols)
    #print("Done 2")
    #Now pad the rows:
    pad.mat.rows <- array( padding(pad.typ, pad.leng.rows*dim(signal.ext)[2]), c(pad.leng.rows,dim(signal.ext)[2]))
    #print("Done 3")
    signal.ext <- rbind(signal.ext,pad.mat.rows)
    #print("Done 4")
    return(signal.ext)
  }
  #Pad cols, chop rows
  if(pad.leng.rows < 0 && pad.leng.cols > 0) {
    #Pad the columns first:
    pad.mat.cols <- array( padding(pad.typ, dims[1]*pad.leng.cols), c(dims[1],pad.leng.cols)) 
    #print(pad.mat.cols)
    signal.ext <- cbind(signal,pad.mat.cols)
    #Now chop the rows:
    signal.ext <- signal.ext[1:(dims[1]+pad.leng.rows),]
    return(signal.ext)
  }
  #Chop cols, pad rows
  if(pad.leng.rows > 0 && pad.leng.cols < 0) {
    #First chop the columns:
    signal.ext <- signal[,1:(dims[2]+pad.leng.cols)]
    #Now pad the rows:
    pad.mat.rows <- array( padding(pad.typ, pad.leng.rows*dim(signal.ext)[2]),c(pad.leng.rows,dim(signal.ext)[2]))
    signal.ext <- rbind(signal.ext, pad.mat.rows)
    return(signal.ext)
  }
  #Chop cols, chop rows
  if(pad.leng.rows < 0 && pad.leng.cols < 0) {
    signal.ext <- signal[,1:(dims[2]+pad.leng.cols)]
    signal.ext <- signal.ext[1:(dims[1]+pad.leng.rows),]
    return(signal.ext)
  }
}


#--------------------------------------------
#Wrapper to construct a dwt object
#NOTE: basis.choice can be a name for a waveslim 
#implemented basis or a list with high-pass and low-pass coefs.
#--------------------------------------------
construct.dwt.obj<-function(dwt.coef.list, basis.choice, boundary.typ) {
  dwt.obj <- dwt.coef.list
  names(dwt.obj) <- c(paste("d", 1:(length(dwt.obj)-1), sep=""), paste("s", (length(dwt.obj)-1), sep=""))
  
  class(dwt.obj) <- "dwt"
  attr(dwt.obj, "wavelet") <- basis.choice
  attr(dwt.obj, "boundary") <- boundary.typ
  
  return(dwt.obj)
}


#--------------------------------------------
#Zero all coefs except a level = level
#Useful for getting idwt at specified level
#NOTE: The s-level (course level) is max level + 1
#--------------------------------------------
zero.all.levels.except<-function(level,dwt.obj) {
  dwt.obj.zeroed<-dwt.obj
  for(i in 1:length(dwt.obj)) {
    if(i!=level) {
      dwt.obj.zeroed[[i]] <- numeric(length(dwt.obj.zeroed[[i]]))
    }
  }
  
  return(dwt.obj.zeroed)
}


#--------------------------------------------
#2D version of above.
#--------------------------------------------
zero.all.levels.except.2d<-function(level,dwt.obj) {
  dwt.obj.zeroed<-dwt.obj
  level.nms<-names(dwt.obj)
  
  #Call the coursest level one number higher. Makes seperation below easier. NEED TO LEARN REGEXPs!!!!!!!!!
  level.nms[length(level.nms)] <- paste("LL",as.character(as.numeric(gsub("\\D", "", level.nms[length(level.nms)]))+1),sep="")
  #print(level.nms[length(level.nms)])
  
  for(i in 1:length(dwt.obj)) {
    level.num<-as.numeric(gsub("\\D", "", level.nms[i]))
    if(level.num!=level) {
      dwt.obj.zeroed[[i]] <- array(0,dim(dwt.obj.zeroed[[i]]))
    }
  }
  
  return(dwt.obj.zeroed)
}


#--------------------------------------------
#Generalized 2D version of above.
#--------------------------------------------
general.zero.all.levels.except.2d<-function(levels,dwt.obj) {
  dwt.obj.zeroed<-dwt.obj
  level.nms<-names(dwt.obj)
  
  #Call the coursest level one number higher. Makes seperation below easier. NEED TO LEARN REGEXPs!!!!!!!!!
  level.nms[length(level.nms)] <- paste("LL",as.character(as.numeric(gsub("\\D", "", level.nms[length(level.nms)]))+1),sep="")
  #print(level.nms[length(level.nms)])
  
  for(i in 1:length(dwt.obj)) {
    level.num<-as.numeric(gsub("\\D", "", level.nms[i]))
    if(!is.element(level.num,levels)) { 
      dwt.obj.zeroed[[i]] <- array(0,dim(dwt.obj.zeroed[[i]]))
    }
  }
  
  return(dwt.obj.zeroed)
}

#--------------------------------------------
#2D Zero out specified sublevels
#Note: each level has three sublevels:
#LH, HL, HH
#The corsest level is LL only
#--------------------------------------------
zero.levels.2d<-function(sublevels.vec,dwt.obj) {
  dwt.obj.zeroed<-dwt.obj
  
  for(i in sublevels.vec) {
    dwt.obj.zeroed[[i]] <- array(0,dim(dwt.obj.zeroed[[i]]))
  }
  
  return(dwt.obj.zeroed)
}