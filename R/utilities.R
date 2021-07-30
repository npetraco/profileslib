#-------------------------------------------------------------------
#Make up a lable vector according to number of samples
#-------------------------------------------------------------------
generate.label.vec<-function(num.samps.vec)
{

lbl.vec<-NULL
for(i in 1:length(num.samps.vec) )
 {
  grpids<-rep(i,num.samps.vec[i])
  lbl.vec<-factor(c(lbl.vec,grpids))
 }
 
return(lbl.vec)
	
}

#-------------------------------------------------------------------
#Counts the number of replicates for each group.
#Should be independent of group naming convention
#-------------------------------------------------------------------
count.group.replicates<-function(arb.lbls)
{

num.samps.vec<-sapply(as.numeric(levels(factor(arb.lbls))), function(x){sum(arb.lbls==x)})
return(num.samps.vec)

}

#-------------------------------------------------------------------
#Generate an indicator matrix from a num samps per group vec
#ie 1 1 1 2 2 =
#1 0 0 0 0
#1 0 0 0 0
#1 0 0 0 0
#0 1 0 0 0
#0 1 0 0 0
#-------------------------------------------------------------------
generate.indicator.mat<-function(num.obs.per.samp.vec)
{

k<-length(num.obs.per.samp.vec)

lbl.indic.mat<-NULL
for(i in 1:k)
 {
  rws<-num.obs.per.samp.vec[i]
  block<-matrix(numeric(rws*k),nrow=rws,ncol=k)
  block[,i]<-rep(1,rws)
  #print(block)
  lbl.indic.mat<-rbind(lbl.indic.mat,block)
 }

return(lbl.indic.mat)
  
}

#-------------------------------------------------------------------
#Pick out groups of observations and form a new X matrix 
#-------------------------------------------------------------------
pick.out.groups<-function(X.mat,all.lbls,grp.picks)
{

pick.out.rows<-NULL
new.grp.lbls<-NULL
for(i in 1:length(grp.picks))
 {
  grp.idxs<-which(as.numeric(all.lbls)==as.numeric(grp.picks[i]))
  pick.out.rows<-c(pick.out.rows,grp.idxs)
  #print(grp.idxs)
  new.grp.lbl<-rep(i,length(grp.idxs))
  new.grp.lbls<-c(new.grp.lbls,new.grp.lbl)	
 }

new.grp.lbls<-factor(new.grp.lbls)
new.X.mat<-X.mat[pick.out.rows,] 

return(list(new.X.mat,new.grp.lbls))
 
}


#--------------
#Pad function
#--------------
padding<-function(typ,leng)
{

if(typ[[1]]=="NAs")
 {
  pad<-rep(NA,leng)
  return(pad)
 }
  
if(typ[[1]]=="zeros")
 {
  pad<-numeric(leng)
  return(pad)
 }

if(typ[[1]]=="uniform")
 {
  pad<-runif(leng, min=(typ[[2]]), max=typ[[3]])
  return(pad)
 }
  
if(typ[[1]]=="normal")
 {
  pad<-rnorm(leng, mean=typ[[2]], sd=typ[[3]])
  return(pad)  
 }

if(typ[[1]]=="constant")
 {
  pad<-rep(typ[[2]],leng)
  return(pad)
 }


}

#----------------------------------------------------
#Replace NAs in a matrix with something else
#----------------------------------------------------
replace.NAs<-function(dmat,padtype)
{
#Logical matrix for NAs:  
na.mat<-is.na(dmat)

#Pick out the NA indices. Had to do it this way bec. which was not cooperating with NAs in dmat....
na.idxs<-which(na.mat==TRUE,arr.ind=TRUE)

#Cast into an indicator matrix
na.mat<-1*na.mat

dmat.new<-dmat

#Temporarily replace NAs with 0s:
dmat.new[na.idxs]<-0

#Make a padding matrix
pad.mat<-na.mat*matrix( padding(padtype, nrow(na.mat)*ncol(na.mat)), nrow=nrow(na.mat), ncol=ncol(na.mat))

#Replace the 0s with the padding:
dmat.new<-dmat.new+pad.mat

return(dmat.new)
}


#----------------------------------------------------
#Index to subscripts function for a rectangular matrix
#----------------------------------------------------
ind2sub<-function(idx,num.row)
{
# rc.ind <- c(row(M)[ind], col(M)[ind] )
r = ((idx-1) %% num.row) + 1
c = floor((idx-1) / num.row) + 1

subscrpt<-c(r,c)

return(subscrpt)
}

#----------------------------------------------------
#Subscripts to index function for a retangular matrix
#----------------------------------------------------
sub2ind<-function(r,c,num.cols)
{
  
idx<-((r-1)*num.cols + c)
 
return(idx)
}

#---------------------------------------------------------------------
#Subscripts to index function for a STRICTLY UPPER TRIANGULAR matrix
#From: http://xlinux.nist.gov/dads/HTML/strictlyUpperTriangMat.html
#Matrix is assumed square => num.rows == num.cols
#---------------------------------------------------------------------
sub2ind.sutri<-function(i,j,num.rows)
{
  
idx<-(2*num.rows - i) * ((i - 1)/2) - i + j
 
return(idx)
}


#---------------------------------------------------------------------
#Index to subscripts LOOK UP TABLE for a STRICTLY UPPER TRIANGULAR matrix
#This is a hack. Need a fomula so we can do this on the fly.
#Matrix is assumed square => num.rows == num.cols
#---------------------------------------------------------------------
ind2sub.sutri.lut<-function(num.rows)
{

lut<-NULL
for(i in 1:num.rows)
 {
  for(j in 1:num.rows)
   {
    if(j>i)
     {
      lut<-rbind(lut,c(i,j,sub2ind.sutri(i,j,num.rows)))
     }
   }
 }

return(lut)
}


#-------------------------------------------------------
#Peter Shenkin's upper (lower) linear-to-lower triangle index2subscripts
#-------------------------------------------------------
upper.ind2sub<-function(kcol,spos) {
  # Finding index (row,column) of element in triangular array
k <- kcol    # number of columns
s <- spos    # position of element
ms <- 1:k  
ctm <-((ms-1)/2)*(2*k-ms+2) # No. of elements until row before ms - ms vector
ctm
# ms2 <- k:1
# ctm2 <- cumsum(ms2)
# ctm2
difr <- s-ctm # Finding element row count up to s element look for negative
difr
difrct <- length(difr[difr>=0])   # Now have number of row for sth element
difrct
if(difr[difrct]==0) {rowd <- difrct-1; cold <- k-rowd+1} else {rowd <- difrct; cold <- s-((rowd-1)/2)*(2*k-rowd+2)}  # Special case for element at end of row
return(c(rowd,cold))
}


#-------------------------------------------------------------------
#Extract the lower triangle of elements in a square matrix 
#-------------------------------------------------------------------
extract.lower.triangle<-function(sqmat)
{

lt.vec<-NULL
for(i in 1:dim(sqmat)[1])
  {
   for(j in 1:dim(sqmat)[2])
     {
      if(i<j)
       {
         lt.vec<-c(lt.vec,sqmat[i,j])
       }	
     }	
  }

return(lt.vec)

}


#------------------------------------------------------------------
#Find right hand side breaks in contiguous runs of 0s in a vector 
#EG: 000000000xxxxxxxx00000000xxxx000000000000000
#    cccccccccd       ccccccccd   ccccccccccccccc
#             Left Break      Left Break
#
#EG: xx000000000xxxxxxxx00000000xxxx000000000000000xxxxx
#      cccccccccd       ccccccccd   cccccccccccccccd
#------------------------------------------------------------------
left.run.breaks<-function(avec) {

  zero.idxs<-which(avec==0)

  forward.conn.vec<-rep("X",length(zero.idxs))
  #forward.conn.vec[1]<-"d"             #d = a break
                                        #c = elem in a contiguous run
  for(i in 1:(length(zero.idxs)-1)) {
    if(zero.idxs[i]+1 == (zero.idxs[i+1])) {
      forward.conn.vec[i]="c"
    } else {
      forward.conn.vec[i]="d"
    }
  }
  left.break.idxs<-(zero.idxs[which(forward.conn.vec=="d")]+1)
  
  return(left.break.idxs) 
}

#--------------------------------------------------------------
#Upsample a vector by inserting num zeros between each element
#--------------------------------------------------------------
# up.sample<-function(avec,num) {
#   uped.avec<-rep(0,length(avec)+length(avec)*num)
#   j<-1
#   for(i in 1:length(avec)) {
#     uped.avec[j] <- avec[i]
#     j <- j+num+1
#   }
#   return(uped.avec)
# }


#------------------------------------------------------------------
#Find left hand side breaks in contiguous runs of 0s in a vector 
#EG: 000000000xxxxxxxx00000000xxxx000000000000000
#                    dcccccccc   dccccccccccccccc
#                    Right Break Right Break
#
#EG: xx000000000xxxxxxxx00000000xxxx000000000000000xxxxx
#     dccccccccc       dcccccccc   dccccccccccccccc
#------------------------------------------------------------------
right.run.breaks<-function(avec) {
  
  zero.idxs<-which(avec==0)
  
  back.conn.vec<-rep("XX",length(zero.idxs))
  #back.conn.vec[1]<-"d"
  #back.conn.vec[length(zero.idxs)]<-"d"
  for(i in 2:(length(zero.idxs))) {
    if(zero.idxs[i]-1 == (zero.idxs[i-1])) {
      back.conn.vec[i]="c"
    } else {
      back.conn.vec[i]="d"
    }
  }
  right.break.idxs<-(zero.idxs[which(back.conn.vec=="d")]-1)
  
  return(right.break.idxs) 
}

#------------------------------------------------------------------
#Returns spans between which are non-zero values
#------------------------------------------------------------------
#find.breaks<-function(avec) {
#  return(cbind(left.run.breaks(avec),right.run.breaks(avec)))
#}