#-------------------------------------------------------------------
#pt = pointer returned by calling file on a digital surf .sur or 
#.pro binary format file
#
#**SEE FORMATS SECTION OF HELP IN MOUNTAINS SOFTWARE FOR 
#INFORMATION ON THE .sur and .pro formats
#
#-------------------------------------------------------------------
read.digital.surf.header<-function(pt)
{

#1
seek(pt, where = 0, rw="r")
code<-readBin(pt, what=character(), size = 12, n = 1, signed = FALSE, endian = "little")

#2
seek(pt, where = 12, rw="r")
format<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#3
seek(pt, where = 14, rw="r")
num.obj<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#4
seek(pt, where = 16, rw="r")
ver.num<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#5
seek(pt, where = 18, rw="r")
stud.typ<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#6********** 
seek(pt, where = 20, rw="r")
name.obj<-readBin(pt, what=character(), size = 30, n = 1, signed = FALSE, endian = "little")
name.obj<-strsplit(name.obj, " ")[[1]]
name.obj<-paste(name.obj[-c(which(name.obj==""),length(name.obj))],collapse="")
#print("name.obj")
#print(name.obj)
#print("")

#7**********
seek(pt, where = 50, rw="r")
name.op<-readBin(pt, what=character(), size = 30, n = 1, signed = FALSE, endian = "little")
#print("name.op")
#print(name.op)
#print("")

#8
seek(pt, where = 80, rw="r")
unused.int1<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#9
seek(pt, where = 82, rw="r")
unused.int2<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#10
seek(pt, where = 84, rw="r")
unused.int3<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#11
seek(pt, where = 86, rw="r")
non.meas.pts.flg<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#12
seek(pt, where = 88, rw="r")
abs.z.flg<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#13************
#This should be a raw byte sequence but is not reading in as such
#Will treat as an 8 byte integer
seek(pt, where = 90, rw="r")
reserv<-readBin(pt, what=integer(), size = 8, n = 1, signed = TRUE, endian = "little")
#print("reserv:")
#print(reserv)
#print("")

#14 Number of bits per point
seek(pt, where = 98, rw="r")
num.bits.pt<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
if(num.bits.pt==32)
 {
  print("WARNING. 32-bit FILE INDICATED!")	
 }
#print(num.bits.pt.a)
#num.bits.pt<-16 #CAUTION!!!!! PATCH, BECAUSE MOUNTAINS SOMETIMES SAVING AS 32!!!!!!!!!!! Why???????

#15
seek(pt, where = 100, rw="r")
min.z<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")

#16
seek(pt, where = 104, rw="r")
max.z<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")

#17 NUMBER OF POINTS PER PROFILE
seek(pt, where = 108, rw="r")
num.pts.line<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")

#18 NUMBER OF PROFILES IN SURFACE
seek(pt, where = 112, rw="r")
num.lines<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")

#19
seek(pt, where = 116, rw="r")
num.pts.total<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")

#20
seek(pt, where = 120, rw="r")
x.inc<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#21
seek(pt, where = 124, rw="r")
y.inc<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#22
seek(pt, where = 128, rw="r")
z.inc<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#23 **************
seek(pt, where = 132, rw="r")
x.label<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
x.label<-strsplit(x.label, " ")[[1]][1]
#print("x.label:")
#print(x.label)
#print("")

#24 ****************
seek(pt, where = 148, rw="r")
y.label<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
y.label<-strsplit(y.label, " ")[[1]][1]
#print("y.label:")
#print(y.label)
#print("")

#25 *****************
seek(pt, where = 164, rw="r")
z.label<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
z.label<-strsplit(z.label, " ")[[1]][1]
#print("z.label")
#print(z.label)
#print("")

#26
seek(pt, where = 180, rw="r")
x.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
x.unit<-strsplit(x.unit, " ")[[1]][1]

#27
seek(pt, where = 196, rw="r")
y.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
y.unit<-strsplit(y.unit, " ")[[1]][1]

#28
seek(pt, where = 212, rw="r")
z.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
z.unit<-strsplit(z.unit, " ")[[1]][1]

#29 ************
seek(pt, where = 228, rw="r")
x.len.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
x.len.unit<-strsplit(x.len.unit, " ")[[1]][1]
#print("x.len.unit:")
#print(x.len.unit)
#print("")

#30 ********
seek(pt, where = 244, rw="r")
y.len.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
y.len.unit<-strsplit(y.len.unit, " ")[[1]][1]
#print("y.len.unit")
#print(y.len.unit)
#print("")

#31 ***********
seek(pt, where = 260, rw="r")
z.len.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
z.len.unit<-strsplit(z.len.unit, " ")[[1]][1]
#print("z.len.unit")
#print(z.len.unit)
#print("")

#32
seek(pt, where = 276, rw="r")
x.ratio<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#33
seek(pt, where = 280, rw="r")
y.ratio<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#34
seek(pt, where = 284, rw="r")
z.ratio<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#35
seek(pt, where = 288, rw="r")
replica<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#36
seek(pt, where = 290, rw="r")
inverted<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#37
seek(pt, where = 292, rw="r")
leveled<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#38
seek(pt, where = 294, rw="r")
unused.sing1<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#39
seek(pt, where = 298, rw="r")
unused.sing2<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#40
seek(pt, where = 302, rw="r")
unused.sing3<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#41
seek(pt, where = 306, rw="r")
meas.time.sec<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#42
seek(pt, where = 308, rw="r")
meas.time.min<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#43
seek(pt, where = 310, rw="r")
meas.time.hr<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#44
seek(pt, where = 312, rw="r")
meas.time.day<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#45
seek(pt, where = 314, rw="r")
meas.time.month<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#46
seek(pt, where = 316, rw="r")
meas.time.yr<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#47
seek(pt, where = 318, rw="r")
meas.time.wkday<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#48
seek(pt, where = 320, rw="r")
meas.dur<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#49
seek(pt, where = 324, rw="r")
unused.int4<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#50
seek(pt, where = 326, rw="r")
unused.int5<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#51
seek(pt, where = 328, rw="r")
unused.sing4<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#52
seek(pt, where = 332, rw="r")
unused.int6<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#53
seek(pt, where = 334, rw="r")
leng.comment.zone<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#54
seek(pt, where = 336, rw="r")
leng.private.zone<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")

#55 *****************
#This probably should be a 128 byte raw as well, but will be treated as a 128 char
#Only type that seems to be able to be this big
seek(pt, where = 338, rw="r")
free.zone<-readBin(pt, what=character(), size = 128, n = 1, signed = FALSE, endian = "little")
#print("free.zone")
#print(free.zone)
#print("")

#56
seek(pt, where = 460, rw="r")
x.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#57
seek(pt, where = 464, rw="r")
y.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#58
seek(pt, where = 468, rw="r")
z.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#59
seek(pt, where = 472, rw="r")
TT.spacing<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#60 NOTE: documented same as above. Weird.
seek(pt, where = 472, rw="r")
TT.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")

#61 **********
seek(pt, where = 476, rw="r")
TT.name<-readBin(pt, what=character(), size = 13, n = 1, signed = FALSE, endian = "little")
#print("T.name")
#print(T.name)
#print("")

#62 ***********
seek(pt, where = 489, rw="r")
TT.step.unit<-readBin(pt, what=character(), size = 13, n = 1, signed = FALSE, endian = "little")
#print("T.step.unit")
#print(T.step.unit)
#print("")

#63 *********
#Comment zone. Assume a char for now.......
seek(pt, where = 489+13, rw="r")
comment.zone<-readBin(pt, what=character(), size = leng.comment.zone, n = 1, signed = FALSE, endian = "little")
#print("comment.zone:")
#print(comment.zone)
#print("")

#64 ************
#Private zone. Also assume a char for now.......
seek(pt, where = 489+13+leng.comment.zone, rw="r")
private.zone<-readBin(pt, what=character(), size = leng.private.zone, n = 1, signed = FALSE, endian = "little")
#print("private.zone")
#print(private.zone)
#print("")

#Short "essential" version:
#header.info<-list(name.obj, stud.typ,num.obj, num.bits.pt, num.pts.line, 
#                  num.lines, num.pts.total, max.z, min.z, x.inc, 
#                  y.inc, z.inc, x.off, y.off, 
#                  z.off, x.unit, y.unit, z.unit, x.len.unit, 
#                  y.len.unit, z.len.unit, leng.comment.zone, 
#                  leng.private.zone)
#names(header.info)<-c("name.obj","stud.typ","num.obj","num.bits.pt", 
#                      "num.pts.line", "num.lines", "num.pts.total", 
#                      "max.z", "min.z", "x.inc", "y.inc", "z.inc", 
#                      "x.off", "y.off", "z.off", "x.unit", "y.unit", 
#                      "z.unit", "x.len.unit", "y.len.unit", 
#                      "z.len.unit", "leng.comment.zone",
#                      "leng.private.zone")

#Long "total info" version:
header.info<-list(code,format,num.obj,ver.num,stud.typ,name.obj,
       name.op,unused.int1,unused.int2,unused.int3,
       non.meas.pts.flg,abs.z.flg,reserv,num.bits.pt,min.z,
       max.z,num.pts.line,num.lines,num.pts.total,x.inc,y.inc,
       z.inc,x.label,y.label,z.label,x.unit,y.unit,z.unit,
       x.len.unit,y.len.unit,z.len.unit,x.ratio,y.ratio,z.ratio,
       replica,inverted,leveled,unused.sing1,unused.sing2,
       unused.sing3,meas.time.sec,meas.time.min,meas.time.hr,
       meas.time.day,meas.time.month,meas.time.yr,meas.time.wkday,
       meas.dur,unused.int4,unused.int5,unused.sing4,unused.int6,
       leng.comment.zone,leng.private.zone,free.zone,x.off,y.off,
       z.off,TT.spacing,TT.off,TT.name,TT.step.unit,comment.zone,
       private.zone)

names(header.info)<-c("code","format","num.obj","ver.num",
     "stud.typ","name.obj",
     "name.op","unused.int1","unused.int2","unused.int3",
     "non.meas.pts.flg","abs.z.flg","reserv","num.bits.pt",
     "min.z","max.z","num.pts.line","num.lines",
     "num.pts.total","x.inc","y.inc","z.inc","x.label",
     "y.label","z.label","x.unit","y.unit","z.unit",
     "x.len.unit","y.len.unit","z.len.unit","x.ratio",
     "y.ratio","z.ratio","replica","inverted","leveled",
     "unused.sing1","unused.sing2","unused.sing3",
     "meas.time.sec","meas.time.min","meas.time.hr",
     "meas.time.day","meas.time.month","meas.time.yr",
     "meas.time.wkday","meas.dur","unused.int4",
     "unused.int5","unused.sing4","unused.int6",
     "leng.comment.zone","leng.private.zone","free.zone",
     "x.off","y.off","z.off","T.spacing","T.off",
     "T.name","T.step.unit","comment.zone","private.zone")
                      
return(header.info)
	
}

#-------------------------------------------------------------------
#Read in the profiles once the header info is known
#Fastest version
#-------------------------------------------------------------------
read.digital.surf.profiles<-function(pt,header.info)
{

#Offset for profiles:
offset<-512 + header.info$leng.comment.zone + header.info$leng.private.zone

#Move the pointer to where the profiles begin:
seek(pt, where = offset, rw="r")

#Read out the profiles. Units should be un microns:
#print(header.info$num.bits.pt)
point.byte.depth<-header.info$num.bits.pt/8
#print(point.byte.depth)
surface3<-readBin(pt, what=integer(), size = point.byte.depth, n = ((header.info$num.pts.line)*(header.info$num.lines)), signed = TRUE, endian = "little")
surface3<-(1000 * header.info$z.inc * t(matrix(surface3, ncol=header.info$num.lines, nrow=header.info$num.pts.line)))

return(surface3)
	
}


#-------------------------------------------------------------------
#Read in a Digital Surf file wrapper. Just specify the path to the
#.sur file. Use fle.choose() to open a file chooser window 
#-------------------------------------------------------------------
read.digital.surf.file<-function(file.path) {
  
  ptr<-file(file.path, "rb") #Open up a connection to the .sur file
  header.info<-read.digital.surf.header(ptr)
  surface.matrix<-read.digital.surf.profiles(ptr,header.info)
  close(ptr)
  
  all.file.info<-c(list(NULL),list(NULL))
  names(all.file.info)<-c("Surface Info","Surface")
  all.file.info[[1]]<-header.info
  all.file.info[[2]]<-surface.matrix
  
  return(all.file.info)
  
}


#-------------------------------------------------------------------
#Faster read in a Digital Surf file wrapper. Calls a faster C++ routine
#via Rcpp. Just specify the path to the
#.sur file. Use fle.choose() to open a file chooser window 
#-------------------------------------------------------------------
read.digital.surf.file2<-function(file.path) {

  print("=============================================")
  print("Working on file:")
  print(file.path)
  
  surface.info<-readSurFile(file.path)
  print(surface.info[[1]]$num.lines)
  #print(surface.info$num.pts.line)
  
  #Unwrap the list of lists output by the above Rcpp function.
  header.info<-rep(list(NULL),63)
  count<-1
  name.list<-NULL
  for(i in 1:4){
    head.list<-surface.info[[i]]
    name.list<-c(name.list,names(head.list))
    for(j in 1:length(head.list)) {
      header.info[[count]] <- head.list[[j]]
      count <- count + 1
    }
  }
  names(header.info)<-name.list
  
  surface.matrix<-header.info[["surface"]]
  
  print(header.info[["num.lines"]])
  
  surface.matrix<-t(matrix(surface.matrix,ncol=header.info[["num.lines"]],nrow=header.info[["num.pts.line"]]))
  all.file.info<-c(list(NULL),list(NULL))
  names(all.file.info)<-c("Surface Info","Surface")
  all.file.info[[1]]<-header.info[-63]
  all.file.info[[2]]<-surface.matrix
  
  return(all.file.info)
  
}


#-------------------------------------------------------------------
#Read in the the raw signed 16/32-bit surface. The header info is known 
#and passed into the function. NOTE: the signed 16/32-bit ints are 
#interally cast (by R) from bytes to readable format.
#NOTE: The Zeiss CSM-700 only recrds signed 16-bit z-heights. SOMETIMES
#Digital Surf's Mountains casts them to 32-bit signed ints when it thinks
#it needs more bit depth when saving to a .sur file.
#-------------------------------------------------------------------
read.digital.surf.rawsurface<-function(pt,header.info)
{
  
  #Offset for profiles:
  offset<-512 + header.info$leng.comment.zone + header.info$leng.private.zone
  
  #Move the pointer to where the profiles begin:
  seek(pt, where = offset, rw="r")
  
  #Read out the raw surface z-heights. Internally they are signed 16/32-bit ints:
  point.byte.depth<-header.info$num.bits.pt/8
  rawsurface<-readBin(pt, what=integer(), size = point.byte.depth, n = ((header.info$num.pts.line)*(header.info$num.lines)), signed = TRUE, endian = "little")
  rawsurface<-(t(matrix(rawsurface, ncol=header.info$num.lines, nrow=header.info$num.pts.line)))
  
  return(rawsurface)
  
}  


#-------------------------------------------------------------------
#Write a matrix of surface data (z-heights) to binary file.
#
#This procedure writes profile data with the digital surf
#header.
#
#This file should be openable not only by R and Mathematica
#but also by Mountains as long as the extension is .sur
#
#This function is most useful for operating on a surface in R
#and importing into Mathematica or MATLAB for further processing
#or Mountains for vewing. 
#
#digital.surf.header should be header info ORIGINALLY read in from
#Mountains so as not to change important information such as the
#microscopes z-increment information 
#-------------------------------------------------------------------

write.binary.digital.surf<-function(surf.mat,digital.surf.header,pt)
{

#print(digital.surf.header)

#1
#seek(pt, where = 0, rw="r")
#code<-readBin(pt, what=character(), size = 12, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$code,pt,size=12)

#2
#seek(pt, where = 12, rw="r")
#format<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$format,pt,size=2)

#3
#seek(pt, where = 14, rw="r")
#num.obj<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$num.obj,pt,size=2)

#4
#seek(pt, where = 16, rw="r")
#ver.num<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$ver.num,pt,size=2)

#5
#seek(pt, where = 18, rw="r")
#stud.typ<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$stud.typ,pt,size=2)


#6********** 
#seek(pt, where = 20, rw="r")
#name.obj<-readBin(pt, what=character(), size = 30, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$name.obj,pt,size=30)


#7**********
#seek(pt, where = 50, rw="r")
#name.op<-readBin(pt, what=character(), size = 30, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$name.op,pt,size=30)


#8
#seek(pt, where = 80, rw="r")
#unused.int1<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.int1,pt,size=2)


#9
#seek(pt, where = 82, rw="r")
#unused.int2<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.int2,pt,size=2)

#10
#seek(pt, where = 84, rw="r")
#unused.int3<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.int3,pt,size=2)

#11
#seek(pt, where = 86, rw="r")
#non.meas.pts.flg<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$non.meas.pts.flg,pt,size=2)

#12
#seek(pt, where = 88, rw="r")
#abs.z.flg<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$abs.z.flg,pt,size=2)

#13************
#This should be a raw byte sequence but is not reading in as such
#Will treat as an 8 byte integer
#seek(pt, where = 90, rw="r")
#reserv<-readBin(pt, what=integer(), size = 8, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$reserv,pt,size=8)

#14 Number of bits per point
#seek(pt, where = 98, rw="r")
#num.bits.pt<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$num.bits.pt,pt,size=2)

#15
#seek(pt, where = 100, rw="r")
#min.z<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$min.z,pt,size=4)

#16
#seek(pt, where = 104, rw="r")
#max.z<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$max.z,pt,size=4)

#17 NUMBER OF POINTS PER PROFILE
#seek(pt, where = 108, rw="r")
#num.pts.line<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$num.pts.line,pt,size=4)

#18 NUMBER OF PROFILES IN SURFACE
#seek(pt, where = 112, rw="r")
#num.lines<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$num.lines,pt,size=4)

#19
#seek(pt, where = 116, rw="r")
#num.pts.total<-readBin(pt, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$num.pts.total,pt,size=4)

#20
#seek(pt, where = 120, rw="r")
#x.inc<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$x.inc,pt,size=4)

#21
#seek(pt, where = 124, rw="r")
#y.inc<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$y.inc,pt,size=4)

#22
#seek(pt, where = 128, rw="r")
#z.inc<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$z.inc,pt,size=4)

#23 **************
#seek(pt, where = 132, rw="r")
#x.label<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$x.label,pt,size=16)

#24 ****************
#seek(pt, where = 148, rw="r")
#y.label<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$y.label,pt,size=16)

#25 *****************
#seek(pt, where = 164, rw="r")
#z.label<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$z.label,pt,size=16)

#26
#seek(pt, where = 180, rw="r")
#x.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$x.unit,pt,size=16)

#27
#seek(pt, where = 196, rw="r")
#y.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$y.unit,pt,size=16)

#28
#seek(pt, where = 212, rw="r")
#z.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$z.unit,pt,size=16)

#29 ************
#seek(pt, where = 228, rw="r")
#x.len.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$x.len.unit,pt,size=16)

#30 ********
#seek(pt, where = 244, rw="r")
#y.len.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$y.len.unit,pt,size=16)

#31 ***********
#seek(pt, where = 260, rw="r")
#z.len.unit<-readBin(pt, what=character(), size = 16, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$z.len.unit,pt,size=16)

#32
#seek(pt, where = 276, rw="r")
#x.ratio<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$x.ratio,pt,size=4)

#33
#seek(pt, where = 280, rw="r")
#y.ratio<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$y.ratio,pt,size=4)

#34
#seek(pt, where = 284, rw="r")
#z.ratio<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$z.ratio,pt,size=4)

#35
#seek(pt, where = 288, rw="r")
#replica<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$replica,pt,size=2)

#36
#seek(pt, where = 290, rw="r")
#inverted<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$inverted,pt,size=2)

#37
#seek(pt, where = 292, rw="r")
#leveled<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$leveled,pt,size=2)

#38
#seek(pt, where = 294, rw="r")
#unused.sing1<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.sing1,pt,size=4)

#39
#seek(pt, where = 298, rw="r")
#unused.sing2<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.sing2,pt,size=4)

#40
#seek(pt, where = 302, rw="r")
#unused.sing3<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.sing3,pt,size=4)

#41
#seek(pt, where = 306, rw="r")
#meas.time.sec<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.time.sec,pt,size=2)

#42
#seek(pt, where = 308, rw="r")
#meas.time.min<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.time.min,pt,size=2)

#43
#seek(pt, where = 310, rw="r")
#meas.time.hr<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.time.hr,pt,size=2)

#44
#seek(pt, where = 312, rw="r")
#meas.time.day<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.time.day,pt,size=2)

#45
#seek(pt, where = 314, rw="r")
#meas.time.month<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.time.month,pt,size=2)

#46
#seek(pt, where = 316, rw="r")
#meas.time.yr<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.time.yr,pt,size=2)

#47
#seek(pt, where = 318, rw="r")
#meas.time.wkday<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.time.wkday,pt,size=2)

#48
#seek(pt, where = 320, rw="r")
#meas.dur<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$meas.dur,pt,size=4)

#49
#seek(pt, where = 324, rw="r")
#unused.int4<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.int4,pt,size=2)

#50
#seek(pt, where = 326, rw="r")
#unused.int5<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.int5,pt,size=2)

#51
#seek(pt, where = 328, rw="r")
#unused.sing4<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.sing4,pt,size=4)

#52
#seek(pt, where = 332, rw="r")
#unused.int6<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$unused.int6,pt,size=2)

#53
#seek(pt, where = 334, rw="r")
#leng.comment.zone<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$leng.comment.zone,pt,size=2)

#54
#seek(pt, where = 336, rw="r")
#leng.private.zone<-readBin(pt, what=integer(), size = 2, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$leng.private.zone,pt,size=2)

#55 *****************
#This probably should be a 128 byte raw as well, but will be treated as a 128 char
#Only type that seems to be able to be this big
#seek(pt, where = 338, rw="r")
#free.zone<-readBin(pt, what=character(), size = 128, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$free.zone,pt,size=128)

#56
#seek(pt, where = 460, rw="r")
#x.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$x.off,pt,size=4)

#57
#seek(pt, where = 464, rw="r")
#y.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$y.off,pt,size=4)

#58
#seek(pt, where = 468, rw="r")
#z.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$z.off,pt,size=4)

#59
#seek(pt, where = 472, rw="r")
#T.spacing<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$T.spacing,pt,size=4)

#60 NOTE: documented same as above. Weird.
#seek(pt, where = 472, rw="r")
#T.off<-readBin(pt, what=single(), size = 4, n = 1, signed = TRUE, endian = "little")
writeBin(digital.surf.header$T.off,pt,size=4)

#61 **********
#seek(pt, where = 476, rw="r")
#T.name<-readBin(pt, what=character(), size = 13, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$T.name,pt,size=13)

#62 ***********
#seek(pt, where = 489, rw="r")
#T.step.unit<-readBin(pt, what=character(), size = 13, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$T.step.unit,pt,size=13)

#63 *********
#Comment zone. Assume a char for now.......
#seek(pt, where = 489+13, rw="r")
#comment.zone<-readBin(pt, what=character(), size = leng.comment.zone, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$comment.zone,pt,size=digital.surf.header$leng.comment.zone)

#64 ************
#Private zone. Also assume a char for now.......
#seek(pt, where = 489+13+leng.comment.zone, rw="r")
#private.zone<-readBin(pt, what=character(), size = leng.private.zone, n = 1, signed = FALSE, endian = "little")
writeBin(digital.surf.header$private.zone,pt,size=digital.surf.header$leng.private.zone)


#Read out the profiles. Units should be un microns:
#surface3<-readBin(pt, what=integer(), size = 2, n = ((header.info$num.pts.line)*(header.info$num.lines)), signed = TRUE, endian = "little")
#surface3<-(1000 * header.info$z.inc * t(matrix(surface3, ncol=header.info$num.lines, nrow=header.info$num.pts.line)))


surf.vec<-as.numeric(matrix(t((1/(1000*digital.surf.header$z.inc))*surf.mat),ncol=digital.surf.header$num.pts.total,nrow=1))
#writeBin(surf.vec,pt,size=(digital.surf.header$num.pts.total*integer()),signed=FALSE)
writeBin(surf.vec,pt,size=2*integer())
#sapply(1:length(surf.vec),function(x){writeBin(surf.vec[x],pt,size=integer())})

#size=(num.pts.stored*double())

}
