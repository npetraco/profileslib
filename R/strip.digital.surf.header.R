#-------------------------------------------------------------------
#In a convoluted way, this strips off the 512+ byte header on
#the digital surf binary files.
#
#It works by reading each 2 or 4-byte chunk (z-height) of the 
#surface data and immediately writing it to the output file as 
#the same 2 or 4-byte chunk.
#
#The purpose of this procedure is to prep a surface or profile
#for compression  
#-------------------------------------------------------------------
strip.digital.surf.header<-function(pt.in,pt.out,header.info)
{
#Number of bytes per point (ie "chunk size"):
num.bytes.pt<-header.info$num.bits.pt/8

#Offset for surface/profile data:
offset<-512 + header.info$leng.comment.zone + header.info$leng.private.zone

#Total number of points ("chunks") to read:
num.chunks<-((header.info$num.lines) * header.info$num.pts.line)

#Move the pointer of the input file to where the profiles begin:
seek(pt.in, where = offset, rw="r")

#Read in a profile (a line) and write it out to file
sapply(1:header.info$num.lines,
 function(dummy)
  {
   dummy;
   val<-readBin(pt.in, what="int", size = 2, n = (2*header.info$num.pts.line), signed = TRUE, endian = "little")
   writeBin(val,pt.out, size=2, endian="little") 
  }   )

print("Done writting header stripped data")

}
