#-------------------------------------------------------------------
#Write profile to binary file.
#
#This procedure writes profile data WITHOUT the digital surf
#header.
#
#It is meant to simply store profile data in a binary file for
#later use
#
#-------------------------------------------------------------------
write.binary.profile<-function(profl,pt.out)
{

num.pts<-length(profl)
writeBin(profl,pt.out,size=(num.pts*double()),endian="little")
	
}


#-------------------------------------------------------------------
#Write a group of profiles 
#binary file.
#
#Note: all profiles are assumes to have the same length.
#
#-------------------------------------------------------------------
write.binary.profiles<-function(profile.dat,file.path)
{

#close(ptr)

num.pro<-nrow(profile.dat)
pro.lens<-rep.int(ncol(profile.dat),num.pro)
num.pts.stored<-sum(pro.lens)

#First open file to write:
ptr<-file(file.path, "wb")

#Store:
#num profiles, 4 byte int
writeBin(num.pro,ptr,size=integer(),endian="little")

#vector of profile lengths:
#Note: these are all the same for now: 
writeBin(pro.lens,ptr,size=(num.pro*integer()),endian="little")

#Put all profiles into one row vector
all.pts<-as.vector(matrix(t(profile.dat),nrow=1,ncol=num.pts.stored))

#Now write the profiles vector
writeBin(all.pts,ptr,size=(num.pts.stored*double()),endian="little")
close(ptr)

}


#-------------------------------------------------------------------
#Write a list of profiles 
#binary file.
#
#Profiles can all be different lengths
#
#ASSUMES ONLY ONE PROFILE IN EACH ELEMENT OF THE LIST!!!!!!
#-------------------------------------------------------------------
write.list.binary.profiles<-function(list.of.profiles,file.path)
{

pro.lengs<-sapply(1:length(list.of.profiles),function(x){length(list.of.profiles[[x]])})

num.pts.stored<-sum(pro.lengs)
num.pro<-length(pro.lengs)
all.pts<-unlist(list.of.profiles)
print(pro.lengs)
print(num.pts.stored)
print(num.pro)
print(length(all.pts))

#Store:
#num profiles, 4 byte int
#length each profile, each a 4 byte int
#profile points, doubles. offset should be 4*(1+num profiles) bytes
ptr<-file(file.path, "wb")
writeBin(num.pro,ptr,size=integer(),endian="little")
writeBin(pro.lengs,ptr,size=(num.pro*integer()),endian="little")
writeBin(all.pts,ptr,size=(num.pts.stored*double()),endian="little")
close(ptr)
  
}

#-------------------------------------------------------------------
#Write a list of profiles composed of blocks of profiles to
#binary file.
#
#Useful to write intra-group aligned profiles to file
#
#-------------------------------------------------------------------
write.block.list.binary.profiles<-function(list.of.profiles,file.path)
{
  pro.lengs<-NULL
  for(i in 1:length(list.of.profiles)) {
    
    pro.mat<-list.of.profiles[[i]]
    num.pro.in.block<-nrow(pro.mat)
    num.pts.in.block<-ncol(pro.mat)
    pro.lengs<-c(pro.lengs,rep(num.pts.in.block,num.pro.in.block))
  }
  
  num.pro<-length(pro.lengs)
  #all.pts<-unlist(list.of.profiles)
  num.pts.stored<-sum(pro.lengs)
  #print(num.pts.stored)
  #Store:
  #num profiles, 4 byte int
  #length each profile, each a 4 byte int
  #profile points, doubles. offset should be 4*(1+num profiles) bytes
  
  print(num.pro)
  print(pro.lengs)
  
  ptr<-file(file.path, "ab") #Write to append
  writeBin(num.pro,ptr,size=integer(),endian="little")
  writeBin(pro.lengs,ptr,size=(num.pro*integer()),endian="little")
  for(i in 1:length(list.of.profiles)) {
    
    pro.mat<-list.of.profiles[[i]]
    
    #Put all profiles into one row vector
    all.pts.profiles<-as.vector(matrix(t(pro.mat),nrow=1,ncol=nrow(pro.mat) * ncol(pro.mat)))
    
    #Now write the block vector
    writeBin(all.pts.profiles,ptr,size=(length(all.pts.profiles)*double()),endian="little")
    
  }
  
  #writeBin(all.pts,ptr,size=(num.pts.stored*double()),endian="little")
  close(ptr)
  
  
}
#-------------------------------------------------------------------
#Write a group of profiles to 
#64-bit binary file to be read by imageJ.
#
#Open the file in imageJ as a RAW image Real-64
#No header in the file is written but the dimensions of the surface 
#must be supplied to imageJ.
#
#Note: all profiles are assumes to have the same length!
#
#-------------------------------------------------------------------
write.imageJ.binary.profiles<-function(profile.dat,file.path)
{

#close(ptr)

num.pro<-nrow(profile.dat)
pro.lens<-rep.int(ncol(profile.dat),num.pro)
num.pts.stored<-sum(pro.lens)

#First open file to write:
ptr<-file(file.path, "wb")

#Store:
#num profiles, 4 byte int
#writeBin(num.pro,ptr,size=integer(),endian="little")

#vector of profile lengths:
#Note: these are all the same for now: 
#writeBin(pro.lens,ptr,size=(num.pro*integer()),endian="little")

#Put all profiles into one row vector
all.pts<-as.vector(matrix(t(profile.dat),nrow=1,ncol=num.pts.stored))

#Now write the profiles vector
writeBin(all.pts,ptr,size=(num.pts.stored*double()),endian="little")
close(ptr)

}


#-------------------------------------------------------------------
#Write selection of mean profiles to binary file.
#
#-------------------------------------------------------------------
write.binary.mean.median.profiles<-function(root.input.path,root.output.path,name.vec,num.vec,typ,sftQ)
{

#close(ptr)

num.grps<-length(name.vec)

mean.pro.lens<-NULL
median.pro.lens<-NULL
#Binary write append mean and median profiles to these file paths:
outfil.mean<-paste(root.output.path,"BIN_FILES/",typ,"_means_unshifted_profiles.bin",sep="")
outfil.median<-paste(root.output.path,"BIN_FILES/",typ,"_medians_unshifted_profiles.bin",sep="")

for(i in 1:num.grps)
 {
  for(j in 1:num.vec[i])
   {
    #infil<-paste(root.path,"SUR_FILES/",name.vec[i],"_SUR_FILES/",name.vec[i],"_",as.character(j),"_",typ,".sur",sep="")
    #infil<-paste(root.path,"SUR_FILES/",name.vec[i],"/",name.vec[i],"_",as.character(j),"_",typ,".sur",sep="")
    infil<-paste(root.input.path,"SUR_FILES/",name.vec[i],as.character(j),"_",typ,".sur",sep="")
    ptr<-file(infil, "rb")
    header.info<-read.digital.surf.header(ptr)
    #hedinfo<-header.info[[1]]
    #print(header.info[[14]])
    #print(header.info)
    surface<-read.digital.surf.profiles(ptr,header.info)
    close(ptr)
    #print(dim(surface))
    
    #Generate mean and median of UNSHIFTED profiles:
    mean.pro<-colMeans(surface)
    #print(mean.pro)
    median.pro<-apply(surface,2,median)
    
    #Record number of points in profiles.
    #Need this info to open the binary file full of processed profiles:
    #NOTE: length(pros) == header.info$num.pts.line
    mean.pro.lens<-c(mean.pro.lens,length(mean.pro))
    median.pro.lens<-c(median.pro.lens,length(median.pro))
    
    #Do means:
    ptr<-file(outfil.mean, "ab")
    write.binary.profile(mean.pro,ptr)
    close(ptr)
    #print("")

    #Do medians:
    ptr<-file(outfil.median, "ab")
    write.binary.profile(median.pro,ptr)
    close(ptr)
    
    print(infil)
    print(paste("Wrote group ", i, "profile ", j))
    
   }	
 }

#result.file.header<-list(mean.pro.lens,median.pro.lens)
#names(result.file.header)<-c("mean.profile.lengths","median.profile.lengths")
num.pts.stored.mean<-sum(mean.pro.lens)
num.pts.stored.median<-sum(median.pro.lens)
num.mean.pro<-length(mean.pro.lens)
num.median.pro<-length(mean.pro.lens)

ptr<-file(outfil.mean, "rb")
all.mean.pts<-readBin(ptr, what=double(), size=8, n=num.pts.stored.mean, signed = TRUE, endian = "little")
close(ptr)

ptr<-file(outfil.median, "rb")
all.median.pts<-readBin(ptr, what=double(), size=8, n=num.pts.stored.median, signed = TRUE, endian = "little")
close(ptr)

#Store:
#num profiles, 4 byte int
#length each profile, each a 4 byte int
#profile points, doubles. offset should be 4*(1+num profiles) bytes
#First reopen file to overwrite it
ptr<-file(outfil.mean, "wb")
writeBin(num.mean.pro,ptr,size=integer(),endian="little")
writeBin(mean.pro.lens,ptr,size=(num.mean.pro*integer()),endian="little")
writeBin(all.mean.pts,ptr,size=(num.pts.stored.mean*double()),endian="little")
close(ptr)

ptr<-file(outfil.median, "wb")
writeBin(num.median.pro,ptr,size=integer(),endian="little")
writeBin(median.pro.lens,ptr,size=(num.median.pro*integer()),endian="little")
writeBin(all.median.pts,ptr,size=(num.pts.stored.median*double()),endian="little")
close(ptr)
	
}

#-------------------------------------------------------------------
#Slight modification of Write selection of mean profiles to binary file.
#
#-------------------------------------------------------------------
write.binary.mean.median.profiles.mod<-function(root.path,name.vec,num.vec)
{

#close(ptr)

num.grps<-length(name.vec)

mean.pro.lens<-NULL
median.pro.lens<-NULL
#Binary write append mean and median profiles to these file paths:
outfil.mean<-paste(root.path,"BIN_FILES/","means_unshifted_profiles.bin",sep="")
outfil.median<-paste(root.path,"BIN_FILES/","medians_unshifted_profiles.bin",sep="")

#print(outfil.mean)
#print(outfil.median)

for(i in 1:num.grps)
 {
  for(j in 1:num.vec[i])
   {
    infil<-paste(root.path,"SUR_FILES/",name.vec[i],"-",as.character(j),sep="")
    ptr<-file(infil, "rb")
    header.info<-read.digital.surf.header(ptr)
    #print(header.info)
    surface<-read.digital.surf.profiles(ptr,header.info)
    close(ptr)
    #print(dim(surface))
    
    #Generate mean and median of UNSHIFTED profiles:
    mean.pro<-colMeans(surface)
    median.pro<-apply(surface,2,median)
    
    #Record number of points in profiles.
    #Need this info to open the binary file full of processed profiles:
    #NOTE: length(pros) == header.info$num.pts.line
    mean.pro.lens<-c(mean.pro.lens,length(mean.pro))
    median.pro.lens<-c(median.pro.lens,length(median.pro))
    
    #Do means:
    ptr<-file(outfil.mean, "ab")
    write.binary.profile(mean.pro,ptr)
    close(ptr)
    #print("")

    #Do medians:
    ptr<-file(outfil.median, "ab")
    write.binary.profile(median.pro,ptr)
    close(ptr)
    
    print(infil)
    print(paste("Wrote group ", i, "profile ", j))
    
   }	
 }

#result.file.header<-list(mean.pro.lens,median.pro.lens)
#names(result.file.header)<-c("mean.profile.lengths","median.profile.lengths")
num.pts.stored.mean<-sum(mean.pro.lens)
num.pts.stored.median<-sum(median.pro.lens)
num.mean.pro<-length(mean.pro.lens)
num.median.pro<-length(mean.pro.lens)

ptr<-file(outfil.mean, "rb")
all.mean.pts<-readBin(ptr, what=double(), size=8, n=num.pts.stored.mean, signed = TRUE, endian = "little")
close(ptr)

ptr<-file(outfil.median, "rb")
all.median.pts<-readBin(ptr, what=double(), size=8, n=num.pts.stored.median, signed = TRUE, endian = "little")
close(ptr)

#Store:
#num profiles, 4 byte int
#length each profile, each a 4 byte int
#profile points, doubles. offset should be 4*(1+num profiles) bytes
#First reopen file to overwrite it
ptr<-file(outfil.mean, "wb")
writeBin(num.mean.pro,ptr,size=integer(),endian="little")
writeBin(mean.pro.lens,ptr,size=(num.mean.pro*integer()),endian="little")
writeBin(all.mean.pts,ptr,size=(num.pts.stored.mean*double()),endian="little")
close(ptr)

ptr<-file(outfil.median, "wb")
writeBin(num.median.pro,ptr,size=integer(),endian="little")
writeBin(median.pro.lens,ptr,size=(num.median.pro*integer()),endian="little")
writeBin(all.median.pts,ptr,size=(num.pts.stored.median*double()),endian="little")
close(ptr)
	
}


#-------------------------------------------------------------------
#Write selection of mean profiles to binary file.
#Special for Julie's naming convention.
#-------------------------------------------------------------------
write.binary.mean.median.profiles.julie<-function(root.input.path,root.output.path,name.vec,typ,sftQ)
{
  
  #close(ptr)
  
  num.grps<-length(name.vec)
  
  mean.pro.lens<-NULL
  median.pro.lens<-NULL
  #Binary write append mean and median profiles to these file paths:
  outfil.mean<-paste(root.output.path,"BIN_FILES/",typ,"_means_unshifted_profiles.bin",sep="")
  outfil.median<-paste(root.output.path,"BIN_FILES/",typ,"_medians_unshifted_profiles.bin",sep="")
  
  for(i in 1:num.grps)
  {
    
    
      infil<-paste(root.input.path,name.vec[i],"_",typ,".sur",sep="")
      #ptr<-file(infil, "rb")
      #header.info<-read.digital.surf.header(ptr)
      #hedinfo<-header.info[[1]]
      #print(header.info[[14]])
      #print(header.info)
      #surface<-read.digital.surf.profiles(ptr,header.info)
      info<-read.digital.surf.file2(infil)
      #close(ptr)
      header.info<-info[[1]]
      surface<-info[[2]]
      
      #print(dim(surface))
      
      #Generate mean and median of UNSHIFTED profiles:
      mean.pro<-colMeans(surface)
      #print(mean.pro)
      median.pro<-apply(surface,2,median)
      
      #Record number of points in profiles.
      #Need this info to open the binary file full of processed profiles:
      #NOTE: length(pros) == header.info$num.pts.line
      mean.pro.lens<-c(mean.pro.lens,length(mean.pro))
      median.pro.lens<-c(median.pro.lens,length(median.pro))
      
      #Do means:
      ptr<-file(outfil.mean, "ab")
      write.binary.profile(mean.pro,ptr)
      close(ptr)
      #print("")
      
      #Do medians:
      ptr<-file(outfil.median, "ab")
      write.binary.profile(median.pro,ptr)
      close(ptr)
      
      print(infil)
      print(paste("Wrote profile ", i))
      
      
  }
  
  #result.file.header<-list(mean.pro.lens,median.pro.lens)
  #names(result.file.header)<-c("mean.profile.lengths","median.profile.lengths")
  num.pts.stored.mean<-sum(mean.pro.lens)
  num.pts.stored.median<-sum(median.pro.lens)
  num.mean.pro<-length(mean.pro.lens)
  num.median.pro<-length(mean.pro.lens)
  
  ptr<-file(outfil.mean, "rb")
  all.mean.pts<-readBin(ptr, what=double(), size=8, n=num.pts.stored.mean, signed = TRUE, endian = "little")
  close(ptr)
  
  ptr<-file(outfil.median, "rb")
  all.median.pts<-readBin(ptr, what=double(), size=8, n=num.pts.stored.median, signed = TRUE, endian = "little")
  close(ptr)
  
  #Store:
  #num profiles, 4 byte int
  #length each profile, each a 4 byte int
  #profile points, doubles. offset should be 4*(1+num profiles) bytes
  #First reopen file to overwrite it
  ptr<-file(outfil.mean, "wb")
  writeBin(num.mean.pro,ptr,size=integer(),endian="little")
  writeBin(mean.pro.lens,ptr,size=(num.mean.pro*integer()),endian="little")
  writeBin(all.mean.pts,ptr,size=(num.pts.stored.mean*double()),endian="little")
  close(ptr)
  
  ptr<-file(outfil.median, "wb")
  writeBin(num.median.pro,ptr,size=integer(),endian="little")
  writeBin(median.pro.lens,ptr,size=(num.median.pro*integer()),endian="little")
  writeBin(all.median.pts,ptr,size=(num.pts.stored.median*double()),endian="little")
  close(ptr)
  
}

#-------------------------------------------------------------------
#Write a RANDOM selection of profiles from groups of surfaces to 
#binary file.
#
#The number of profiles to select from each surface is 1 for now.
#
#-------------------------------------------------------------------
write.random.binary.profiles<-function(root.path,name.vec,num.vec,typ,sftQ,eraseQ)
{

#close(ptr)

num.grps<-length(name.vec)

pro.lens<-NULL
#Binary write append selected profiles to this file path:
outfil<-paste(root.path,"BIN_FILES/",typ,"_RANDOM_unshifted_profiles.bin",sep="")

if(eraseQ==TRUE)
 {
  system(paste("rm -f",outfil)) 	
 }


for(i in 1:num.grps)
 {
  for(j in 1:num.vec[i])
   {
    infil<-paste(root.path,"SUR_FILES/",name.vec[i],"_SUR_FILES/",name.vec[i],"_",as.character(j),"_",typ,".sur",sep="")
    ptr<-file(infil, "rb")
    header.info<-read.digital.surf.header(ptr)
    #print(header.info)
    surface<-read.digital.surf.profiles(ptr,header.info)
    close(ptr)
    #print(dim(surface))
    
    #Sample an UNSHIFTED profile from a surface:
    choice<-sample(nrow(surface),size=1)
    chosen.pro<-surface[choice,]
        
    #Record number of points in profiles.
    #Need this info to open the binary file full of processed profiles:
    #NOTE: length(pros) == header.info$num.pts.line
    pro.lens<-c(pro.lens,length(chosen.pro))
    
    #Write append profile to output file:
    ptr<-file(outfil, "ab")
    write.binary.profile(chosen.pro,ptr)
    close(ptr)
    #print("")
    
    print(infil)
    print(paste("Wrote group ", i, "profile ", j, "profile#:", choice))
    
   }	
 }

#Read the chosen profiles back into memory:
#result.file.header<-list(pro.lens)
#names(result.file.header)<-c("profile.lengths")
num.pts.stored<-sum(pro.lens)
num.pro<-length(pro.lens)

ptr<-file(outfil, "rb")
all.pts<-readBin(ptr, what=double(), size=8, n=num.pts.stored, signed = TRUE, endian = "little")
close(ptr)

#Store:
#num profiles, 4 byte int
#length each profile, each a 4 byte int
#profile points, doubles. offset should be 4*(1+num profiles) bytes
#First reopen file to overwrite it
ptr<-file(outfil, "wb")
#Write the header info first:
writeBin(num.pro,ptr,size=integer(),endian="little")
writeBin(pro.lens,ptr,size=(num.pro*integer()),endian="little")
#Now write the profiles back:
writeBin(all.pts,ptr,size=(num.pts.stored*double()),endian="little")
close(ptr)

}

#-------------------------------------------------------------------
#Write N sets (num.samples) of RANDOMLY chosen profiles from a set of surfaces. 
#
#This is useful for preforming the same recogniton task on many randomly
#chosen profiles which can characterize a surface. 
#
#Given a set of n surfaces each with n_i profiles, a random profile (row)
#is drawn for each surface. This random profile is written to (binary)
#file. This is done N times producing N files each with n randomly 
#chosen profile, one from each surface.  
#
#CAUTION. Since files are written to disk, there could be a lot of
#disk consumed for large samples, say N > 200.
#-------------------------------------------------------------------
sample.random.profile.sets<-function(input.path, output.path, name.vec, num.vec, typ, num.samples)
{

num.grps<-length(name.vec)  

#Open output files
file.ptrs<-lapply(1:num.samples,function(x){file(paste(output.path,"randomly_selected_profiles_",as.character(x),".bin",sep=""),"ab")})
#print(file.ptrs)

#Open each surface file (one at a time) and select num.samples (ie N) random profiles from each surface
pro.lens<-NULL #will hold lengths of profiles selected from each surface
for(i in 1:num.grps) #Loop over the groups of surfaces:
 {
  for(j in 1:num.vec[i]) #Loop over the replicates in each group:
   {
    infil<-paste(input.path, name.vec[i],as.character(j),"_",typ,".sur",sep="")
    ptr<-file(infil, "rb")
    header.info<-read.digital.surf.header(ptr)
    #print(header.info)
    surface<-read.digital.surf.profiles(ptr,header.info)
    close(ptr)
    #print(dim(surface))

    #Sample (with replacement) the requested number of UNSHIFTED profiles from a surface:
    choice<-sample(nrow(surface),size=num.samples,replace=TRUE)
    chosen.profiles<-surface[choice,]

    #Record number of points in a set of profiles. 
    #Each set should have the same length becaues all the profiles are from the same surface
    #Need this info to open the binary file full of processed profiles:
    #NOTE: length(pros) == header.info$num.pts.line
    pro.lens<-c(pro.lens,dim(chosen.profiles)[2])
    #print(pro.lens)

    #Write append profile to output file:
    lapply(1:num.samples,function(x){write.binary.profile(chosen.profiles[x,],file.ptrs[[x]])})
    #lapply(1:num.samples,function(x){print(file.ptrs[[x]])})
    
    #print(infil)
    print(paste("Group", i, "replicate", j, ", selected profile#'s:"))
    print(choice)
   }  
 }

#Close the output files temporarily
junk<-lapply(1:num.samples,function(x){close(file.ptrs[[x]])})
    
#Read the chosen profiles back into memory:
#Necessary to attach the headers
num.pts.stored<-sum(pro.lens)
num.pro<-length(pro.lens)

#Read out the points stored for each file. Prep to add headers.
file.ptrs<-lapply(1:num.samples,function(x){file(paste(output.path,"randomly_selected_profiles_",as.character(x),".bin",sep=""),"rb")})
list.all.pts<-lapply(1:num.samples,function(x){NULL})
for(i in 1:num.samples)
 {
  all.pts<-readBin(file.ptrs[[i]], what=double(), size=8, n=num.pts.stored, signed = TRUE, endian = "little")
  list.all.pts[[i]]<-all.pts
 }
#Annoying, but close the output files temporarily AGAIN!
junk<-lapply(1:num.samples,function(x){close(file.ptrs[[x]])})

#Finally open up the output files for writing over, this time with the headers
file.ptrs<-lapply(1:num.samples,function(x){file(paste(output.path,"randomly_selected_profiles_",as.character(x),".bin",sep=""),"wb")})
for(i in 1:num.samples)
 {
  #all.pts<-0
  all.pts<-list.all.pts[[i]]
  
  #Write the header info first:
  writeBin(num.pro,file.ptrs[[i]],size=integer(),endian="little")
  writeBin(pro.lens,file.ptrs[[i]],size=(num.pro*integer()),endian="little")
  #Now write the profiles back:
  writeBin(all.pts,file.ptrs[[i]],size=(num.pts.stored*double()),endian="little")
  print(paste("**Wrote profile sample file:",i))
 }
#Close the output files, now with headers.
junk<-lapply(1:num.samples,function(x){close(file.ptrs[[x]])})

}


#-------------------------------------------------------------------
#Read in processed binary profiles
#-------------------------------------------------------------------
read.binary.profiles<-function(file.path)
{

#close(ptr)

#open the file:
ptr<-file(file.path, "rb")

#Read out the first number = #of profiles:
seek(ptr, where = 0, rw="r")
num.profiles<-readBin(ptr, what=integer(), size = 4, n = 1, signed = TRUE, endian = "little")
#print(num.profiles)

#Read out the int vector of profile lengths:
seek(ptr, where = 4, rw="r")
profile.lengths <- readBin(ptr, what=integer(), size = 4, n = num.profiles, signed = TRUE, endian = "little")
#print(profile.lengths)

#Read out the int vector of profile lengths:
#seek(ptr, where =(4 + 4*num.profiles), rw="r")

#Pointer should now be at the beginning of the profile data
max.len<-max(profile.lengths)

#print(num.profiles)
#print(profile.lengths)
#print(max.len)

all.profils<-NULL
for(i in 1:num.profiles)
 {
  profil<-readBin(ptr, what=double(), size=8, n=profile.lengths[i], signed = TRUE, endian = "little")
  #print(c(length(profil),profile.lengths[i]))
  profil<-c(profil,rep(NA,(max.len-length(profil))))
  all.profils<-rbind(all.profils,profil)
  
  #plot(1:length(profil),profil,typ="l")
  #par(new=T)
  
 }
close(ptr) 

rownames(all.profils)<-NULL

return(all.profils)
 
}


#-------------------------------------------------------------------
#A much faster version of the above. Read in processed binary profiles
#-------------------------------------------------------------------
read.binary.profiles2<-function(file.path)
{
  
  #close(ptr)
  
  #open the file:
  ptr<-file(file.path, "rb")
  
  #Read out the first number = #of profiles:
  seek(ptr, where = 0, rw="r")
  num.profiles<-readBin(ptr, what="integer", size = 4, n = 1, signed = TRUE, endian = "little")
  #print(num.profiles)
  
  #Read out the int vector of profile lengths:
  seek(ptr, where = 4, rw="r")
  profile.lengths <- readBin(ptr, what="integer", size = 4, n = num.profiles, signed = TRUE, endian = "little")
  #print(profile.lengths)
  
  #Read out the int vector of profile lengths:
  #seek(ptr, where =(4 + 4*num.profiles), rw="r")
  
  #Pointer should now be at the beginning of the profile data
  max.len<-max(profile.lengths)
  
  #print(num.profiles)
  #print(profile.lengths)
  #print(max.len)
  
  all.profils<-array(0.0,c(num.profiles,max.len))
  for(i in 1:num.profiles)
  {
    profil<-readBin(ptr, what="double", size=8, n=profile.lengths[i], signed = TRUE, endian = "little")
    #print(c(length(profil),profile.lengths[i]))
    profil<-c(profil,rep(NA,(max.len-length(profil))))
    #all.profils<-rbind(all.profils,profil)
    all.profils[i,] <- profil
    
    #plot(1:length(profil),profil,typ="l")
    #par(new=T)
    
  }
  close(ptr) 
  
  rownames(all.profils)<-NULL
  
  return(all.profils)
  
}