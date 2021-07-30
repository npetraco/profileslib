#--------------------------------------------------
#Measure the "distance"" between two profiles via 
#the ccf
#--------------------------------------------------
ccfmax.distance<-function(prof1,prof2,maxlag,printQ=TRUE)
{
  
croscor<-ccf(na.omit(prof1),na.omit(prof2),lag.max=maxlag,plot=printQ)

maxcorr<-max(croscor[[1]])
maxcorr.idx<-which(croscor[[1]]==max(croscor[[1]]))
lag.at.maxcorr<-croscor$lag[maxcorr.idx]

if(printQ==TRUE)
 {
  print(paste("Maximum of ccf: ", maxcorr, " at lag: ", lag.at.maxcorr)) 
 }

return(c(maxcorr,lag.at.maxcorr))
  
}

#--------------------------------------------------
#Compute ccf max "similarity" scores
#This slow function is taken from the chemometric 
#routenes replacing cor with ccf
#--------------------------------------------------
ccfmax.compare<-function(dat,lbls,actual.lbls,lagmax,score.warning.cutoff,outfil.path)
{
 
#num.grps<-nlevels(lbls)

#IMPORTANT: Assumes group labels are numerical, ordered and start from 1
grp.vec<-as.numeric(levels(lbls))

same.scores<-NULL
diff.scores<-NULL
mscgs<-NULL
for(i in 1:length(grp.vec))
  {
   #Label name for the test group:	
   grp<-grp.vec[i]

   #Get the indices of the test group obs. vecs:
   #Use these for spitting out actual group labels later
   grp.idxs<-which(lbls==grp)
      
   #Labels names of the other groups:
   not.grp<-grp.vec[-i]

   #Get indices of other groups.
   #Unfortunately this is slow: 
   not.grp.idxs<-which(lbls!=grp)
   not.lbl<-as.numeric(lbls[not.grp.idxs])
   #print(length(not.lbl))

   #Compute the ccf max scores between data vecs from the same group
   #Grab vecs from the same group
   tmp.same<-pick.out.groups(dat,lbls,c(grp))
   dat.same<-as.matrix(tmp.same[[1]])

   #Compute ccf max scores between vecs.
   #Note, need more than 1 observation in the test group
   if(length(grp.idxs)>1)
    {
     ###################
     #Surgery needed here
     #
     #tmp.same<-cor(t(dat.same))
     tmp.same<-NULL
     for(j in 1:nrow(dat.same))
      {
       for(k in 1:nrow(dat.same))
        {
         if(j<k) #Don't doulbe count
          {
           croscor<-ccf(dat.same[j,],dat.same[k,],na.action=na.omit,lag.max=lagmax,plot=FALSE)
           tmp.same<-c(tmp.same,max(croscor[[1]]))
           print(paste("SAME Grp:",i,"#",j,"vs. #",k,max(croscor[[1]])))
           #which(croscor[[1]]==max(croscor[[1]]))
           #croscor[[1]][which(croscor[[1]]==max(croscor[[1]]))]
           #croscor$lag[which(croscor[[1]]==max(croscor[[1]]))]
          }	
        }	
      }
     #
     ###################
     #grp.grp.scores<-extract.lower.triangle(tmp.same)
     grp.grp.scores<-tmp.same
     
     #Store the ccf max scores between profiles of the same group.
     same.scores<-c(same.scores,grp.grp.scores)
    }

     
   #Separate out data of other groups   
   tmp.diff<-pick.out.groups(dat,lbls,not.grp)
   dat.diff<-as.matrix(tmp.diff[[1]])
   #print(nrow(dat.diff))
   #print("")

   for(j in 1:nrow(dat.same))
     {
      for(k in 1:nrow(dat.diff))
        {
         #print(paste(grp,"",not.lbl[k] ))
         if(grp < not.lbl[k])
          {
           #print(paste("GOOD!","Compare grp:",grp,"to grp:",not.lbl[k]))
           ########################
           #
           #Surgery needed here
           #
           #grp.notgrp.score<-cor(dat.same[j,],dat.diff[k,])
           croscor<-ccf(dat.same[j,],dat.diff[k,],na.action=na.omit,lag.max=lagmax,plot=FALSE)
           grp.notgrp.score<-max(croscor[[1]])
           print(paste("DIFF Grp:",i,"#",j,"vs. #",k,max(croscor[[1]])))
           #
           ########################
           if(grp.notgrp.score>score.warning.cutoff)
             {
              mscgs<-rbind(mscgs,paste("UH OH!:",actual.lbls[grp.idxs[j]],"vs.",actual.lbls[not.grp.idxs[k]],"CCF-MAX=",grp.notgrp.score))
             }
           diff.scores<-c(diff.scores,grp.notgrp.score) 
           #print(grp.notgrp.score)
          }
#         if(grp > not.lbl[k])
#          {
#           print(paste("BAD!","Grp:",grp,"greater than grp:",not.lbl[k],"DON'T COMPARE AGAIN"))
#          }
#         if(grp == not.lbl[k])
#          {
#           print(paste("********************************REALLY BAD!","Grp:",grp,"is the same as grp:",not.lbl[k],"IT'S NOT SUPPOSED TO BE!!!!!!"))
#          }
          	
        }      
     }

   print(paste("Done with group:",i))

  }

#print(mscgs)
#write.table(mscgs,file=paste(outfil.path,"Warning_Messages_",as.character(corr.warning.cutoff),".txt",sep=""))
#print(paste(outfil.path,"Warning_Messages_",as.character(corr.warning.cutoff),".txt",sep=""))

return(list(same.scores,diff.scores))  
	
}


#--------------------------------------------------
#Compute correlation coefs "similarity" scores but
#give more diagnostic info
#
#--------------------------------------------------
correlation.compare.maxccf<-function(dmat, lbls, actual.lbls, maxlag, sames.warning.cutoff, difs.warning.cutoff, outfil.path)
{

k<-nlevels(lbls) #num grps
n<-nrow(dmat)    #num obs
#Total number of comparisons:
total.num.comp<-choose(n,2)

#obs. indices of unique pairs of comparisons
comp.idxs<-combn(n,2) #note, conbinations are dumped out in column-wise format

info.mat<-matrix(rep("x",total.num.comp*8),nrow=total.num.comp,ncol=8) #initialize a character mat to hold info
#info.mat<-matrix(rep("x",10*6),nrow=10,ncol=6)
print(paste(" Starting correlation score computations:",date()))
tim<-system.time(
for(i in 1:ncol(comp.idxs))
#for(i in 1:10)
{
 idxa<-comp.idxs[1,i]
 idxb<-comp.idxs[2,i]
 grpa<-lbls[idxa]
 grpb<-lbls[idxb]
 actual.grpa<-actual.lbls[idxa]
 actual.grpb<-actual.lbls[idxb]
 sscore.info<-ccf(dmat[idxa,],dmat[idxb,], lag.max=maxlag, type="correlation", plot=FALSE, na.action=na.omit) #compute the similarity scores here
 sscore<-max(sscore.info$acf)
 info.vec<-c(as.character(grpa),as.character(grpb),as.character(actual.grpa), as.character(actual.grpb), grpa==grpb, as.character(idxa),as.character(idxb),sscore)
 info.mat[i,]<-info.vec
 #print(info.vec)
}
)
print(paste("Done with similarity score computations:",date()))
print(tim)
rownames(info.mat)<-NULL
colnames(info.mat)<-c("Grp a","Grp b","Name a","Name b", "Same GrpQ","Obs Idx a","Obs Idx b","Sim score")
#print(info.mat)

num.same.comp<-sum(as.logical(info.mat[,5]))
num.dif.comp<-total.num.comp-num.same.comp

same.idxs<-which(as.logical(info.mat[,5])==TRUE) #pull out same/different group comparison scores
dif.idxs<-which(as.logical(info.mat[,5])==FALSE) 

#Process comparisons between obs of the SAME group
same.mat<-info.mat[same.idxs,] 
#print(same.mat)
same.scores<-as.numeric(same.mat[,8])

problem.idxs.same<-which(same.scores<=sames.warning.cutoff) #pull out the problem scores
problem.same.mat<-same.mat[problem.idxs.same,]
#print(problem.same.mat)

#Process comparisons between obs between DIFFERENT groups
dif.mat<-info.mat[dif.idxs,] 
#print(dif.mat)
dif.scores<-as.numeric(dif.mat[,8])

problem.idxs.dif<-which(dif.scores>=difs.warning.cutoff) #pull out the problem scores
problem.dif.mat<-dif.mat[problem.idxs.dif,]
#print(nrow(problem.dif.mat))

print(paste("Total number of comparisons:           ",total.num.comp))
print(paste("Number of SAME group comparisons:      ",num.same.comp,", ",length(problem.idxs.same)/num.same.comp*100,"% a problem; below threshold: ",sames.warning.cutoff, sep=""))
print(paste("Number of DIFFERENT group comparisons: ",num.dif.comp,", ",length(problem.idxs.dif)/num.dif.comp*100,"% a problem; above threshold: ",difs.warning.cutoff, sep=""))

#Write scores to file for troubleshooting:
write.table(problem.dif.mat,file=paste(outfil.path,"Problem_Observations_Between_Different_Groups",as.character(difs.warning.cutoff),".txt",sep=""))
write.table(problem.same.mat,file=paste(outfil.path,"Problem_Observations_Within_Same_Groups",as.character(sames.warning.cutoff),".txt",sep=""))

write.table(dif.mat,file=paste(outfil.path,"All_Score_Info_Between_Different_Groups",".txt",sep=""))
write.table(same.mat,file=paste(outfil.path,"All_Score_Info_Within_Same_Groups",".txt",sep=""))

colnames(info.mat)<-c("Grp a","Grp b","Name a","Name b", "Same GrpQ","Obs Idx a","Obs Idx b","Sim score")

info.mat<-data.frame(info.mat[,1:4], as.logical(info.mat[,5]), as.numeric(info.mat[,6:8]) )
problem.dif.mat<-data.frame(problem.dif.mat[,1:4], as.logical(problem.dif.mat[,5]), as.numeric(problem.dif.mat[,6:8]) )
problem.same.mat<-data.frame(problem.same.mat[,1:4], as.logical(problem.same.mat[,5]), as.numeric(problem.same.mat[,6:8]) )

return(list(same.scores, dif.scores, info.mat, problem.same.mat, problem.dif.mat))

}


#BELOW BROKEN!!!!!!!!!
#MISSES SOME SCORES!!!!!!!!
#3/20/11, NEEDS FIXING
#------------------------------------------------------
#Compute correlation coefs "similarity" scores faster??
#------------------------------------------------------
correlation.compare2<-function(dat,lbls,actual.lbls,same.corr.warning.cutoff,diff.corr.warning.cutoff)
{

#Compute the correlation scores between OBSERVATION VECTORS.
#This correlation matrix is n x n NOT p x p.
#This function simply parses this matrix:
obs.cor<-cor(t(dat))

#Unique, ordered names of all groups, starting from 1:
grp.vec<-as.numeric(levels(lbls))

#Parse obs.cor matrix: 
same.scores<-NULL
diff.scores<-NULL
for(i in 1:length(grp.vec))
 {
  #Label name for the test group:	
  grp<-grp.vec[i]

  #Get the indices of the test group obs. vecs:
  grp.idxs<-which(lbls==grp)

  #Labels names of the other groups:
  not.grp<-grp.vec[-i]
  
  #Just the names of the not-test groups in the test group column of the 
  #lower triangle of obs.cor matrix:
  not.grp.needed<-which(not.grp<grp)

  #Check that there is something below the grp-grp block and that there 
  #is more than one element in the test group:
  if((length(not.grp.needed)>0) & length(grp.idxs)>1)
   {
   	#--------------Parse a test block (same-same) of scores-----------------
   	
   	#Grab a block of test group scores:
   	same.same.block<-obs.cor[grp.idxs,grp.idxs]
   	
   	#See if any fall below the warning thresh (same-same-problem):
   	ssp.idxs<-which(same.same.block<same.corr.warning.cutoff,arr.ind=TRUE)
   	
   	#Needed to pick out the actual scores:
    ssp.idxs.lin<-which(same.same.block<same.corr.warning.cutoff)
    same.same.block.vec<-as.numeric(same.same.block)
    
    #Put the indices and scores together:
    ssp.idxs<-cbind(ssp.idxs,same.same.block.vec[ssp.idxs.lin])
    
    #Just pick out the unique scores:
   	ssp.idxs<-ssp.idxs[which((ssp.idxs[,1]<ssp.idxs[,2])==TRUE),]
   	ssp.idxs<-as.matrix(ssp.idxs)
   	if(dim(ssp.idxs)[2]==1)
   	 {
   	  ssp.idxs<-t(ssp.idxs)	
   	 }
   	   	
   	#Pull out the unique scores (lower triangle):
    grp.grp.scores<-extract.lower.triangle(same.same.block)
    
    #Store in the test group score vector:
    same.scores<-c(same.scores,grp.grp.scores)
    #------------------------------------------------------------------------

   	#----Parse a test-not/test rectangular block (different-same) of scores--
    
    #Don't go below the last row of obs.cor matrix: 
    if((max(grp.idxs)+1)<nrow(obs.cor))
     {
     	
      #Grab a (rectangle) block of diff-same scores from the lower triangle
      #of obs.cor matrix: 
      grp.notgrp.score<-obs.cor[(max(grp.idxs)+1):nrow(obs.cor),grp.idxs]
      rownames(grp.notgrp.score)<-actual.lbls[(max(grp.idxs)+1):nrow(obs.cor)]
      colnames(grp.notgrp.score)<-actual.lbls[grp.idxs]
      #print(grp.notgrp.score)
      
      #See if any scores are above the warning thresh (diff-same-problem):
      dsp.idxs<-which(grp.notgrp.score>diff.corr.warning.cutoff,arr.ind=TRUE)

      #Needed to pick out the actual scores:
      dsp.idxs.lin<-which(grp.notgrp.score>diff.corr.warning.cutoff)
      diff.same.block.vec<-as.numeric(grp.notgrp.score)
      
      #Put the indices and scores together:
      dsp.idxs<-as.matrix(cbind(dsp.idxs,diff.same.block.vec[dsp.idxs.lin]))
      #print(dim(dsp.idxs))
      if(dim(dsp.idxs)[2]==1)
   	   {
   	    dsp.idxs<-t(dsp.idxs)	
   	   }
      
      #Flatten out the score rectangle into a vector and store:
      diff.scores<-c(diff.scores,as.numeric(grp.notgrp.score))
     } 
    #------------------------------------------------------------------------

    #Printing
   	#Print out only those indices/scores with a problem:
    if((length(ssp.idxs)>0) | (length(dsp.idxs)>0)  )
     {
      print("============================================================")
      tmp<-as.character(actual.lbls[grp.idxs] )
     }
   	if(length(ssp.idxs)>0)
   	 {
      col1s<-sapply(ssp.idxs[,1], function(x){tmp[x] } ) 
      col2s<-sapply(ssp.idxs[,2], function(x){tmp[x] } ) 
      ssp.idxs<-cbind(col1s,col2s,ssp.idxs[,3])
      colnames(ssp.idxs)<-c("Replicate 1#","Replicate 2#","Corr Score" )

   	  print(paste("Warning! SMALL corr WITHIN some observations of this group:"))
      print("-------------------")
      print("Correlation Scores:")
      print("-------------------")
   	  print(ssp.idxs)
   	  
   	 }
    
    if(length(dsp.idxs)>0)
     {
   	  col1d<-rownames(dsp.idxs)
   	  col2d<-sapply(dsp.idxs[,2], function(x){tmp[x] } )
   	  dsp.idxs<-cbind(col2d,col1d,dsp.idxs[,3])
   	  colnames(dsp.idxs)<-c("This Group#","Other Group#","Corr Score")
   	  rownames(dsp.idxs)<-NULL

      print(paste("Warning! LARGE corr BETWEEN some observations of this group and other groups:"))
   	  print("-------------------")
      print("Correlation Scores:")
      print("-------------------")
      print(dsp.idxs)

     }
    if((length(ssp.idxs)>0) | (length(dsp.idxs)>0)  )
     {
      print("============================================================")
      print("")
     }

        
   }
   	
 }

return(list(same.scores,diff.scores)) 

}