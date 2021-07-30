#-------------------------------------------------------------------
#tile.coords = Preprocessed tile mark lists wrt a landmark
#
#There should be -1s at the top, bottom and inbetween the tiles sets 
#from each exemplar.
#This prevents return of proximate tiles from the exemplar
#-------------------------------------------------------------------
proximate.tiles<-function(tile.coords,pattern.names,tolerance,min.num.matches,printQ=0)
{

#Get indices of tiles within each exemplar:
patt.ends<-which(tile.coords[,1]==-1)
num.patts<-length(patt.ends)-1
print(paste("There are ",num.patts," exemplars",sep=""))

for(i in 1:num.patts)
 {
  all.close.tiles<-NULL
  keepQ.mat<-NULL
  for(j in 1:num.patts)
   {
    if(i!=j) #Enter Exemplar patterns i vs j:
     {
      #Grab tile coords for exemplar i:	
      x1<-dat[(patt.ends[i]+1):(patt.ends[i+1]-1),3]
      
      #Grab tile coords for exemplar j:
      x2<-dat[(patt.ends[j]+1):(patt.ends[j+1]-1),3]
      
      #Compute proximities of tiles between exemplars:
      diffs<-outer(x1,x2,FUN="-")
      
      #Translate output distances to query if proximities are within tolerance:
      ldiffs<-as.numeric(abs(diffs)<=tolerance)
      ldiffs<-matrix(ldiffs,nrow=dim(diffs)[1],ncol=dim(diffs)[2])
      close.tiles<-which(ldiffs==1,arr.ind=TRUE)

      #Print out proximate tiles between exemplars:
      colnames(close.tiles)<-c(as.character(pattern.names[i,]),as.character(pattern.names[j,]))
      #print(close.tiles)
      all.close.tiles<-c(all.close.tiles,list(close.tiles))
      
      #This tell which tiles in an exemplar are proximal to at least one other tile from
      #another exemplar:
      keepQ<-colSums(t(ldiffs))>=1      
      keepQ.mat<-rbind(keepQ.mat,keepQ)
     }	
   }#end for j
  
  #This tells which tiles from an exemplar are proximal to enough other tiles from other 
  #exemplars (min.num.matches of them) to bother keeping for comparisons:   
  res<-as.matrix(colSums(keepQ.mat)>=min.num.matches)
  colnames(res)<-as.character(pattern.names[i,])

  res2<-cbind(1:dim(res)[1],res)
  print(paste("======================",as.character(pattern.names[i,]),"======================"))
  print("Keep Tiles:")
  print(res2[which(res2[,2]!=0),1])   	
  
  if(printQ==1) print(paste("======================",as.character(pattern.names[i,]),"======================"))
  
  #Here, go through all close tile lists for an exemplar
  sapply(which(res==TRUE), #"vector" over the tiles of an exemplar that appear enough 
   function(x)
    {
     sapply(1:length(all.close.tiles), #"vector" over list of all.close.tiles for an exemplar
      function(y)
       {
       	#Grab a matrix of close.tiles
       	tmp<-all.close.tiles[[y]]
       	
       	#Knock out those tiles that are not being "vectored" over 
       	tmp2<-cbind(tmp[,1]*(tmp[,1]==x), tmp[,2]*(tmp[,1]==x))
       	
       	#Get rid of all the (0,0)s and retain only the tile indices
        tmp2<-matrix(which(tmp2!=0),ncol=2,byrow=TRUE)
        
        #Put the exemplar names being compaired on the columns
       	colnames(tmp2)<-colnames(tmp)
       	
       	if(printQ==1) print(tmp2)
       	
       }   ) #end sapply2
    }   ) #end sapply1
  
 }#end for i

}


#-------------------------------------------------------------------
#This function just finds tiles that have at least the number of 
#specified proxilal tiles
#
#It reurns these tiles for processing to determine exactly which
#of the other exemplars tiles they are close to.
#-------------------------------------------------------------------
keep.tiles<-function(tile.coords,pattern.names,tolerance,min.num.matches,printQ=0)
{

#Get indices of tiles within each exemplar:
patt.ends<-which(tile.coords[,1]==-1)
num.patts<-length(patt.ends)-1
print(paste("There are ",num.patts," exemplars",sep=""))

keep.mat<-NULL
for(i in 1:num.patts)
 {
  all.close.tiles<-NULL
  keepQ.mat<-NULL
  for(j in 1:num.patts)
   {
    if(i!=j) #Enter Exemplar patterns i vs j:
     {
      #Grab tile coords for exemplar i:	
      x1<-dat[(patt.ends[i]+1):(patt.ends[i+1]-1),3]
      
      #Grab tile coords for exemplar j:
      x2<-dat[(patt.ends[j]+1):(patt.ends[j+1]-1),3]
      
      #Compute proximities of tiles between exemplars:
      diffs<-outer(x1,x2,FUN="-")
      
      #Translate output distances to query if proximities are within tolerance:
      ldiffs<-as.numeric(abs(diffs)<=tolerance)
      ldiffs<-matrix(ldiffs,nrow=dim(diffs)[1],ncol=dim(diffs)[2])
      close.tiles<-which(ldiffs==1,arr.ind=TRUE)

      #Print out proximate tiles between exemplars:
      colnames(close.tiles)<-c(as.character(pattern.names[i,]),as.character(pattern.names[j,]))
      #print(close.tiles)
      all.close.tiles<-c(all.close.tiles,list(close.tiles))
      
      #This tell which tiles in an exemplar are proximal to at least one other tile from
      #another exemplar:
      keepQ<-colSums(t(ldiffs))>=1      
      keepQ.mat<-rbind(keepQ.mat,keepQ)
     }	
   }#end for j
  
  #This tells which tiles from an exemplar are proximal to enough other tiles from other 
  #exemplars (min.num.matches of them) to bother keeping for comparisons:   
  res<-as.matrix(colSums(keepQ.mat)>=min.num.matches)
  colnames(res)<-as.character(pattern.names[i,])

  res2<-cbind(1:dim(res)[1],res)  
  tiles.to.be.kept<-res2[which(res2[,2]!=0),1]
  exemplar.name<-rep(as.character(pattern.names[i,]),length(tiles.to.be.kept))
  exemplar.number<-rep(i,length(tiles.to.be.kept))
  keep.sub.mat<-cbind(exemplar.name,exemplar.number,tiles.to.be.kept)
  keep.mat<-rbind(keep.mat,keep.sub.mat)

  if(printQ==1)
   {
    print(keep.sub.mat)
   }

 }

return(keep.mat)
 
}


#-------------------------------------------------------------------
#
#-------------------------------------------------------------------
close.tiles<-function(tile.coords,tile.info)
{

#FIRST DROP -1s from tile COORDS

#Get indices of tiles within each exemplar:
patt.ends<-which(tile.coords[,1]==-1)
num.each.pattern<-diff(patt.ends,lag=1)-1
patt.ends<-patt.ends[-1]

indices<-NULL
for(i in 1:length(patt.ends))
 {
  indices<-c(indices,c(rep(i,num.each.pattern[i]),-1))
 }
indices<-c(-1,indices)

for(i in 1:nrow(tile.info))
 {
  query.rows<-which(indices==tile.info[i,2])
  #print(query.rows)
  #print(tile.coords[query.rows,2])
  print(query.rows[as.numeric(tile.info[i,3]) ])
 }


#print(cbind(indices,tile.coords[,2]))

#print(length(indices))
#print(length(patt.ends))
#print(length(num.each.pattern))


kept.tile.lbls<-cbind(as.numeric(tile.info[,2]), as.numeric(tile.info[,3]) )

#print(kept.tile.lbls)

}


#-------------------------------------------------------------------
#
#-------------------------------------------------------------------
close.tiles2<-function(tile.coords,pattern.names,tile.info,tol)
{

#Get indices of tiles within each exemplar:
patt.ends<-which(tile.coords[,1]==-1)

#Drop the -1s between the patterns
tile.coords.mod<-tile.coords[-patt.ends,]

#Determine the number of tiles for each pattern
num.each.pattern<-diff(patt.ends,lag=1)-1

#Drop the first element. Not needed
patt.ends<-patt.ends[-1]

#For each pattern, determine the indices of thier tiles
indices<-NULL
for(i in 1:length(patt.ends))
 {
  indices<-c(indices,rep(i,num.each.pattern[i]))
 }
tile.coords.mod<-cbind(indices,tile.coords.mod)
colnames(tile.coords.mod)<-c("pattern#","patch#","tile#","coord")

for(i in 1:nrow(tile.info))
 {
  #Pick out the set of tiles beloning to the pattern of the query idx
  query.rows<-which(indices==as.numeric(tile.info[i,2]))
   
  #Find index and coord of the tile of interest in pattern
  query.idx<-query.rows[as.numeric(tile.info[i,3]) ]
  query.tile.coords.info<-tile.coords.mod[query.idx,]
  
  #Drop the whole set of tiles beloning to the pattern of the query idx
  tmp.tile.coords.mod<-tile.coords.mod[-query.rows,]
  
  #Compute distance between query.coord and remaining tile.coords.mod
  abs.tile.proximities<-abs( (tmp.tile.coords.mod[,4]) - as.numeric(query.tile.coords.info[4]) )
  compare.tile.info<-tmp.tile.coords.mod[which(abs.tile.proximities<=tol),]

  #Print query pattern/tile info out nicely:
  query.pattern.name<-as.character(pattern.names[ as.numeric(query.tile.coords.info[1]), ] )
  
  all.query.pattern.info<-data.frame( c(query.pattern.name, query.tile.coords.info[2:4]) )
  #all.query.pattern.info<-cbind( c(query.pattern.name, query.tile.coords.info[2:4]) )
  
  names(all.query.pattern.info)<-c("pattern name","patch#","tile#","coord")
  print("COMPARE:") 
  print(all.query.pattern.info)
  print("",quote=FALSE)
  
  print("TO THESE:")
  compare.pattern.names<-as.character(pattern.names[as.numeric(compare.tile.info[,1]),])
  print(cbind(compare.pattern.names,compare.tile.info[,2:4]))
  
  print("------------------------------------------")
 }

}