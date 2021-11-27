probmultinom<-function(count,prob){
logprob<-0
index<-which(prob!=0)
for (i in 1:length(index)){
logprob<-logprob+count[index[i]]*log(prob[index[i]])
}
return(list("likelihood"=logprob))
}
#merging non-adjacent bins#
mutationname<-read.csv("/data/nameof96mutation.csv",header=FALSE,stringsAsFactors=FALSE)
mutationname<-as.matrix(mutationname)
merge_before<-read.csv("/data/sample0_segmentation.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
storagelocation<-matrix(nrow=dim(merge_before)[1],ncol=2)
storagelocation[1,]<-c(1,merge_before[1,2])
for (i in 2:dim(storagelocation)[1]){
  storagelocation[i,1]<-sum(merge_before[1:(i-1),2])+1
  storagelocation[i,2]<-sum(merge_before[1:i,2])
}
sample<-read.csv("/data/sample0.txt",sep=" ",header=FALSE,stringsAsFactors=FALSE)
likelihoodmatrix<-matrix(ncol=dim(storagelocation)[1],nrow=1)
for (i in 1:dim(storagelocation)[1]){
  countmatrix<-rep(0,96)
  for (m in 1:96){
    countmatrix[m]<-length(which(sample[c(storagelocation[i,1]:storagelocation[i,2]),6]==mutationname[m]))
  }
  likelihoodmatrix[i]<-probmultinom(c(countmatrix),prob=c(countmatrix)/sum(countmatrix))$likelihood
}
if (dim(storagelocation)[1]==1){
clustermatrix<-1
}
if (dim(storagelocation)[1]>=2){
clustering<-list(NULL)
#clustering is a list with its j-th element contains the index of segments in the j-th cluster#
length(clustering)<-dim(storagelocation)[1]
for (j in 1:dim(storagelocation)[1]){
  clustering[[j]]<-c(j)
}
#mergematrix is an upper triangular matrix with (j,k) element denotes the similarity of the j-th segment and the k-th segment#
mergematrix<-matrix(nrow=dim(storagelocation)[1],ncol=dim(storagelocation)[1])
countmatrix1<-countmatrix2<-countmatrixmerge<-rep(0,96)
#now lambda=5, larger lambda induces smaller number of clusters#
lambda<-5
thres<-5*log(dim(sample)[1])
for (j in 1:(dim(storagelocation)[1]-1)){
  for (k in (j+1):dim(storagelocation)[1]){
    for (m in 1:96){
      countmatrix1[m]<-length(which(sample[c(storagelocation[j,1]:storagelocation[j,2]),6]==mutationname[m]))
      countmatrix2[m]<-length(which(sample[c(storagelocation[k,1]:storagelocation[k,2]),6]==mutationname[m]))
      countmatrixmerge[m]<-length(which(sample[c(union(c(storagelocation[j,1]:storagelocation[j,2]),c(storagelocation[k,1]:storagelocation[k,2]))),6]==mutationname[m]))
    }
    mergelikelihood<-probmultinom(c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
    mergematrix[j,k]<-likelihoodmatrix[j]+likelihoodmatrix[k]-mergelikelihood
  }
}
count<-0
while (length(which(mergematrix<=thres))){
  #minlocation[1] and minlocation[2] are the most similar segments, if the similarity is smaller than threshold, merge these 2 segments and update the similarity matrix mergematrix#
  count<-count+1
  minlocation<-c(which(mergematrix==mergematrix[which.min(mergematrix)],arr.ind=T))
  if (length(minlocation)>2){
    halflength<-length(minlocation)/2
    minlocation<-c(minlocation[1],minlocation[1+halflength])
  }
  clustering[[minlocation[1]]]<-c(clustering[[minlocation[1]]],clustering[[minlocation[2]]])
  clustering[[minlocation[2]]]<-NULL
  nowlength<-dim(mergematrix)[1]
  if (nowlength==2){
    mergematrix<-NA
  }
  if (nowlength>2){
  mergematrix<-mergematrix[,-(minlocation[2])]
  mergematrix<-mergematrix[-(minlocation[2]),]
  likelihoodmatrix<-likelihoodmatrix[-(minlocation[2])]   
  if (minlocation[1]!=1&minlocation[2]!=nowlength){
    locationunion1<-c()
    for (k in 1:length(clustering[[minlocation[1]]])){
      locationunion1<-c(locationunion1,c(storagelocation[clustering[[minlocation[1]]][k],1]:storagelocation[clustering[[minlocation[1]]][k],2]))
    }
    for (m in 1:96){
      countmatrix1[m]<-length(which(sample[c(locationunion1),6]==mutationname[m]))
    }
    likelihoodmatrix[minlocation[1]]<-probmultinom(c(countmatrix1),prob=c(countmatrix1)/sum(countmatrix1))$likelihood
   for (j in 1:(minlocation[1]-1)){
    locationunion2<-c()
    for (k in 1:length(clustering[[j]])){
      locationunion2<-c(locationunion2,c(storagelocation[clustering[[j]][k],1]:storagelocation[clustering[[j]][k],2]))
    }
    for (m in 1:96){
      countmatrix2[m]<-length(which(sample[c(locationunion2),6]==mutationname[m]))
      countmatrixmerge[m]<-length(which(sample[c(union(locationunion1,locationunion2)),6]==mutationname[m]))
    }
    mergelikelihood<-probmultinom(c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
    mergematrix[j,minlocation[1]]<-likelihoodmatrix[j]+likelihoodmatrix[minlocation[1]]-mergelikelihood
   }
    for (j in (minlocation[1]+1):dim(mergematrix)[1]){
      locationunion2<-c()
      for (k in 1:length(clustering[[j]])){
        locationunion2<-c(locationunion2,c(storagelocation[clustering[[j]][k],1]:storagelocation[clustering[[j]][k],2]))
      }
      for (m in 1:96){
        countmatrix2[m]<-length(which(sample[c(locationunion2),6]==mutationname[m]))
        countmatrixmerge[m]<-length(which(sample[c(union(locationunion1,locationunion2)),6]==mutationname[m]))
      }
      mergelikelihood<-probmultinom(c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
      mergematrix[minlocation[1],j]<-likelihoodmatrix[j]+likelihoodmatrix[minlocation[1]]-mergelikelihood  
    }
  }
  if (minlocation[1]==1&minlocation[2]!=nowlength){
    locationunion1<-c()
    for (k in 1:length(clustering[[minlocation[1]]])){
      locationunion1<-c(locationunion1,c(storagelocation[clustering[[minlocation[1]]][k],1]:storagelocation[clustering[[minlocation[1]]][k],2]))
    }
    for (m in 1:96){
      countmatrix1[m]<-length(which(sample[c(locationunion1),6]==mutationname[m]))
    }
    likelihoodmatrix[minlocation[1]]<-probmultinom(c(countmatrix1),prob=c(countmatrix1)/sum(countmatrix1))$likelihood
    for (j in (minlocation[1]+1):dim(mergematrix)[1]){
      locationunion2<-c()
      for (k in 1:length(clustering[[j]])){
        locationunion2<-c(locationunion2,c(storagelocation[clustering[[j]][k],1]:storagelocation[clustering[[j]][k],2]))
      }
      for (m in 1:96){
        countmatrix2[m]<-length(which(sample[c(locationunion2),6]==mutationname[m]))
        countmatrixmerge[m]<-length(which(sample[c(union(locationunion1,locationunion2)),6]==mutationname[m]))
      }
      mergelikelihood<-probmultinom(c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
      mergematrix[minlocation[1],j]<-likelihoodmatrix[j]+likelihoodmatrix[minlocation[1]]-mergelikelihood  
    }
  }
  }
}
#numcluster is the number of clusters after merging/clustering segments#
numcluster<-length(clustering)
numincluster<-rep(0,numcluster)
for (m in 1:numcluster){
  numincluster[m]<-length(clustering[[m]])
}
#clustermatrix is a matrix with the m-th row containing segment index in the m-th cluster#
#due to different numbers of segments in each cluster, in order to write in matrix form, we add 0 in rows whose corresponding clusters are not containing many segments#
#for example, if the m-th row (the m-th cluster) contains the largest number of segments: max(numincluster), the i-th row contains only 2 segments: i_1 and i_2, then the i-th row of "clustermatrix" contains i_1, i_2 and 0 with count max(numincluster)-2 #
clustermatrix<-matrix(nrow=numcluster,ncol=max(numincluster))
for (m in 1:numcluster){
  if (numincluster[m]!=max(numincluster)){
    clustermatrix[m,c(1:numincluster[m])]<-c(clustering[[m]])
    clustermatrix[m,c((numincluster[m]+1):max(numincluster))]<-rep(0,(max(numincluster)-numincluster[m])) 
  }
  if (numincluster[m]==max(numincluster)){
    clustermatrix[m,]<-c(clustering[[m]])
  }
}
}
write.csv(clustermatrix,"/data/sample0_cluster.csv",row.names=F)


