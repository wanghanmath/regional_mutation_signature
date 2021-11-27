mutationname<-read.csv("/data/nameof96mutation.csv",header=FALSE,stringsAsFactors=FALSE)
lambda<-5
mutationname<-as.matrix(mutationname)
#calculate log likelihood in a bin#
probmultinom<-function(count,prob){
logprob<-0
index<-which(prob!=0)
for (i in 1:length(index)){
logprob<-logprob+count[index[i]]*log(prob[index[i]])
}
return(list("likelihood"=logprob))
}
#if the location-th bin and (location+1)-th bin are merged, calculate the log likelihood of the newly formed bin#
mergeforthdir<-function(location){
  countmatrixmerge<-matrix(ncol=96,nrow=1)
  countmatrix1<-matrix(ncol=96,nrow=1)
  countmatrix2<-matrix(ncol=96,nrow=1)
  for (k in 1:96){
    countmatrix1[k]<-length(which(sample[(storagelocation[location,1]:storagelocation[location,2]),6]==mutationname[k]))
    countmatrix2[k]<-length(which(sample[(storagelocation[(location+1),1]:storagelocation[(location+1),2]),6]==mutationname[k]))
    countmatrixmerge[k]<-length(which(sample[(storagelocation[location,1]:storagelocation[(location+1),2]),6]==mutationname[k]))
  }
  mergelikelihood<-probmultinom(count=c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(count=c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
  return(list("likelihood"=mergelikelihood))
}
#read the mutation file#
  binsize<-30 #we choose bin size to be 30#
  sampleoverall<-read.csv("/data/sample0.txt",sep=" ",header=F,stringsAsFactors=F)
  chromosome<-as.character(unique(sampleoverall[,2]))
  storageoverall<-matrix(ncol=100,nrow=1)
  for (chr in 1:length(chromosome)){
  sample<-sampleoverall[which(sampleoverall[,2]==chromosome[chr]),]
  bincount<-floor(dim(sample)[1]/binsize)
  #bincount is the number of bins when cutting the chromosome using 30-length bins#
  if (bincount==0){
  countmatrix<-matrix(ncol=96,nrow=1)
    for (k in 1:96){
      countmatrix[k]<-length(which(sample[c(1:dim(sample)[1]),6]==mutationname[k]))
    }
  storagelocation<-matrix(c(1,dim(sample)[1]),ncol=2)
  storagechr<-matrix(c(chromosome[chr],(storagelocation[1,2]-storagelocation[1,1]+1)),ncol=2)
  storagelocation[1,1]<-sample[storagelocation[1,1],3]
  storagelocation[1,2]<-sample[storagelocation[1,2],3]
  storagelocationnew<-cbind(storagechr,storagelocation)
  storagelocationnew<-cbind(storagelocationnew,countmatrix)
  }
if (bincount>=2){
  storagelocation<-matrix(ncol=2,nrow=bincount)
  for (bin in 1:(bincount-1)){
    storagelocation[bin,1]<-(bin-1)*binsize+1
    storagelocation[bin,2]<-bin*binsize
  }
  storagelocation[bincount,1]<-(bincount-1)*binsize+1
  storagelocation[bincount,2]<-dim(sample)[1]
  likelihoodmatrix<-matrix(ncol=bincount,nrow=1)
  mergeinformation<-matrix(ncol=bincount,nrow=3)
  }
  if (bincount==1){
    countmatrix<-matrix(ncol=96,nrow=1)
    for (k in 1:96){
      countmatrix[k]<-length(which(sample[c(1:dim(sample)[1]),6]==mutationname[k]))
    }
	storagelocation<-matrix(c(1,dim(sample)[1]),ncol=2)
	storagechr<-matrix(c(chromosome[chr],(storagelocation[1,2]-storagelocation[1,1]+1)),ncol=2)
	storagelocation[1,1]<-sample[storagelocation[1,1],3]
	storagelocation[1,2]<-sample[storagelocation[1,2],3]
  storagelocationnew<-cbind(storagechr,storagelocation)
	storagelocationnew<-cbind(storagelocationnew,countmatrix)
  }
  if (bincount==2){
    for (j in 1:(bincount-1)){
      countmatrix<-matrix(ncol=96,nrow=1)
      for (k in 1:96){
        countmatrix[k]<-length(which(sample[c(((j-1)*binsize+1):(j*binsize)),6]==mutationname[k]))
      }
      likelihoodmatrix[j]<-probmultinom(count=c(countmatrix),prob=c(countmatrix)/binsize)$likelihood
    }
    for (k in 1:96){
      countmatrix[k]<-length(which(sample[c(((bincount-1)*binsize+1):dim(sample)[1]),6]==mutationname[k]))
    }
    likelihoodmatrix[bincount]<-probmultinom(count=c(countmatrix),prob=c(countmatrix)/sum(countmatrix))$likelihood
    mergeinformation[1,]<-c(1:bincount)
    mergeinformation[2,]<-rep(1,bincount)
	  countmatrixmerge<-matrix(ncol=96,nrow=1)
    countmatrix1<-matrix(ncol=96,nrow=1)
    countmatrix2<-matrix(ncol=96,nrow=1)
    for (k in 1:96){
      countmatrix1[k]<-length(which(sample[c(((bincount-2)*binsize+1):((bincount-1)*binsize)),6]==mutationname[k]))
      countmatrix2[k]<-length(which(sample[c(((bincount-1)*binsize+1):dim(sample)[1]),6]==mutationname[k]))
      countmatrixmerge[k]<-length(which(sample[c(((bincount-2)*binsize+1):dim(sample)[1]),6]==mutationname[k]))
    }
    mergelikelihood<-probmultinom(count=c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(count=c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
    mergeinformation[3,(bincount-1)]<-likelihoodmatrix[bincount-1]+likelihoodmatrix[bincount]-mergelikelihood
    if (mergeinformation[3,1]<(lambda*log(dim(sample)[1]))){
      storagelocation<-t(as.matrix(c(1,dim(sample)[1])))
      countmatrix<-matrix(ncol=96,nrow=1)
      for (k in 1:96){
        countmatrix[k]<-length(which(sample[c(1:dim(sample)[1]),6]==mutationname[k]))
      }
	storagechr<-matrix(c(chromosome[chr],(storagelocation[1,2]-storagelocation[1,1]+1)),ncol=2)
	storagelocation[1,1]<-sample[storagelocation[1,1],3]
	storagelocation[1,2]<-sample[storagelocation[1,2],3]
  storagelocationnew<-cbind(storagechr,storagelocation)
	storagelocationnew<-cbind(storagelocationnew,countmatrix)
    }
 if (mergeinformation[3,1]>=lambda*log(dim(sample)[1])){
      countmatrix<-rbind(countmatrix1,countmatrix2)
	  storagechr<-matrix(ncol=2,nrow=2)
	  storagechr[,1]<-rep(chromosome[chr],2)
	  storagechr[1,2]<-storagelocation[1,2]-storagelocation[1,1]+1
	  storagechr[2,2]<-storagelocation[2,2]-storagelocation[2,1]+1
	  storagelocation[1,1]<-sample[storagelocation[1,1],3]
	  storagelocation[1,2]<-sample[storagelocation[1,2],3]
	  storagelocation[2,1]<-sample[storagelocation[2,1],3]
	  storagelocation[2,2]<-sample[storagelocation[2,2],3]
    storagelocationnew<-cbind(storagechr,storagelocation)
	  storagelocationnew<-cbind(storagelocationnew,countmatrix)
    }
  }
  if (bincount>=3){
  for (j in 1:(bincount-1)){
    countmatrix<-matrix(ncol=96,nrow=1)
    for (k in 1:96){
      countmatrix[k]<-length(which(sample[c(((j-1)*binsize+1):(j*binsize)),6]==mutationname[k]))
    }
    likelihoodmatrix[j]<-probmultinom(count=c(countmatrix),prob=c(countmatrix)/binsize)$likelihood
  }
  for (k in 1:96){
    countmatrix[k]<-length(which(sample[c(((bincount-1)*binsize+1):dim(sample)[1]),6]==mutationname[k]))
  }
  likelihoodmatrix[bincount]<-probmultinom(count=c(countmatrix),prob=c(countmatrix)/sum(countmatrix))$likelihood
  mergeinformation[1,]<-c(1:bincount)
  mergeinformation[2,]<-rep(1,bincount)
  #the 3rd row ,i-th column of mergeinformation contains similarity of the i-th bin and (i+1)-th bin, if this value is small, smaller than lambda*log(dim(sample)[1]), merge the i-th and (i+1)-th bin#
  for (j in 1:(bincount-2)){
    countmatrixmerge<-matrix(ncol=96,nrow=1)
    countmatrix1<-matrix(ncol=96,nrow=1)
    countmatrix2<-matrix(ncol=96,nrow=1)
    for (k in 1:96){
      countmatrix1[k]<-length(which(sample[c(((j-1)*binsize+1):(j*binsize)),6]==mutationname[k]))
      countmatrix2[k]<-length(which(sample[c((j*binsize+1):(j*binsize+binsize)),6]==mutationname[k]))
      countmatrixmerge[k]<-length(which(sample[c(((j-1)*binsize+1):(j*binsize+binsize)),6]==mutationname[k]))
    }
    mergelikelihood<-probmultinom(count=c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(count=c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
    mergeinformation[3,j]<-likelihoodmatrix[j]+likelihoodmatrix[j+1]-mergelikelihood
  }
  for (k in 1:96){
    countmatrix1[k]<-length(which(sample[c(((bincount-2)*binsize+1):((bincount-1)*binsize)),6]==mutationname[k]))
    countmatrix2[k]<-length(which(sample[c(((bincount-1)*binsize+1):dim(sample)[1]),6]==mutationname[k]))
    countmatrixmerge[k]<-length(which(sample[c(((bincount-2)*binsize+1):dim(sample)[1]),6]==mutationname[k]))
  }
  mergelikelihood<-probmultinom(count=c(countmatrix1),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood+probmultinom(count=c(countmatrix2),prob=c(countmatrixmerge)/sum(countmatrixmerge))$likelihood
  mergeinformation[3,(bincount-1)]<-likelihoodmatrix[bincount-1]+likelihoodmatrix[bincount]-mergelikelihood
  mergeinformationnew<-mergeinformation
#the i-th element in the likelihoodinformation is the log likelihood of the i-th bin#
  likelihoodinformation<-likelihoodmatrix
  mergelocation<-c()
  iternumber<-dim(storagelocation)[1]-2
  while(length(which(mergeinformation[3,]<(lambda*log(dim(sample)[1]))))>0&length(mergelocation)<iternumber){
    bbb<-0
    aa<-which(mergeinformation[3,c(1:(dim(mergeinformation)[2]-1))]<(lambda*log(dim(sample)[1])))
    location<-aa[which.min(mergeinformation[3,aa])]
    mergelocation<-c(mergelocation,location)
    if (location!=(dim(as.matrix(mergeinformation))[2]-1)&location!=1){
      bbb<-1
      mergeinformation[2,location]<-mergeinformation[2,location]+mergeinformation[2,location+1]
      storagelocation[location,2]<-storagelocation[(location+1),2]
      countmatrix<-matrix(ncol=96,nrow=1)
      for (k in 1:96){
        countmatrix[k]<-length(which(sample[c(storagelocation[location,1]:storagelocation[location,2]),6]==mutationname[k]))
      }
      likelihoodinformation[location]<-probmultinom(count=c(countmatrix),prob=c(countmatrix)/sum(countmatrix))$likelihood
      mergeinformation<-mergeinformation[,-(location+1)]
      likelihoodinformation<-likelihoodinformation[-(location+1)]
      storagelocation<-storagelocation[-(location+1),]
      mergeinformation[3,location]<-likelihoodinformation[location]+likelihoodinformation[location+1]-mergeforthdir(location)$likelihood
      mergeinformation[3,(location-1)]<-likelihoodinformation[location-1]+likelihoodinformation[location]-mergeforthdir(location-1)$likelihood
    }
    if (location==1){	
      bbb<-1
      mergeinformation[2,location]<-mergeinformation[2,location]+mergeinformation[2,location+1]	
      storagelocation[location,2]<-storagelocation[(location+1),2]
      countmatrix<-matrix(ncol=96,nrow=1)
      for (k in 1:96){
        countmatrix[k]<-length(which(sample[c(storagelocation[location,1]:storagelocation[location,2]),6]==mutationname[k]))
      }
      likelihoodinformation[location]<-probmultinom(count=c(countmatrix),prob=c(countmatrix)/sum(countmatrix))$likelihood
      mergeinformation<-mergeinformation[,-(location+1)]
      likelihoodinformation<-likelihoodinformation[-(location+1)]
      storagelocation<-storagelocation[-(location+1),]
      mergeinformation[3,location]<-likelihoodinformation[location]+likelihoodinformation[location+1]-mergeforthdir(location)$likelihood
    }
 if (location==(dim(as.matrix(mergeinformation))[2]-1)&bbb==0){
      mergeinformation[2,location]<-mergeinformation[2,location]+mergeinformation[2,(location+1)]	
      storagelocation[location,2]<-storagelocation[(location+1),2]
      countmatrix<-matrix(ncol=96,nrow=1)
      for (k in 1:96){
        countmatrix[k]<-length(which(sample[c(storagelocation[location,1]:storagelocation[location,2]),6]==mutationname[k]))
      }
      likelihoodinformation[location]<-probmultinom(count=c(countmatrix),prob=c(countmatrix)/sum(countmatrix))$likelihood
      mergeinformation<-mergeinformation[,-(location+1)]
      likelihoodinformation<-likelihoodinformation[-(location+1)]
      storagelocation<-storagelocation[-(location+1),]
      mergeinformation[3,location]<-NA
      mergeinformation[3,(location-1)]<-likelihoodinformation[location-1]+likelihoodinformation[location]-mergeforthdir(location-1)$likelihood
    }
  }
  countmatrixnew<-matrix(ncol=96,nrow=dim(storagelocation)[1])
  for (location in 1:dim(storagelocation)[1]){
    for (k in 1:96){
      countmatrixnew[location,k]<-length(which(sample[(storagelocation[location,1]:storagelocation[location,2]),6]==mutationname[k]))
    }
  }
  storagechr<-matrix(ncol=2,nrow=dim(storagelocation)[1])
  storagechr[,1]<-c(rep(chromosome[chr],dim(storagelocation)[1]))
  for (j in 1:dim(storagechr)[1]){
  storagechr[j,2]<-as.numeric(storagelocation[j,2])-as.numeric(storagelocation[j,1])+1
  }
  for (j in 1:dim(storagechr)[1]){
  storagelocation[j,1]<-sample[storagelocation[j,1],3]
  storagelocation[j,2]<-sample[storagelocation[j,2],3]
  }
  storagelocationnew<-cbind(storagechr,storagelocation)
  storagelocationnew<-cbind(storagelocationnew,countmatrixnew)
  }
  storageoverall<-rbind(storageoverall,storagelocationnew)
  }
  storageoverall<-storageoverall[c(2:dim(storageoverall)[1]),]
  colnames(storageoverall)<-c("chr","num_of_mutation","location_start","location_end",c(mutationname)) 
  write.csv(storageoverall,"/data/sample0_segmentation.csv",row.names=F)
