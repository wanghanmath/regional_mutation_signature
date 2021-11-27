#in each segment/cluster, calculate posterior probability that a certain mutation is attributed to each mutation signature#
posterior_seglevel<-function(seg_post){
posterior<-array(dim=c(96,67,dim(seg_post)[1]))
for (j in 1:dim(seg_post)[1]){
sample_sigweight1<-seg_post[j,c(1:67)]
sample_sigweight1<-as.numeric(as.matrix(sample_sigweight1))
sample_sigweight1[c(1:67)]<-sample_sigweight1[c(1:67)]/sum(sample_sigweight1[c(1:67)])
if (length(which(sample_sigweight1==1))>0){
	loc<-which(sample_sigweight1==1)
	loc_diff<-setdiff(c(1:67),loc)
	posterior[,loc,j]<-1
    posterior[,c(loc_diff),j]<-0
}
if (length(which(sample_sigweight1==1))==0){
mutationvector<-rep(0,96)
for (k in 1:96){
for (m in 1:67){
posterior[k,m,j]<-sample_sigweight1[m]*sigprob[k,m]
}
mutationvector[k]<-sum(posterior[k,,j])
posterior[k,,j]<-posterior[k,,j]/mutationvector[k]
}
}
}
row.names<-c(mutationname)
column.names<-c(colnames(sigprob))
matrix.names<-c(1:dim(posterior)[3])
dimnames(posterior)<-list(row.names,column.names,matrix.names)
return(posterior)
}

#example#
sigprob<-read.table("/data/hanwang/signature/COSMIC_v3_SBS_GRCh37.txt",header=T)
mutationname<-as.character(sigprob[,1])
sigprob<-sigprob[,c(2:68)]
seg_post<-read.csv("/data/sample0_prop_segment.csv",header=TRUE,sep=",")
mergesample<-read.csv("/data/sample0_segmentation.csv",header=T,sep=",")
mutsample<-read.csv("/data/sample0.txt",sep=" ",header=FALSE,stringsAsFactors=FALSE)
mutsample<-as.matrix(mutsample)
clustermatrix<-read.csv("/data/sample0_cluster.csv",header=T)
clusterprop<-read.csv("/data/sample0_prop_cluster.csv",header=T,sep=",")
mergenummatrix<-matrix(ncol=2,nrow=dim(mergesample)[1])
mergenummatrix[1,1]<-1
for (i in 1:dim(mergenummatrix)[1]){
	mergenummatrix[i,2]<-sum(mergesample[c(1:i),2])
}
for (i in 2:dim(mergenummatrix)[1]){
	mergenummatrix[i,1]<-mergenummatrix[i-1,2]+1
}
posterior<-posterior_seglevel(seg_post)
posterior_clus<-posterior_seglevel(clusterprop)
likelihoodmatrix<-matrix(ncol=7,nrow=dim(mutsample)[1])
colnames(likelihoodmatrix)<-c("chr","location","mutation","segmentID","clusterID","sig_seglevel","sig_clustlevel")
for (i in 1:dim(mergenummatrix)[1]){
	posterior_seg<-posterior[,,i]
	clusterid<-which(clustermatrix==i,arr.ind=TRUE)[1]
	for (j in mergenummatrix[i,1]:mergenummatrix[i,2]){
		likelihoodmatrix[j,1:2]<-c(mutsample[j,2:3])
		likelihoodmatrix[j,3]<-mutsample[j,6]
		likelihoodmatrix[j,4]<-i
		likelihoodmatrix[j,5]<-clusterid
		mut_loc<-which(mutationname==mutsample[j,6])
		mut_maxloc<-which.max(posterior_seg[mut_loc,])
		sig_num<-colnames(sigprob)[mut_maxloc]
        likelihoodmatrix[j,6]<-sig_num
	}
}
for (i in 1:dim(likelihoodmatrix)[1]){
	clusid<-likelihoodmatrix[i,5]
	posterior_clusid<-posterior_clus[,,clusid]
	mut_loc<-which(mutationname==mutsample[i,6])
	mut_maxloc<-which.max(posterior_clusid[mut_loc,])
	sig_num<-colnames(sigprob)[mut_maxloc]
    likelihoodmatrix[i,7]<-sig_num
}
write.csv(likelihoodmatrix,"/data/sample0_likelihood.csv",row.names=F)



