#calculate the proportion of each mutation signature in each genomic segment, COSMIC V3#
source("/data/deconstructSigs/datasets.R")
source("/data/deconstructSigs/deconstructSigs.R")
source("/data/deconstructSigs/globalVariables.R")
source("/data/deconstructSigs/golden_section_search.R")
source("/data/deconstructSigs/internal.scripts.R")
source("/data/deconstructSigs/mut.to.sigs.input.R") 
source("/data/deconstructSigs/normalize.data.R")
source("/data/deconstructSigs/plotting.R")
source("/data/deconstructSigs/vcf.to.sigs.input.R") 
source("/data/deconstructSigs/whichSignatures.R")
library("BSgenome.Hsapiens.UCSC.hg19")
load("/data/signatures.genome.cosmic.v3.may2019.rda")
load("/data/sample.mut.ref.rda")
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,sample.id = "Sample",chr = "chr",pos = "pos",ref = "ref",alt = "alt")
mergesample<-read.csv("/data/sample0_segmentation.csv",header=TRUE,sep=",")
mergesamplenew<-as.matrix(mergesample[,c(5:100)])
rownames(mergesamplenew)<-c(1:dim(mergesamplenew)[1])
mutationmatrix<-mergesamplenew
colnames(mutationmatrix)<-colnames(sigs.input)
rownames(mutationmatrix)<-c(1:dim(mutationmatrix)[1])
mutationmatrix<-as.data.frame(mutationmatrix)
proportionmatrix<-matrix(ncol=68,nrow=dim(mergesamplenew)[1])
for (j in 1:dim(mutationmatrix)[1]){
signature<-whichSignatures(tumor.ref =mutationmatrix,signatures.ref = signatures.genome.cosmic.v3.may2019,sample.id = j,contexts.needed = TRUE,tri.counts.method = 'default')$weights
supp<-1-sum(signature)
signature<-t(as.matrix(c(signature,supp)))
name_proportion<-colnames(as.matrix(signature))
name_proportion[68]<-c("error term")
signaturenew<-matrix(ncol=68,nrow=1)
for (m in 1:68){
signaturenew[m]<-signature[m][[1]]
}
proportionmatrix[j,]<-c(as.matrix(signaturenew))
}
colnames(proportionmatrix)<-name_proportion
write.csv(proportionmatrix,"/data/sample0_prop_segment.csv",row.names=F)
