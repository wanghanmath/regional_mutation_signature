# regional_mutation_signature
To account for the heterogeneity of mutation signatures along the genome, we implement genome segmentation and clustering based on BIC-seq algorithm, and then assign a mutation signature to each mutation based on the maximum likelihood method.
segmentation.R and cluster.R are used for segmentation and clustering of genome respectively.
proportion_seg.R and proportion_clu.R are used to calculate the proportion of each mutation signature in each segmentation or cluster respectively.
likelihood.R is used to assign a mutation signature to each mutation along the genome.
