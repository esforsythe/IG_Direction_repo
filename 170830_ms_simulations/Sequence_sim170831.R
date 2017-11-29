#This is a script for simulating sequence evolution during speciation and introgression
#I'll aslo use this script to test my test for directionality of IG.

#Set the working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/IG_Direction/170830_ms_simulations")

install.packages("phyclust")
library(phyclust)
library(parallel)
library(plyr)
install.packages("seqinr")
library(seqinr)

#Make function for replicating sequence simulation
full_sim<-function(){

#simulate an alignment in which no IG occurs
#set.seed(1234)
#Set parameters for no IG
ret.ms <- ms(nsam = 4, nreps = 1, opts = "-T -t 50 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -r 5 5000")

tree.anc <- read.tree(text = ret.ms[3])
#plot.phylo(tree.anc)
#set.seed(123)
#Use seqgen to simulate sequences
seqs<-seqgen(opts = "-mHKY -l5000 -s 0.01", newick.tree = ret.ms[3])

#Extract each sequence
seq1<-strsplit(c2s(substring(text = (seqs[grep("s1", seqs)]), 11, 5010)), split = "")
seq2<-strsplit(c2s(substring(text = (seqs[grep("s2", seqs)]), 11, 5010)), split = "")
seq3<-strsplit(c2s(substring(text = (seqs[grep("s3", seqs)]), 11, 5010)), split = "")
seq4<-strsplit(c2s(substring(text = (seqs[grep("s4", seqs)]), 11, 5010)), split = "")

#Bind sequences into a matrix and convert to DNAbin format
dnabin<-as.DNAbin(rbind(seq1, seq2, seq3, seq4)) 

#calculate distance matrix
dist_mat<-as.matrix(dist.dna(dnabin))

#Extract distances from matrix
S1sp<-dist_mat[2, 3]
S2sp<-dist_mat[1, 2]
S3sp<-dist_mat[1, 3]

### Now simulate IG from P3 to P2
#set.seed(1234)
#Here I'll try to make the f = 0.6
ret.msIG <- ms(nsam = 4, nreps = 1, opts = "-T -t 50 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -es 0.1 2 0.4 -ej 0.2 5 3 -r 5 5000")

tree.ancIG <- read.tree(text = ret.msIG[3])
#plot.phylo(tree.ancIG)
#set.seed(123)

#Use seqgen to simulate sequences
seqsIG<-seqgen(opts = "-mHKY -l5000 -s 0.01", newick.tree = ret.msIG[3])

#Extract each sequence
seqIG1<-strsplit(c2s(substring(text = (seqsIG[grep("s1", seqsIG)]), 11, 5010)), split = "")
seqIG2<-strsplit(c2s(substring(text = (seqsIG[grep("s2", seqsIG)]), 11, 5010)), split = "")
seqIG3<-strsplit(c2s(substring(text = (seqsIG[grep("s3", seqsIG)]), 11, 5010)), split = "")
seqIG4<-strsplit(c2s(substring(text = (seqsIG[grep("s4", seqsIG)]), 11, 5010)), split = "")

#Bind sequences into a matrix and convert to DNAbin format
dnabinIG<-as.DNAbin(rbind(seqIG1, seqIG2, seqIG3, seqIG4)) 

#calculate distance matrix
dist_matIG<-as.matrix(dist.dna(dnabinIG))

#Extract distances from matrix
S1ig<-dist_matIG[2, 3]
S2ig<-dist_matIG[1, 2]
S3ig<-dist_matIG[1, 3]

#Now plot all of the distances (3 to 2 IG)
#plot(c(S1sp, S1ig, S2sp, S2ig, S3sp, S3ig))

### Now simulate IG from P2 to P3
#set.seed(1234)

#Here I'll try to make the f = 0.6
ret.msIG2_3 <- ms(nsam = 4, nreps = 1, opts = "-T -t 50 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -es 0.1 3 0.4 -ej 0.2 5 2 -r 5 5000")

tree.ancIG2_3 <- read.tree(text = ret.msIG2_3[3])
#plot.phylo(tree.ancIG2_3)
#set.seed(123)

#Use seqgen to simulate sequences
seqsIG2_3<-seqgen(opts = "-mHKY -l5000 -s 0.01", newick.tree = ret.msIG2_3[3])

#Extract each sequence
seqIG2_3_1<-strsplit(c2s(substring(text = (seqsIG2_3[grep("s1", seqsIG2_3)]), 11, 5010)), split = "")
seqIG2_3_2<-strsplit(c2s(substring(text = (seqsIG2_3[grep("s2", seqsIG2_3)]), 11, 5010)), split = "")
seqIG2_3_3<-strsplit(c2s(substring(text = (seqsIG2_3[grep("s3", seqsIG2_3)]), 11, 5010)), split = "")
seqIG2_3_4<-strsplit(c2s(substring(text = (seqsIG2_3[grep("s4", seqsIG2_3)]), 11, 5010)), split = "")

#Bind sequences into a matrix and convert to DNAbin format
dnabinIG2_3<-as.DNAbin(rbind(seqIG2_3_1, seqIG2_3_2, seqIG2_3_3, seqIG2_3_4)) 

#calculate distance matrix
dist_matIG2_3<-as.matrix(dist.dna(dnabinIG2_3))

#Extract distances from matrixs
S1ig2_3<-dist_matIG2_3[2, 3]
S2ig2_3<-dist_matIG2_3[1, 2]
S3ig2_3<-dist_matIG2_3[1, 3]

#Now plot all of the distances (2 to 3 IG)
#plot(c(S1sp, S1ig2_3, S2sp, S2ig2_3, S3sp, S3ig2_3))

return(c(S1sp, S2sp, S3sp, S1ig, S2ig, S3ig, S1ig2_3, S2ig2_3, S3ig2_3))
}

sim_rep_out<-replicate(1000, full_sim())

stats_out_df<-data.frame(matrix(unlist(sim_rep_out), nrow=1000, byrow=TRUE))

names(stats_out_df)=c("S1sp", "S2sp", "S3sp", "S1igto2", "S2igto2", "S3igto2", "S1igto3", "S2igto3", "S3igto3")


#Create boxplot for 3 to 2 
boxplot(stats_out_df$S1sp, stats_out_df$S1igto2, 
        stats_out_df$S2sp, stats_out_df$S2igto2, 
        stats_out_df$S3sp, stats_out_df$S3igto2, outline=FALSE)


#Stats for 3 to 2 boxplot
wilcox.test(stats_out_df$S1sp, stats_out_df$S1igto2)
wilcox.test(stats_out_df$S1sp, stats_out_df$S2igto2)
wilcox.test(stats_out_df$S1sp, stats_out_df$S3igto2)

#Create a boxplot for 2 to 3
boxplot(stats_out_df$S1sp, stats_out_df$S1igto3, 
        stats_out_df$S2sp, stats_out_df$S2igto3, 
        stats_out_df$S3sp, stats_out_df$S3igto3, outline=FALSE)

#stats for 2 to 3
wilcox.test(stats_out_df$S1sp, stats_out_df$S1igto3)
wilcox.test(stats_out_df$S2sp, stats_out_df$S2igto3)
wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto3)


### BI-DIRECTIONAL
#create a boxplot where 500 reps are from 3 to 2 and 500 reps are from 2 to 3
boxplot(stats_out_df$S1sp, c(sample(stats_out_df$S1igto2, 500),sample(stats_out_df$S1igto3, 500)), 
        stats_out_df$S2sp, c(sample(stats_out_df$S2igto2, 500),sample(stats_out_df$S2igto3, 500)), 
        stats_out_df$S3sp, c(sample(stats_out_df$S3igto2, 500),sample(stats_out_df$S3igto3, 500)),
        outline=FALSE)

wilcox.test(stats_out_df$S1sp, c(sample(stats_out_df$S1igto2, 500),sample(stats_out_df$S1igto3, 500)))
wilcox.test(stats_out_df$S2sp, c(sample(stats_out_df$S2igto2, 500),sample(stats_out_df$S2igto3, 500)))
wilcox.test(stats_out_df$S3sp, c(sample(stats_out_df$S3igto2, 500),sample(stats_out_df$S3igto3, 500)))



install.packages("ggplot2")
library(ggplot2)

means<-c(mean(stats_out_df$S1sp), mean(stats_out_df$S1igto2), mean(stats_out_df$S2sp), mean(stats_out_df$S2igto2), mean(stats_out_df$S3sp), mean(stats_out_df$S3igto2))
SDs<-c(sd(stats_out_df$S1sp), sd(stats_out_df$S1igto2), sd(stats_out_df$S2sp), sd(stats_out_df$S2igto2), sd(stats_out_df$S3sp), sd(stats_out_df$S3igto2))
maxs<-means+(2*SDs)
mins<-means-(2*SDs)

plot_data<-rbind.data.frame(means, SDs, maxs, mins)
names(plot_data)<-c("S1sp", "S1ig", "S2sp", "S2ig", "S3sp", "s3ig")
rownames(plot_data)<-c("mean", "SD", "max", "min")

plot_data_t<-t(plot_data)

plot_data_t_df<-as.data.frame(plot_data_t)

ggplot(plot_data_t_df, aes(x=rownames(plot_data_t_df), y=mean))+ geom_point() + 
  geom_errorbar(width=.1, aes(ymin=min, ymax=max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



#Plot means of 3 to 2 IG simulation
plot(c(mean(stats_out_df$S1sp), mean(stats_out_df$S1igto2), mean(stats_out_df$S2sp), mean(stats_out_df$S2igto2), mean(stats_out_df$S3sp), mean(stats_out_df$S3igto2)))


#Plot means of 2 to 3 IG simulation
plot(c(mean(stats_out_df$S1sp), mean(stats_out_df$S1igto3), mean(stats_out_df$S2sp), mean(stats_out_df$S2igto3), mean(stats_out_df$S3sp), mean(stats_out_df$S3igto3)))

