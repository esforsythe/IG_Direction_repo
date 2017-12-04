#This is a script for simulating sequence evolution during speciation and introgression
#I'll aslo use this script to test my test for directionality of IG.

#Set the working directory
setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/IG_Direction/IG_Direction_repo/170830_ms_simulations")

install.packages("phyclust")
library(phyclust)
library(parallel)
library(plyr)
install.packages("seqinr")
library(seqinr)

#Start giant for loop in order to iterate through different values of time of IN (tIG)
#time of IG parameter values to test
tIG<-c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

#initialize an empty dataframe for output from the for loop
datalist = list()

### Beginning of FOR LOOP ###
for (x in 1:length(tIG)){

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
#note: changed last -ej from -ej 0.2 5 3 to -ej 0.1 5 3
#tIG<-0.1
#note: when running this in a loop, change the tIG[7] below to tIG[x]
#old opts: "-T -t 50 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -es 0.1 2 0.4 -ej 0.1 5 3 -r 5 5000"
opts_32<-paste("-T -t 50 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -es", tIG[7], "2 0.4 -ej", tIG[7], "5 3 -r 5 5000"  )
ret.msIG <- ms(nsam = 4, nreps = 1, opts=opts_32)

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

#note: when running this in a loop, change the tIG[7] below to tIG[x]

#Here I'll try to make the f = 0.6
#Old opts: "-T -t 50 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -es 0.1 3 0.4 -ej 0.1 5 2 -r 5 5000"
opts_23<-paste("-T -t 50 -I 4 1 1 1 1 -ej 1 2 1 -ej 2 3 1 -ej 3 4 1 -es", tIG[7], "3 0.4 -ej", tIG[7], "5 2 -r 5 5000"  )
ret.msIG2_3 <- ms(nsam = 4, nreps = 1, opts =opts_23 )

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

#Extract distances from matricies
S1ig2_3<-dist_matIG2_3[2, 3]
S2ig2_3<-dist_matIG2_3[1, 2]
S3ig2_3<-dist_matIG2_3[1, 3]

#Now plot all of the distances (2 to 3 IG)
#plot(c(S1sp, S1ig2_3, S2sp, S2ig2_3, S3sp, S3ig2_3))

return(c(S1sp, S2sp, S3sp, S1ig, S2ig, S3ig, S1ig2_3, S2ig2_3, S3ig2_3))
}
#set the number of replicates to use
reps<-10000
sim_rep_out<-replicate(reps, full_sim())

stats_out_df<-data.frame(matrix(unlist(sim_rep_out), nrow=reps, byrow=TRUE))

names(stats_out_df)=c("S1sp", "S2sp", "S3sp", "S1igto2", "S2igto2", "S3igto2", "S1igto3", "S2igto3", "S3igto3")

#### Vary the proportion that were introgressed in each direction
#Make a variable "prop3_2", which is the proportion of the introgressed genes
#that were introgressed from 3 to 2. The proportion of introgressed gene from 2 to 3 = 1-prop3_2

#initiate an empty matrix
matS2<-matrix(, ncol = 11, nrow = 9)
matS3<-matrix(, ncol = 11, nrow = 9)

prop3_2<-c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

for(y in 1:length(prop3_2)){
#total measurement from uni or bidirectional introgressed genes.
total_S2<-c(sample(stats_out_df$S2igto2, prop3_2[y]*reps),sample(stats_out_df$S2igto3, ((1-(prop3_2[y]))*reps)))
total_S3<-c(sample(stats_out_df$S3igto2, prop3_2[y]*reps),sample(stats_out_df$S3igto3, ((1-(prop3_2[y]))*reps)))

#subsample the totals to test different numbers of introgressed vs speciation genes
#different value of the proportion of the genome that was introgressed
propIG<-c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

#start a nested loop for proportion of genome introgressed
for(z in 1:length(propIG)){
#Assign the final dataset for speciation genes
sp_finalS2<-sample(stats_out_df$S2sp, ((1-(propIG[z]))*reps))
sp_finalS3<-sample(stats_out_df$S3sp, ((1-(propIG[z]))*reps))

#Assign the final dataset for IG genes
ig_finalS2<-sample(total_S2, propIG[z]*reps)
ig_finalS3<-sample(total_S3, propIG[z]*reps)

# 32 = -1
# 23 = 1

#Ask which hypothesis S2 is consistent with
if(wilcox.test(sp_finalS2,ig_finalS2)$p.value < 0.05 &
   median(sp_finalS2) < median(ig_finalS2)){
  S2_concl<-(-1)
}else if(wilcox.test(sp_finalS2,ig_finalS2)$p.value < 0.05 &
         median(sp_finalS2) > median(ig_finalS2)) {
  S2_concl<-0
}else{S2_concl<-1}

#Ask which hypothesis S3 is consistent with
if(wilcox.test(sp_finalS3,ig_finalS3)$p.value < 0.05 &
   median(sp_finalS3) < median(ig_finalS3)){
  S3_concl<-0
}else if(wilcox.test(sp_finalS3,ig_finalS3)$p.value < 0.05 &
         median(sp_finalS3) > median(ig_finalS3)) {
  S3_concl<-1
}else{S3_concl<-(-1)}

#store results to a matrix of results
matS2[z,y]<-S2_concl
matS3[z,y]<-S3_concl
}
}
#clean up S2 and S3 matrix
row.names(matS2)<-propIG
colnames(matS2)<-prop3_2

row.names(matS3)<-propIG
colnames(matS3)<-prop3_2

#Now combine the two results
mat_comb<-matrix(, ncol = 11, nrow = 9)
for(y in 1:length(prop3_2)){
  for(z in 1:length(propIG)){
  if(matS2[z,y]==matS3[z,y]){
    mat_comb[z,y]<-matS2[z,y]
  }else{
    mat_comb[z,y]<-0
  }
  }
}

row.names(mat_comb)<-propIG
colnames(mat_comb)<-prop3_2

#install.packages(gplots)
#library(gplots)

mat_comb_new<-mat_comb[ nrow(mat_comb):1, ]

heatmap.2(mat_comb_new, Rowv =NULL, Colv = NULL, col=c("red", "white", "dark gray"),
          sepwidth=c(0.01,0.02),sepcolor="black",colsep=1:ncol(mat_comb_new),rowsep=1:nrow(mat_comb_new),
          key = FALSE, trace = "none", main = c(reps,"genes sampled"),
          xlab = "proportion of introgressed genes introgressed from 3 to 2",
          ylab = "proportion of sampled genes introgressed"
          )












### CREATE CODE FOR RESULTS ###
# 3 to 2
#Assign an object called "S2_igto2"
#32= consistent with IG from 3 to 2; 23=consistent with IG from 2 to 3; 0= not consistent with either. 
if(wilcox.test(stats_out_df$S2sp, stats_out_df$S2igto2)$p.value < 0.05 &
   median(stats_out_df$S2sp) < median(stats_out_df$S2igto2)){
     S2_igto2<-32
   }else if(wilcox.test(stats_out_df$S2sp, stats_out_df$S2igto2)$p.value < 0.05 &
        median(stats_out_df$S2sp) > median(stats_out_df$S2igto2)) {
          S2_igto2<-0
        }else{S2_igto2<-23}

#Assign an object called "S3_igto2"
#32= consistent with IG from 3 to 2; 23=consistent with IG from 2 to 3; 0= not consistent with either. 
if(wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto2)$p.value >= 0.05){
  S3_igto2<-32
}else if(wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto2)$p.value < 0.05 &
         median(stats_out_df$S3sp) > median(stats_out_df$S2igto2)) {
  S3_igto2<-23
}else{S3_igto2<-0}

# 2 to 3
#Assign an object called "S2_igto3"
#32= consistent with IG from 3 to 2; 23=consistent with IG from 2 to 3; 0= not consistent with either. 
if(wilcox.test(stats_out_df$S2sp, stats_out_df$S2igto3)$p.value < 0.05 &
   median(stats_out_df$S2sp) < median(stats_out_df$S2igto3)){
  S2_igto3<-32
}else if(wilcox.test(stats_out_df$S2sp, stats_out_df$S2igto3)$p.value < 0.05 &
         median(stats_out_df$S2sp) > median(stats_out_df$S2igto3)) {
  S2_igto3<-0
}else{S2_igto3<-23}

#Assign an object called "S3_igto3"
#32= consistent with IG from 3 to 2; 23=consistent with IG from 2 to 3; 0= not consistent with either. 
if(wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto3)$p.value >= 0.05){
  S3_igto3<-32
}else if(wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto3)$p.value < 0.05 &
         median(stats_out_df$S3sp) > median(stats_out_df$S2igto3)) {
  S3_igto3<-23
}else{S3_igto3<-0}

# Bidirectional IG #
#take half of 3 to 2 and half of 2 to 3
bidiS2<-c(sample(stats_out_df$S2igto2, 500),sample(stats_out_df$S2igto3, 500))
#Assign an object called "S2_bidi"
#32= consistent with IG from 3 to 2; 23=consistent with IG from 2 to 3; 0= not consistent with either. 
if(wilcox.test(stats_out_df$S2sp, bidiS2)$p.value < 0.05 &
   median(stats_out_df$S2sp) < median(bidiS2)){
  S2_bidi<-32
}else if(wilcox.test(stats_out_df$S2sp, bidiS2)$p.value < 0.05 &
         median(stats_out_df$S2sp) > median(bidiS2)) {
  S2_bidi<-0
}else{S2_bidi<-23}

#take half of 3 to 2 and half of 2 to 3
bidiS3<-c(sample(stats_out_df$S3igto2, 500),sample(stats_out_df$S3igto3, 500))
#Assign an object called "S3_bidi"
#32= consistent with IG from 3 to 2; 23=consistent with IG from 2 to 3; 0= not consistent with either. 
if(wilcox.test(stats_out_df$S3sp, bidiS3)$p.value >= 0.05){
  S3_bidi<-32
}else if(wilcox.test(stats_out_df$S3sp, bidiS3)$p.value < 0.05 &
         median(stats_out_df$S3sp) > median(bidiS3)) {
  S3_bidi<-23
}else{S3_bidi<-0}

#Return results (when added as a function)
#This isn't working
datalist[[x]]<-(c(tIG[x], S2_igto2, S2_igto3, S2_bidi, S3_igto2, S3_igto3, S3_bidi,
                  wilcox.test(stats_out_df$S2sp, stats_out_df$S2igto2)$p.value,
                  wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto2)$p.value,
                  wilcox.test(stats_out_df$S2sp, stats_out_df$S2igto3)$p.value,
                  wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto3)$p.value,
                  wilcox.test(stats_out_df$S2sp, bidiS2)$p.value,
                  wilcox.test(stats_out_df$S3sp, bidiS3)$p.value))

#### END OF GIANT FOR LOOP
}

#Combine results into dataframe
data_mat<-do.call(rbind, datalist)

colnames(data_mat)=c("IG_time", "S2_3to2", "S2_2to3", "S2_bidi", "S3_3to2", "S3_2to3", "S3_bidi", "p_S2_3to2", "p_S2_2to3", "p_S2_bidi", "p_S3_3to2", "p_S3_2to3", "p_S3_bidi")

data_framed<-as.data.frame(data_mat)





#Plot P values
install.packages("ggplot2")
library(ggplot2)

ggplot(data_framed, aes(x=IG_time, y=p_S2_3to2))+ geom_point() + scale_y_log10()
  geom_errorbar(width=.1, aes(ymin=min, ymax=max)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


class(data_framed)












   wilcox.test(stats_out_df$S3sp, stats_out_df$S3igto2)

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

