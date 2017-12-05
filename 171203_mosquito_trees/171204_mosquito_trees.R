#R script for analysing the topologies of mosquito gene trees and calculated node depths 
#in order to impliment our test of IG directionality

setwd("/Users/esforsythe/Documents/Beilstiein_lab_research/BIOINFORMATICS/IG_Direction/IG_Direction_repo/171203_mosquito_trees") 

install.packages("ape")
library(ape)
install.packages("phangorn")
library(phangorn)
install.packages("phytools")
library(phytools)

#Load trees and store the total number of trees
trees<-read.tree("catfile_mosquitotrees_171204")

#IMPORTANT!!!!
#Use this switch to look at only autozomal loci!
#There are 228 X chrom loci and they're at the end of the "catfile"
trees<-trees[1:(length(trees)-228)]


#Assign one tree for testing
#tree<-trees[[500]]

#Root tree by outgroup
#root_tree<-root(tree, "O", resolve.root=TRUE, edgelabel=TRUE)


skeeter_depths<-function(tree){

#midpoint root the tree
#chose this method for rooting because it allows me to calculate of node depths involving species R
root_tree<-midpoint(tree, node.labels="support")

#Turn tree into chronogram
#Note that lambda=0 here. I may need to explore this parameter
chrono_tree<-chronopl(root_tree, lambda = 0)

#now drop extra tips from tree
drop_tree<-drop.tip(chrono_tree, c("A", "G", "C"))

#Now return the topology of the tree
if(is.monophyletic(drop_tree, c("L", "Q"))){
  top<-"LQ"
}else if(is.monophyletic(drop_tree, c("R", "Q"))){
    top<-"RQ"
}else if(is.monophyletic(drop_tree, c("R", "L"))){
    top<-"RL"
}else{
    top<-"NA"
  }

#Get node depths for trees
D2<-1-node.depth.edgelength(drop_tree)[getMRCA(phy=drop_tree, c("L", "Q"))]

D1<-1-node.depth.edgelength(drop_tree)[getMRCA(phy=drop_tree, c("R", "Q"))]

D3<-1-node.depth.edgelength(drop_tree)[getMRCA(phy=drop_tree, c("R", "L"))]

return(c(top, D1, D2, D3))

}

#############################
###### END FUNCTION  ########
#############################

#run the script
output<-lapply(trees, skeeter_depths) 

#clean up output
#convert the output from a list to a dataframe
output_df <- data.frame(matrix(unlist(output), nrow=length(trees), byrow=TRUE), stringsAsFactors = FALSE)
names(output_df) <- c("Tops", "D1", "D2" ,"D3")

#subset dataset
Sub_LQ<-subset(output_df, output_df$Top=="LQ")
Sub_RQ<-subset(output_df, output_df$Top=="RQ")
Sub_RL<-subset(output_df, output_df$Top=="RL")


#Box plot for "T2" nodes from each tree
boxplot(
  as.numeric(Sub_LQ$D2), as.numeric(Sub_RQ$D1), as.numeric(Sub_RL$D3)
)

#LQ vs RQ
wilcox.test(as.numeric(Sub_LQ$D2), as.numeric(Sub_RQ$D1))$p.value

#LQ vs RL
wilcox.test(as.numeric(Sub_LQ$D2), as.numeric(Sub_RL$D3))$p.value

#LQ vs RQ
wilcox.test(as.numeric(Sub_RL$D3), as.numeric(Sub_RQ$D1))$p.value

### NOW DO IG directionality test!
boxplot(
  as.numeric(Sub_LQ$D1), as.numeric(Sub_RQ$D1), 
  as.numeric(Sub_LQ$D2), as.numeric(Sub_RQ$D2), 
  as.numeric(Sub_LQ$D3), as.numeric(Sub_RQ$D3)
)

wilcox.test(as.numeric(Sub_LQ$D1), as.numeric(Sub_RQ$D1))$p.value
wilcox.test(as.numeric(Sub_LQ$D2), as.numeric(Sub_RQ$D2))$p.value
wilcox.test(as.numeric(Sub_LQ$D3), as.numeric(Sub_RQ$D3))$p.value
