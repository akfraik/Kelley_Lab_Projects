###########################################################################################################################
# 6. Power Analysis (SubSeq)
###########################################################################################################################

## Navigate into the directory your gene count matrix is in
#cd $COUNTDIR
#module load r
#R

### Or download your gene count matrix file onto your local computer and read it in as a csv file 
# adjust the working directory in the read.csv line to include the path to your locally downloaded file
# install_github("jdstorey/qvalue")
# install_github("StoreyLab/subSeq")
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("limma","edgeR","DESeq2","DEXSeq","pasilla"))

### Load the library containing the required libraries for Subseq and figure generation
library(devtools)
library(edgeR)
library(subSeq)
library(ggplot2)

### Read in the gene count matrix 
# Define your grouping factors of interest (i.e: numeric coding for sex, 1 = female, 2 = male stored in object LIST)
Devil_Data<-read.csv("gene_count_matrix.csv",row.names=1)
groupings<-factor(c(1,2,1,2,1,1,2,1,2,1,1,2,2,1,2,1,2,2,2,1))
DGE_Object<-DGEList(counts=Devil_Data,group=groupings)

# Define your sample IDs, in this case we stored in object Sample)IDs
Sample_Names<-c(Sample_IDs)
colnames(DGE_Object)<-Sample_Names

### Define the proportion of the provided data you want to subsample
# Then define the methods you are interested in comparing, in this case we only used EdgeR
# Make sure to specify the number of points in your interval you want to map (length.out=20, so 20)
proportions<-10^seq(-2,0,length.out=20)
subsamples<-subsample(DGE_Object,proportions,method=c("edgeR","DESeq2"),treatment=DGE_Object$samples$group)
options(width=40)
subsamples

### Summarize the power analysis output for each proportion of read depth in csv file
subsamples.summary<-summary(subsamples)
subsamples.summary
write.csv(subsamples.summary,"subsamples_summary_12_18.csv")

### Create a four paneled figure showing the relationship between sequencing depth and
# 1: number of significant DE genes detected
# 2: strength of pearson's correlate value (grouping factor VS cpm)
# 3: MSE
# 4: rFDP or rate of false discovery
pdf("Figure_subsamples_1.pdf")
par(mfrow=c(2,2))
ggplot(subsamples.summary,aes(x=depth,y=significant,col=method)) + geom_line()
ggplot(subsamples.summary,aes(x=depth,y=pearson,col=method)) + geom_line()
ggplot(subsamples.summary,aes(x=depth,y=MSE,col=method)) + geom_line()
ggplot(subsamples.summary,aes(x=depth,y=rFDP,col=method)) + geom_line()
dev.off()


### We focused on number of significant DE genes and pearon's correlate in our figures for the MS
# Panel a
pdf("Power_Analysis_Figure_Panel_1.pdf")
p<-ggplot(subsamples.summary,aes(x=depth,y=significant,color='black'))
p<-p+geom_line(size=1.5)
p<-p+theme_classic()
p<-p+theme(panel.grid.major=element_blank())
p<-p+theme(panel.grid.minor=element_blank())
p<-p+theme(axis.line=element_line(size=1,color='black'))
p<-p+theme(axis.ticks=element_line(color='black'))
p<-p+theme(axis.title.x=element_text(size=17,color='black'))
p<-p+theme(axis.title.y=element_text(size=17,color='black'))
p<-p+theme(axis.text.x=element_text(size=15,color='black'))
p<-p+theme(axis.text.y=element_text(size=15,color='black'))
p<-p+labs(x="Sequencing Depth",y="Number of Significantly Differentially Expressed Genes")
p
dev.off()

# Panel b 
pdf("Power_Analysis_Figure_Panel_2.pdf")
p<-ggplot(subsamples.summary,aes(x=depth,y=pearson,col=method))
p<-p+geom_line(size=1.5)
p<-p+theme_classic()
p<-p+theme(panel.grid.major=element_blank())
p<-p+theme(panel.grid.minor=element_blank())
p<-p+theme(axis.line=element_line(size=1,color='black'))
p<-p+theme(axis.ticks=element_line(color='black'))
p<-p+theme(axis.title.x=element_text(size=17,color='black'))
p<-p+theme(axis.title.y=element_text(size=17,color='black'))
p<-p+theme(axis.text.x=element_text(size=15,color='black'))
p<-p+theme(axis.text.y=element_text(size=15,color='black'))
p<-p+labs(x="Sequencing Depth",y="Pearson's Correlate")
p
dev.off()


