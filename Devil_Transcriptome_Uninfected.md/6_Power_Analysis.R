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
Sex<-factor(c("F","M","F","M","F","F","M","F","M","F","F","M","M","F","M","F","M","M","M","F"))
DGE_Object<-DGEList(counts=Devil_Data,group=Sex)

# Define your sample IDs, in this case we stored in object Sample)IDs
Sample_Names<-c("1-1","1-2","1-3","1-4","2-1","2-3","2-4","2-5","2-6","3-1","3-3","3-6","4-2","4-3","4-4","4-5","4-6","5-1","5-3","5-4")
colnames(DGE_Object)<-Sample_Names

#Filter out genes not expressed in all individuals then create the model matrix or the design for what you are testing
filtered_genes<-rowSums(cpm(DGE_Object)>0) > 0
DGE_Object<-DGE_Object[filtered_genes, , keep.lib.sizes=FALSE]
test<-DGE_Object$samples$group
mm<-model.matrix(~test)

## Try going through this step by step like you would for EdgeR
#glmQLFit and glmQLFTest functions, which are alternatives to glmFit and glmLRT 
#replace the chisquare approximation to the likelihood ratio statistic with a quasi-likelihood F-test
#resulting in more conservative and rigorous type I error rate control methods
MyMethod<-
    function(count.matrix, mm, column=2, group=NULL) {
       if (is.null(group)) {
            group = rep(1, ncol(count.matrix))
        }
        mm = mm
        d = edgeR::DGEList(counts=count.matrix, group=group)
        d = edgeR::calcNormFactors(d)
        d = edgeR::estimateDisp(d, mm)
        fit = edgeR::glmQLFit(d,mm,dispersion=d$common.dispersion,robust=TRUE)
        lrt = edgeR::glmQLFTest(fit)
        ret = lrt$table[, c("logFC", "PValue")]
        colnames(ret) = c("coefficient", "pvalue")
        ret
    }

### Define the proportion of the provided data you want to subsample
# Then define the methods you are interested in comparing, in this case we only used EdgeR
# Make sure to specify the number of points in your interval you want to map (length.out=20, so 20)
proportions<-10^seq(-2,0,length.out=10)
subsamples<-subsample(Devil_Data,proportions,method=c(MyMethod),mm=mm)
subsamples.summary<-summary(subsamples,FDR.level=0.05,p.adjust.method="fdr")
subsamples.summary

### Summarize the power analysis output for each proportion of read depth in csv file
write.csv(subsamples.summary,"Subseq_output.csv")

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


### We wanted to customize some settings a bit
## So, here is my code for my figure
pdf("Alt_Fig_2.pdf")
par(mfrow=c(2,2))
ggplot(subsamples.summary,aes(x=depth,y=significant,col=method)) + geom_line() + theme(axis.title.x=element_text(size=17,color='black')) + theme(axis.title.y=element_text(size=17,color='black')) + theme(axis.text.x=element_text(size=15,color='black')) + theme(axis.text.y=element_text(size=15,color='black')) + labs(x="Sequencing Depth",y="Number of Significantly Differentially Expressed Genes")
ggplot(subsamples.summary,aes(x=depth,y=pearson,col=method)) + geom_line() + theme(axis.title.x=element_text(size=17,color='black')) + theme(axis.title.y=element_text(size=17,color='black')) + theme(axis.text.x=element_text(size=15,color='black')) + theme(axis.text.y=element_text(size=15,color='black')) + labs(x="Sequencing Depth",y="Pearson's Correlate")
ggplot(subsamples.summary,aes(x=depth,y=MSE,col=method)) + geom_line() + theme(axis.title.x=element_text(size=17,color='black')) + theme(axis.title.y=element_text(size=17,color='black')) + theme(axis.text.x=element_text(size=15,color='black')) + theme(axis.text.y=element_text(size=15,color='black')) + labs(x="Sequencing Depth",y="Mean-squared error (MSE)")
ggplot(subsamples.summary,aes(x=depth,y=rFDP,col=method)) + geom_line() + theme(axis.title.x=element_text(size=17,color='black')) + theme(axis.title.y=element_text(size=17,color='black')) + theme(axis.text.x=element_text(size=15,color='black')) + theme(axis.text.y=element_text(size=15,color='black')) + labs(x="Sequencing Depth",y="Relative false discovery proportion (rFDP)")
dev.off()
