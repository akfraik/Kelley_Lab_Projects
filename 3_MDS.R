
###########################################################################################################################
# 3. Multidimenional Scaling Analysis
###########################################################################################################################

## Navigate into the directory your gene count matrix is in
#cd $COUNTDIR
#module load r
#R

## Or download your gene count matrix file onto your local computer and read it in as a csv file 
# adjust the working directory in the read.csv line to include the path to your locally downloaded file

## Install the necessary packages and then comment them out 
# install.packages("edgeR")
# source("https://www.bioconductor.org/biocLite.R")
# BiocInstaller::biocLite(packages)
# biocLite('edgeR') 

## Load the library containing the required packages for MDS 
library('edgeR')

## Read in the csv file containing your gene count matrix
# Fill in $LIST with the "groupings" information or the factors your are comparing expression against
x<-read.csv("gene_count_matrix.csv",row.names=1)
groupings<-factor(c($LIST))
DGE_Object<-DGEList(counts=x,group=groupings)

## Filters out individuals in which genes are not expressed in at least 2 individuals
# and those genes with a gene count out less than 6
filtered_genes<-rowSums(cpm(DGE_Object)>1) >= 2
DGE_Object<-DGE_Object[filtered_genes, , keep.lib.sizes=FALSE]
DGE_Object<-calcNormFactors(DGE_Object)
DGE_Object$samples

# Create a dataframe with the normalized library counts 
data.frame<-cpm(DGE_Object,normalized.lib.sizes=TRUE,log=TRUE,prior.count=1)

## Plot the MDS of the top 1,000 genes or the top genes differentially expressed across all samples
# color by your groupings factor and create labels with you sample information
# Sample_Names <- Names of your samples (i.e: 1:20)
# colors <- colors for grouping factors
# pch <- "pch" or vector of values for the shape you want for a grouping factor
pdf(file="Top_1000_Sample_Names.pdf")
plotMDS(DGE_Object,top=1000,col=colors[groupings],labels=Sample_Names,gene.selection="common")
legend(3,1,legend=levels(groupings),cex=1,pch=points,col=colors,pt.bg=colors,ncol=1)
legend.text=TRUE
dev.off()

## Plot the MDS without the sample names, just have points with colors and shapes representing grouping factors
# use the top 500, 1000, 2000, 10000 genes and save a different jpeg
values<-c(500,1000,2000,10000)
for ( i in values ) {
jpeg(paste(i,'MDS_Plots_sn.jpg',sep= '_'))
plotMDS(DGE_Object,top=1000,col=colors[groupings],pch=points[groupings],bg=colors[groupings],gene.selection="common")
legend(3,1,legend=levels(groupings),cex=1,pch=points,col=colors,pt.bg=colors,ncol=1)
legend.text=TRUE
dev.off()
}
