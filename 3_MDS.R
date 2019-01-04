
###########################################################################################################################
# 3. Multidimenional Scaling Analysis
###########################################################################################################################

## Navigate into the directory your gene count matrix is in
cd $COUNTDIR
module load r
R

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
