###########################################################################################################################
# 5. Weighted Gene Co-Expression Network Analysis
###########################################################################################################################

## Navigate into the directory your gene count matrix is in
#cd $COUNTDIR
#module load r
#R

## Or download your gene count matrix file onto your local computer and read it in as a csv file 
# adjust the working directory in the read.csv line to include the path to your locally downloaded file
#install.packages("BiocManager") 
#BiocManager::install("WGCNA") 

### Load the library containing the required packages for EdgeR & WGCNA
library(WGCNA)
library(flashClust)
enableWGCNAThreads()
library('edgeR')
options(stringsAsFactors = FALSE)

### Read in the csv file containing your gene count matrix
# Fill in $LIST with the "groupings" information or the factors your are comparing expression against
Devil_Data<-read.csv("gene_count_matrix.csv",row.names=1)
DGE_Object<-DGEList(counts=Devil_Data)

# Use this to check to make sure that your sample information is correct and that your library sizes were included
# Then we filtered out genes not expressed in any individuals
filtered_genes<-rowSums(cpm(DGE_Object)>0) > 0
DGE_Object<-DGE_Object[filtered_genes, , keep.lib.sizes=FALSE]

### Normalize by library size for RNA composition by finding scaling factors for library sizes that minimize the log-fold changes 
# between the samples- kind of like an effective library size
# Log transformed the library-size normalzied gene expression count data converting output to cpms
# Finally Transpsoed data frame so it would be appropriate for required WGCNA input
DGE_Object<-calcNormFactors(DGE_Object)
data.frame<-cpm(DGE_Object,normalized.lib.sizes=TRUE,log=TRUE,prior.count=1)
mydata.trans<-as.data.frame(t(Devil_Data))

### Formatted the transformed data frame to include the sample ID information encompassed in our LIST variable
# Also read in our phenotypic data which was a csv file containing the grouping factors of interest, population and sex
# data for each sample 
# The csv file was formatted so that the first column contained the sample IDs and the second and third grouping factor data 
rownames(mydata.trans)<-c(LIST)
mydata.trans<-as.data.frame(mydata.trans)
pheno_data<-read.csv("Phenotype_Data.csv")
d1<-data.frame(lapply(mydata.trans,function(x) as.numeric(as.character(x))),check.names=F,row.names=rownames(mydata.trans))

#Next, we removed all the genes that have 0 gene copies
gsg<-goodSamplesGenes(d1,verbose=3)
gsg$allOK
if (!gsg$allOK)
{
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:",paste(names(d1)[!gsg$goodGenes],collapse=", ")))
mydata<-d1[gsg$goodSamples,gsg$goodGenes]
}

#Then we clustered our samples to see if there were any outliers to filter out and generate the soft threshold values
sampleTree<-hclust(dist(mydata),method="average")
powers<-c(c(1:10),seq(from=12,to=20,by=2))
sft<-pickSoftThreshold(mydata,powerVector=powers,verbose=5)

# Called the network topology analysis function and plotted the scale-free topology fit index as a function of the soft-thresholding power
pdf(file="soft_thresholding_power.pdf")
cex1=0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main=paste("Scale independence"),ylim=c(-0.5,0.9))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.7,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,cex=0.9,col="red")
dev.off()

### At this point, we stopped to look at the figure to see derive the soft thresholding power, which was the
# lowest power for which the scale-free topology fit index reaches 0.70 
softPower<-7
adjacency<-adjacency(mydata,power=softPower)

### To minimize effects of noise and spurious associations transform the adjacency into TOM
# TOM: Topological Overlap Matrix is used to minimize effects of noise and spurious associations, 
# we transform the adjacency into and calculate the corresponding dissimilarity
TOM<-TOMsimilarity(adjacency)
dissTOM<-1-TOM
geneTree<-hclust(as.dist(dissTOM),method="average")
pdf(file="TOM clustering tree.pdf")
plot(geneTree,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity",labels=FALSE,hang=0.04)
dev.off()

### Then we generate and color the modules and set the minimum module size relatively high
minModuleSize<-30
dynamicMods<-cutreeDynamic(dendro=geneTree,distM=dissTOM,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=minModuleSize)
table(dynamicMods)
dynamicColors<-labels2colors(dynamicMods)
d2<-data.frame(lapply(mydata,function(x) as.numeric(as.character(x))),check.names=F,row.names=rownames(mydata))

### Figure showing dendogram colors with the colors and labels beneath it
sizeGrWindow(8,6)
pdf(file="Module_Dendogram.pdf")
plotDendroAndColors(geneTree,dynamicColors,"Dynamic Tree Cut",dendroLabels=FALSE,hang=0.03,addGuide=TRUE,
guideHang=0.05,main="Gene dendrogram and module colors")
dev.off()

#### Calculate eigengenes values for the modules and generate a tree to look at clustering patterns
# THen plot the tree output 
MEList<-moduleEigengenes(d2,colors=dynamicColors)
MEs<-MEList$eigengenes
MEDiss<-1-cor(MEs)
METree<-hclust(as.dist(MEDiss),method="average")
sizeGrWindow(7, 6)
pdf(file="Plot of Clustered Eigengenes.pdf")
plot(METree,main="Clustering of module eigengenes",xlab="",sub="")
MEDissThres<-0.25
abline(h=MEDissThres,col="red")
dev.off()

### In order to see what the merging did to the module colors we plot the gene dendogram again
# Call an automatic merging function
merge<-mergeCloseModules(d2,dynamicColors,cutHeight=MEDissThres,verbose=3)
mergedColors<-merge$colors
mergedMEs<-merge$newMEs
sizeGrWindow(12, 9)
pdf(file="Merged_Gene_Dendrogram.pdf",wi=9,he=6)
plotDendroAndColors(geneTree,cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05)
dev.off()

### For analysis of the relation of these modules to specific traits of interest
# i.e. sex and geography (our grouping factors)
moduleColors<-mergedColors
colorOrder<-c("grey",standardColors(50))
moduleLabels<-match(moduleColors,colorOrder)-1
MEs<-mergedMEs
MEs0<-moduleEigengenes(mydata,moduleColors)$eigengenes
rownames(MEs0)<-rownames(d2)

##Quantifying module-trait associations
#Define numbers of genes and samples
nGenes<-ncol(d2)
nSamples<-nrow(d2)
d3<-data.frame(lapply(MEs0,function(x) as.numeric(as.character(x))),check.names=F,row.names=rownames(MEs0))

### We want to create a new file in which we have bound our recoded population values into our data traits
# So we want to create columns in which each sample that is of a certain value receives a 1 or 0
# i.e: sex: F = 1, M = 0
# i.e: geography: individual from ARV=1, FRY or WPP = 0
# Recode your variables accordingly, but we are using population and sex
# We named the variable containing the dataframe testTraits 
test<-cbind(pheno_data,ARV,FRY,WPP)
testTraits<-test[,-c(1,2,4,6)]
colnames(testTraits)<-c('ARV','FRY','WPP','Sex')
datTraits<-data.frame(lapply(testTraits,function(x) as.numeric(as.character(x))),check.names=F,row.names=rownames(testTraits))

### Sample dendogram with the heat map of expression values based on grouping facotrs (i.e: sex and geography) 
# Start by reclustering the samples then convert the traits to a color representation
sampleTree2<-hclust(dist(mydata),method="average")
traitColors<-numbers2colors(datTraits,signed=FALSE)
pdf("Dendogram_Heat_Map.pdf")
plotDendroAndColors(sampleTree2,traitColors,groupLabels=names(datTraits),main="Sample dendrogram and trait heatmap")
dev.off()

##### The following code looks for gene modules significantly associated with one of our grouping 
# factors of interest, in this case sex 
# This code was also rerun replacing sex with ARV, FRY and WPP (our populations of interest)
MEs<-orderMEs(d3)
moduleTraitCor<-cor(MEs,sex,use="p")
moduleTraitPvalue<-corPvalueStudent(moduleTraitCor,nSamples)

### Figure that will generate a heat map looking at the correlations between each grouping factor and module
pdf(file="Module_Trait_Relationship_Sex.pdf")
textMatrix<-paste(signif(moduleTraitCor,2),"\n(",signif(moduleTraitPvalue, 1),")",sep="")
dim(textMatrix)<-dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix=moduleTraitCor,xLabels=names(FRY),yLabels=names(MEs),ySymbols=names(MEs),
colorLabels=FALSE,colors=blueWhiteRed(50),textMatrix=textMatrix,setStdMargins=FALSE,cex.text=0.5,
zlim=c(-1,1),main=paste("Module-trait relationships"))
dev.off()
