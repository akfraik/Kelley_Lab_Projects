###########################################################################################################################
# 4. EdgeR Differential Expression Analysis
###########################################################################################################################

## Navigate into the directory your gene count matrix is in
#cd $COUNTDIR
#module load r
#R

## Or download your gene count matrix file onto your local computer and read it in as a csv file 
# adjust the working directory in the read.csv line to include the path to your locally downloaded file
## Load the library containing the required packages for EdgeR
library('edgeR')
library('statmod')
library('limma')

## Read in the csv file containing your gene count matrix
# Fill in $LIST with the "groupings" information or the factors your are comparing expression against
x<-read.csv("gene_count_matrix.csv",row.names=1)
groupings<-factor(c("FRY_F","FRY_M","WPP_F","WPP_M","FRY_F","WPP_F","WPP_M","ARV_F","ARV_M","FRY_F","WPP_F","ARV_M","FRY_M","WPP_F","WPP_M","ARV_F","ARV_M","FRY_M","WPP_M","ARV_F"))
DGE_Object<-DGEList(counts=x,group=groupings)

# Use this to check to make sure that your sample information is correct
# and that your library sizes were included
# Then we filtered out genes not expressed in any individuals
DGE_Object$samples
filtered_genes<-rowSums(cpm(DGE_Object)>0) > 0
DGE_Object<-DGE_Object[filtered_genes, , keep.lib.sizes=FALSE]

### Normalize by library size for RNA composition by finding scaling factors for library sizes that minimize the log-fold changes 
# between the samples- kind of like an effective library size
DGE_Object<-calcNormFactors(DGE_Object)
DGE_Object$samples

#Now, for the MDS plots Reprint all of the sample labels
Sample_Names<-c("1-1","1-2","1-3","1-4","2-1","2-3","2-4","2-5","2-6","3-1","3-3","3-6","4-2","4-3","4-4","4-5","4-6","5-1","5-3","5-4")

#Color by population
colors<-c('blue','blue','hotpink','hotpink','limegreen','limegreen')

#Make points by sex
points<-c(1,2,1,2,1,2)

#Now set up the design matrix in which you are establishing all
#of your variables of interest
Population<-factor(c("FRY","FRY","WPP","WPP","FRY","WPP","WPP","ARV","ARV","FRY","WPP","ARV","FRY","WPP","WPP","ARV","ARV","FRY","WPP","ARV"))
Sex<-factor(c("F","M","F","M","F","F","M","F","M","F","F","M","M","F","M","F","M","M","M","F"))
data.frame(Sample=colnames(DGE_Object),Population,Sex)

### Create the model matrix where you are using population as the replicate and the 
# sex as the variable we are interested in our additive linear model using "females" as our null or group of comparison
### Model matrix starting point- more models can be tested using the following code
# Just replacing the "design" varaible with a new model to look at DE patterns between other grouping factors of interest such as 
# design 2 (line 86), design 3 (line 87), design 4 (line 88)
design<-model.matrix(~Sex)
rownames(design)<-colnames(DGE_Object)

#####For experiments with multiple factors, edgeR uses Cox-Reid adjusted
#likelihood methods in estimating dispersion
#This produces both the common and tagwise dispersion parameters for comaprison
DGE_Object<-estimateDisp(DGE_Object,design)
DGE_Object$common.dispersion
BCV<-sqrt(DGE_Object$common.dispersion)

#Write BCV plot to output not blocking by population
jpeg("BCV_plot_design.jpg")
plotBCV(DGE_Object)
dev.off()

#plot the fit of the QL dispersions around the dispersion trend
fit<-glmQLFit(DGE_Object,design,robust=TRUE)
jpeg("QLDplot_design.jpg")
plotQLDisp(fit)
dev.off()

#Test the output to find the top 10 DE genes between sex
#Positive LF changes means higher expression in M {coefficient SexM}
test<-glmQLFTest(fit)
design_output<-topTags(test,n=10)
summary(decideTests(test))
write.csv(design_output,"design_output.csv")

### We also tested the following model where we ran population as the blocked effect
# The following models can be run using the above code by replacing
# "design" with your desired model
# In our case we substitituted "design" with design_2, design_3 and design_4 using the above code
# design_2: DE with additive model looking at differences between sex while blocking for population
# design_3: DE with additive linear model with 0 intercept - no "null" group of comparison
# design_4: DE with additive model looking at differences between population while blocking for sex
design_2<-model.matrix(~Population+Sex)
design_3<-model.matrix(~0+Population)
design_4<-model.matrix(~Sex+Population)
