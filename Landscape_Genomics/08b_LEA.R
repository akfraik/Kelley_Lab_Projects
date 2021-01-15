###########################################################################################################################
## 8b. LEA (R script)
###########################################################################################################################

## First step set the working directory and load the required package
 #source("https://bioconductor.org/biocLite.R")
 #biocLite("LEA")
 #browseVignettes("LEA")
 library(LEA)

## Read in file/convert the file from PED to LEA format
LEA_Variable<-ped2lfmm("ped_Variable_disease.ped","/LEA/Variable/ped_Variable_disease.lfmm")

## Perfom Tracy-Widom tests on all eigenvalues and create file: tuto.tracyWidom - tracy-widom test information.
 # display the p-values for the Tracy-Widom tests.
tw<-tracy.widom(LEA_Variable)
tw$pvalues[1:10]

## Plot the percentage of variance explained by each component
PCA_Plotted_Variable<-plot(tw$percentage)

## Store PCA output for Variable
pdf(file="/LEA/Variable/PCA_Variable.pdf")
plot(tw$percentage,xlim=c(0,15))
dev.off()

## Now, convert our PCA LEA file to a genotype file
 #Runs with K between 1 and 10 with cross-entropy and 10 repetitions.
GENO_Pre<-ped2geno("ped_Variable_disease.ped","/LEA/Variable/ped_Variable_disease.geno")
project<-snmf("/LEA/Variable/ped_Variable_disease.geno",K=1:10,entropy=TRUE,repetitions=10,project="new")

## Get the cross-entropy of each run for K = 6
 # then select the run with the lowest cross-entropy
ce<-cross.entropy(project,K=6)
best<-which.min(ce)

## Plot CE criterion of all runs of the project
pdf(file="/LEA/Variable/Figure_CE_Variable_disease.pdf")
plot(project,lwd=5,col="red",pch=1)
dev.off()

## Write the environmental data to a file with the main options, K: (the number of latent factors)
 #Runs with K = 6 and 5 repetitions of 6000 iterations including 3000 iterations for burnin.
C<-read.table("/XTX/Veg_Variable")
write.env(C,"/LEA/Variable/env_Variable.env")
project<-lfmm("/LEA/Variable/beagle_rap_pre.lfmm","env_Variable.env",K=6,repetitions=5,project="new")

## Get the zscores of each run for K = 6 and combine using the Stouffer method
zs_pre<-z.scores(project,K=6)
zs.stouffer<-apply(zs,MARGIN=1,median)

## A genomic inflation factor was computed and calculate adjusted p-values
lambda<-median(zs^2)/.456
cp.values<-pchisq(zs.stouffer^2/lambda,df=1,lower=FALSE)

## List of candidate loci
for (alpha in c(.05,.1,.15,.2)) {
+     # expected FDR
+     print(paste("expected FDR:", alpha))
+     L = length(cp.values)
+     # return a list of candidates with an expected FDR of alpha.
+     w = which(sort(cp.values) < alpha * (1:L) / L)
+     candidates = order(cp.values)[w]
+     # estimated FDR and True Positif
+     estimated.FDR = length(which(candidates <= 350))/length(candidates)
+     estimated.TP = length(which(candidates > 350))/50
+ print(paste("FDR:",estimated.FDR,"True Positive:",estimated.TP)) +}

## The project displays the list of candidate loci for the "project"
show(project)
summary(project)

# Get the z-scores, p-values and log transformed p-values for the 2nd run for K = 6
z<-z.scores(project,K=6,run=5)
p<-p.values(project,K=6,run=5)
mp<-mlog10p.values(project,K=6,run=5)

## Then create a document with all of these data
final_outlier_set_Variable<-cbind(project,z,p,mp)
write.csv(final_outlier_set_Variable,"/LEA/Variable/final_outlier_set_Variable.csv")
