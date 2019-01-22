
###########################################################################################################################
# 6. Run PCAdapt (Outlier test) with selected K values
###########################################################################################################################

## Take the output from the sript you wrote to produce the PCAdapt input files 
 # Then set the working directory and load the required package
setwd("/Plink_Input")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)
library(pcadapt)

## Read in input files for choosing K for the pre-DFTD populations
 # Choosing K for the pre-DFTD populations
pre<-read.pcadapt("/Plink_Input/ped_pre_disease.ped",type="ped")
data<-pcadapt(pre,K=20)
pdf(file="/PCAdapt/Pre_DFTD_K1_20.pdf")
KCH<-plot(data,option="screeplot")
dev.off()

## Once you have chosen a K value then plug it into the equation below
 # our K value was optimized at K=2
K2_Pre<-pcadapt(pre,K=2)
pdf(file="/PCAdapt/K2_Pre")
summary(Pre_DFTD_K2)
dev.off()

## Outliers are deemed significant if the q values < alpha value (significance threshold)
qval<-qvalue(K2_Pre$pvalues)$qvalues
alpha<-0.01
outliers<-which(qval<alpha)
write.table(outliers,file="/PCAdapt/K2_Pre.csv",quote=FALSE)
write.csv(qval,file="/PCAdapt/K2_Pre_All.csv",quote=FALSE)

#######################################################################################################
## Read in input files for choosing K for the post-DFTD populations
post<-read.pcadapt("/Plink_Input/ped_post_disease.ped,type="ped")
data<-pcadapt(post,K=20)
pdf(file="/PCAdapt/Post_DFTD_K1_20.pdf")
KCH<-plot(data,option="screeplot")
dev.off()

## Once you have chosen a K value then plug it into the equation below
 # our K value was optimized at K=2
K2_Post<-pcadapt(post,K=2)
pdf(file="/PCAdapt/K2_Post")
summary(Post_DFTD_K2)
dev.off()

## Outliers are deemed significant if the q values < alpha value (significance threshold)
qval<-qvalue(K2_post$pvalues)$qvalues
alpha<-0.01
outliers<-which(qval<alpha)
write.table(outliers,file="/PCAdapt/K2_Post.csv",quote=FALSE)
write.csv(qval,file="/PCAdapt/K2_Post_All.csv",quote=FALSE)
