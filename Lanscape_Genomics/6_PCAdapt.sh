#!/usr/bin/env Rscript
#R

## Comment this line out if you run this script more then once
mkdir PCAdapt

## Take the output from the sript you wrote to produce the PCAdapt input files 
 # Then set the working directory and load the required package
setwd("/Plink_Input")
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)
library(pcadapt)

## Read in input files
 # Choosing K
x<-"beagle_rap_pre.ped"
pre<-read.pcadapt(x,type="ped")
data<-pcadapt(pre,K=20)

## Saving plot in pdf
pdf(file="/PCAdapt/Pre_DFTD_K1_20.pdf")
KCH<-plot(data,option="screeplot")
dev.off()

## Once you have chosen a K value then plug it into the equation below
 # our K value was optimized at K=2
K2<-pcadapt(pre,K=2)
pdf(file="/PCAdapt/Rap_2_Pre")
summary(Rap_2)
dev.off()
 # Outliers are deemed significant if the q values < alpha value (significance threshold)
qval<-qvalue(Rap_2$pvalues)$qvalues
alpha<-0.01
outliers<-which(qval<alpha)
write.table(outliers,file="/PCAdapt/K_2_Pre.csv",quote=FALSE)

#K=6
Rap_6<-pcadapt(pre,K=6)
pdf(file="/project/devils_alexf/try3/results/Rapture/PCAdapt/Rap_6_Pre")
summary(Rap_6)
dev.off()

#outliers have to be less then alpha values
qval<-qvalue(Rap_6$pvalues)$qvalues
alpha<-0.01
outliers<-which(qval<alpha)
write.table(outliers,file="/project/devils_alexf/try3/results/Rapture/PCAdapt/Rap_6_Pre.csv",quote=FALSE)
write.csv(qval,file="/project/devils_alexf/try3/results/Rapture/PCAdapt/Rap_Pre_All.csv",quote=FALSE)

#K=4
Rap_4<-pcadapt(pre,K=4)
pdf(file="/project/devils_alexf/try3/results/Rapture/PCAdapt/Rap_4_Pre")
summary(Rap_4)
dev.off()

#outliers have to be less then alpha values
qval<-qvalue(Rap_4$pvalues)$qvalues
alpha<-0.01
outliers<-which(qval<alpha)
write.table(outliers,file="/PCAdapt/Rap_4_Pre.csv",quote=FALSE)
