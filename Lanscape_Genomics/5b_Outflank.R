#!/usr/bin/env Rscript
#R

## Take the output from the sript you wrote to produce the outflank input files 
 # Then set the working directory and load the required package
setwd("/Outflank_Input")
library(OutFLANK)

## Read in the required snpmat, popnames and locus name files to make a datamatrix containing all of the data
snpmat<-read.csv('outflank_rap_pre.snpmat',sep='\t',as.is=TRUE,header=FALSE)
popnames<-scan('outflank_rap_pre.popnames',what='character')
locusnames<-scan('outflank_rap_pre.locusnames',what='character')
Rap_Dat<-MakeDiploidFSTMat(snpmat,locusnames,popnames)
head(Rap_Dat)

## Convert data matrix to data frame and trim the FST distribution to create your "null"
 # The "trim factor" parameters of this matrix allow you to set FST threshold cut offs
 # Number of samples is the number of putative populations/sampling locations
Out_Rap<-OutFLANK(FstDataFrame=Rap_Dat,LeftTrimFraction=0.05,RightTrimFraction=0.05,Hmin=0.05,NumberOfSamples=6,qthreshold=0.05)

## Subset by those that are outliers
Results<-Out_Rap[["results"]]
Flagged<-Results[,c("OutlierFlag")]
Out_New<-Results[which(Results$OutlierFlag=='TRUE'),]

## Produce a table of the outliers 
table(Results$GoodH)
sum(Flagged)
Out_Rap$numberLowFstOutliers
Out_Rap$numberHighFstOutliers
head(Out_Rap$Results)
sink()

## Write your output of all SNPs pairwise FST values to csv file
write.table(Rap_Dat,file="outflank_rap_pre.csv",row.names=FALSE,na="",col.names=TRUE,sep="\t",quote=FALSE)

#write output of statistical outliers to csv file
write.table(Out_Rap,file="outflank_rap_pre_summary.csv",row.names=FALSE,na="",col.names=TRUE,sep="\t",quote=FALSE)

## Create figure 
pdf("plot_rap_pre.pdf")
OutFLANKResultsPlotter(Out_New,withOutliers=TRUE,NoCorr=TRUE,Hmin=0.05,binwidth=0.005,RightZoomFraction=0.05,titletext=NULL)
dev.off()
