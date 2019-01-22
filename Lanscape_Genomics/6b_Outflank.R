###########################################################################################################################
## 5b. Outflank R Script
###########################################################################################################################

## First step, is to load the required package
library(OutFLANK)

## Read in the required snpmat, popnames and locus name files to make a datamatrix containing all of the data
snpmat<-read.csv('outflank_Variable.snpmat',sep='\t',as.is=TRUE,header=FALSE)
popnames<-scan('outflank_Variable.popnames',what='character')
locusnames<-scan('outflank_Variable.locusnames',what='character')
Variable_Data<-MakeDiploidFSTMat(snpmat,locusnames,popnames)
head(Variable_Data)

## Convert data matrix to data frame and trim the FST distribution to create your "null"
 # The "trim factor" parameters of this matrix allow you to set FST threshold cut offs
 # Number of samples is the number of putative populations/sampling locations
Variable_Out<-OutFLANK(FstDataFrame=Variable_Data,LeftTrimFraction=0.05,RightTrimFraction=0.05,Hmin=0.05,NumberOfSamples=6,qthreshold=0.05)

## Subset by those that are outliers
Results<-Variable_Out[["results"]]
Flagged<-Results[,c("OutlierFlag")]
Variable_Out_Subset<-Results[which(Results$OutlierFlag=='TRUE'),]

## Produce a table of the outliers 
table(Results$GoodH)
sum(Flagged)
Variable_Out$numberLowFstOutliers
Variable_Out$numberHighFstOutliers
head(Variable_Out$Results)
sink()

## Write your output of all SNPs pairwise FST values to csv file
write.table(Variable_Data,file="Outflank_Output/Variable_Data.csv",row.names=FALSE,na="",col.names=TRUE,sep="\t",quote=FALSE)

#write output of statistical outliers to csv file
write.table(Variable_Out,file="Outflank_Output/Variable_Data_Summary.csv",row.names=FALSE,na="",col.names=TRUE,sep="\t",quote=FALSE)

## Create figure 
pdf("plot_Variable.pdf")
OutFLANKResultsPlotter(Variable_Out_Subset,withOutliers=TRUE,NoCorr=TRUE,Hmin=0.05,binwidth=0.005,RightZoomFraction=0.05,titletext=NULL)
dev.off()
