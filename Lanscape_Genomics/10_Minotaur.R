#!/usr/bin/env Rscript
#R

## Install the required packages and read in the librarires
 # install.packages('MINOTAUR')
 # require(devtools)
 #' install_github('NESCent/MINOTAUR',build_vignettes=TRUE)
  library(MINOTAUR)
  library(rrcovNA)
  observations_pre<-read.csv("/Minotaur/Minotaur_Pre.csv",header=TRUE)
  observations_post<-read.csv("/Minotaur/Minotaur_Pre.csv",header=TRUE)  
  
#################################################Pre-DFTD########################################################################
## Mahalanobis distance is a multi-dimensional measure of the number of standard deviations that a point  
 # lies from the mean of the distribution
  Md<-Mahalanobis(observations,column.nums=c(2:4,6:12,14:20))
  pdf(file="/Minotaur/Minotaur_Figure_All_Pre_DFTD_Md")
  hist(Md,bins=100)
  dev.off()
  Md<-cbind(observations[,1],Md)
  write.csv(Md,"/Minotaur/Minotaur_Md_Pre_DFTD.csv", quote=FALSE)

## Now look at ordination among just the outlier tests: XtX, Outflank and PCAdapt
  Md_FST_Pre_DFTD<-Mahalanobis(observations_pre,column.nums=c(18:20))
  pdf(file="/Minotaur/Minotaur_Fst_Pre_Md")
  hist(Md_FST_Pre_DFTD,bins=100)
  dev.off()
  Md_Pre_DFTD<-cbind(observations_pre[,1],Md_FST_Pre_DFTD)
  write.csv(Md_FST_Pre_DFTD,"/Minotaur/Minotaur_Fst_Md_Pre.csv",quote=FALSE)

## Just GEAs- Bayenv2 and LFMM output
  Md_GEA_Pre_DFTD<-Mahalanobis(observations_pre,column.nums=c(2:17))
  pdf(file="/Minotaur/Md_GEA_Pre_DFTD")
  hist(Md_GEA_Pre_DFTD,bins=100)
  Md_GEA_Pre_DFTD<-cbind(observations_pre[,1],Md_GEA_Pre_DFTD)
  dev.off()
  write.csv(Md_GEA_Pre_DFTD,"Minotaur/Md_GEA_Pre_DFTD.csv",quote=FALSE)
  
#################################################Post-DFTD########################################################################  
## Mahalanobis distance is a multi-dimensional measure of the number of standard deviations that a point  
 # lies from the mean of the distribution
  Md<-Mahalanobis(observations_post,column.nums=c(2:4,6:12,14:20))
  pdf(file="/Minotaur/Minotaur_Figure_All_Post_DFTD_Md")
  hist(Md,bins=100)
  dev.off()
  Md<-cbind(observations_post[,1],Md)
  write.csv(Md,"/Minotaur/Minotaur_Md_Post_DFTD.csv", quote=FALSE)

## Now look at ordination among just the outlier tests: XtX, Outflank and PCAdapt
  Md_FST_Pre_DFTD<-Mahalanobis(observations_post,column.nums=c(18:20))
  pdf(file="/Minotaur/Minotaur_Fst_Post_Md")
  hist(Md_FST_Post_DFTD,bins=100)
  dev.off()
  Md_Post_DFTD<-cbind(observations_post[,1],Md_FST_Post_DFTD)
  write.csv(Md_FST_Post_DFTD,"/Minotaur/Minotaur_Fst_Md_Post.csv",quote=FALSE)

## Just GEAs- Bayenv2 and LFMM output
  Md_GEA_Pre_DFTD<-Mahalanobis(observations_post,column.nums=c(2:17))
  pdf(file="/Minotaur/Md_GEA_Post_DFTD")
  hist(Md_GEA_Post_DFTD,bins=100)
  Md_GEA_Post_DFTD<-cbind(observations_post[,1],Md_GEA_Pre_DFTD)
  dev.off()
  write.csv(Md_GEA_Post_DFTD,"Minotaur/Md_GEA_Post_DFTD.csv",quote=FALSE)
