###########################################################################################################################
## 6d. Population genetics analyses
# Create the input files for Adegenet
###########################################################################################################################

#### Set your working directory
setwd("/Users/workingdirectory")

#### Load all the necessary R packages and install if necessary
library(adegenet)
library(hierfstat)
library(pegas)

#### Name columns of the dataframe for the loci names
prefix<-"locus"
test<-rep(seq(1:71320))
loci.names<-paste(prefix,test,sep="_")

#### First calculate FIS for all pre-dam removal populations
# Load in the pre-dam removal adegenet file you created in scritp 6c_Population_Genetic_Statistics.sh 
elw_pre<-read.csv("Adegenet_Pre.csv",header=FALSE)

#### Now make the format correct for converting a data frame to a genind format
elw<-elw[c(-3,-4,-5,-6)]
elw<-elw[c(2,1,3:142642)]
colnames(elw)[1]<-"ID"
colnames(elw)[2]<-"Population"
pre_locus<-elw[,-c(1,2)]
pre_locus<-data.frame(mapply(paste0,pre_locus[c(T,F)],pre_locus[c(F,T)]))
colnames(pre_locus)[c(1:71320)]<-loci.names
pre_locus[pre_locus=="00"]<-NA

#### Now make the pre and post files
population_pre<-rep("Pre",times=567)
population_post<-rep("Post",times=558)
ind_pre<-as.character(elw$ID)
ind_post<-as.character(elw$ID)
population<-c(population_pre,population_post)

#### Merging it together in a proper format
pre_df<-df2genind(pre_locus,ind.names=ind_pre,pop=population_pre,ploidy=2,ncode=1)

#### Now calculate basic pop gen statistics for the dataframe
stats_elw<-basic.stats(pre_df,diploid=TRUE,digits=2)
write.table(stats_elw$Fis,"Pre_Elw_All_Pops.csv")

#### Now bootstrap the confidence intervals
FIS_Boot_CI_Pre<-boot.ppfis(pre_df,nboot=1000,quant=c(0.025,0.975))
FIS_Boot_CI_Pre<-as.matrix(FIS_Boot_CI_Pre)
write.table(FIS_Boot_CI_Pre,"Pre_Elw_All_Pops_FIS_CI.csv")
elw_fis<-as.data.frame(stats_elw$Fis)

#### Second calculate FIS for all post-dam removal life-history cohorts
# Load in the post-dam removal adegenet file you created in scritp 6c_Population_Genetic_Statistics.sh 
elw_post<-read.csv("Adegenet_Post.csv",header=FALSE)

#### Now make the format correct for converting a data frame to a genind format
elw<-elw[c(-3,-4,-5,-6)]
elw<-elw[c(2,1,3:142642)]
colnames(elw)[1]<-"ID"
colnames(elw)[2]<-"Cohort"
post_locus<-elw[,-c(1,2)]
post_locus<-data.frame(mapply(paste0,pre_locus[c(T,F)],pre_locus[c(F,T)]))
colnames(post_locus)[c(1:71320)]<-loci.names
post_locus[pre_locus=="00"]<-NA
post_df<-df2genind(post_locus,ind.names=ind_post,pop=population_post,ploidy=2,ncode=1)

#### Calculate basic pop gen statistics for the dataframe
stats_elw<-basic.stats(post_df,diploid=TRUE,digits=2)
write.table(stats_elw$Fis,"Post_Elw_All_Cohorts.csv")
FIS_Boot_CI_Post<-boot.ppfis(post_df,nboot=1000,quant=c(0.025,0.975))
FIS_Boot_CI_Post<-as.matrix(FIS_Boot_CI_Post)
write.table(FIS_Boot_CI_Post,"Post_Elw_All_Cohorts_FIS_CI.csv")

#### Finallly, calculate FIS for pre and post-dam removal data sets
# Now make the format correct for converting a data frame to a genind format
elw<-rbind(elw_pre,elw_post)
elw<-elw[c(-3,-4,-5,-6)]
elw<-elw[c(2,1,3:142642)]
colnames(elw)[1]<-"ID"
colnames(elw)[2]<-"Time"
locus<-elw[,-c(1,2)]
locus<-data.frame(mapply(paste0,locus[c(T,F)],locus[c(F,T)]))
colnames(locus)[c(1:71320)]<-loci.names
locus[locus=="00"]<-NA
df<-df2genind(locus,ind.names=ind,pop=population,ploidy=2,ncode=1)

#### Calculate basic pop gen statistics for the dataframe
stats_elw<-basic.stats(df,diploid=TRUE,digits=2)
write.table(stats_elw$Fis,"Elw_Pre_Post.csv")

#### Now bootstrap confidence intervals
FIS_Boot_CI<-boot.ppfis(df,nboot=1000,quant=c(0.025,0.975))
FIS_Boot_CI<-as.matrix(FIS_Boot_CI)
write.table(FIS_Boot_CI,"Elw_Pre_Post_FIS_CI.csv")
