###########################################################################################################################
## 7a. Genetic Stock Identification
# Generate SNP Input files
###########################################################################################################################

#### Set your working directory
setwd("/Users/workingdirectory")

#### Load all the necessary R packages and install if necessary
library(tidyverse)
library(dplyr)
library(rrcovNA)
library(car)

#### Take the individual based assignments from DAPC and keep only individuals that had concordant
# DAPC group assignment and sampled population
# Let's start wiht K = 4
DAPC_Results<-read.csv("~/DAPC_Output")

#### Merge the sample ID values from the metadata with the genotype file in the same format "1_Sample ID"
pop_code<-read.csv("/METADATA_FILE")
pop_code<-pop_code[c("sample.id","Location","site.x","Time","GSI")]
test<-pop_code[pop_code$Time=="Pre",]
test$sample.id<-paste("1",test$sample.id,sep="_")

#### I used a couple of different methods to generate the SNPs and the loci for GSI analysis
# First, I'm going to remove all individuals that don't assign to the population they were sampled from
DAPC_Results$group<-as.factor(DAPC_Results$group)
levels(DAPC_Results$group)<-gsub("1","BD",levels(DAPC_Results$group))
levels(DAPC_Results$group)<-gsub("2","SBLR",levels(DAPC_Results$group))
levels(DAPC_Results$group)<-gsub("3","AD",levels(DAPC_Results$group))
levels(DAPC_Results$group)<-gsub("4","ID",levels(DAPC_Results$group))

#### Consolidate DAPC results into a dataframe for downstream analyses
DAPC_Results<-DAPC_Results[c("sample.id","PCA.pc.1","PCA.pc.2","PCA.pc.3","PCA.pc.4","group","LD1","LD2","LD3")]
dataframe<-merge(DAPC_Results,test,by="sample.id")

#### Let's look at the individuals whose DAPC group doesn't match their location then remove those locations
no_match<-filter(dataframe,group != Location)
match<-filter(dataframe,group == Location)
pop_code<-read.csv("/Users/alexandra.fraik/Desktop/pop_code_subset.csv")
colnames(pop_code)[1]<-"Indiv"
pop_code$sample.id<-paste("1",pop_code$Indiv,sep="_")
match<-match[c("sample.id")]
match_df<-merge(match,pop_code,by="sample.id")
match_df<-subset(match_df,select=-c(sample.id))
match_df$GSI_Time<-"Pre"

#### I'm also going to make a "no match" file which will include the pre individuals that didn't match their population assignment
# and the post-dam removal individuals that will be used for future analyses
post<-pop_code[pop_code$Time=="Post",]
no_match<-no_match[c("sample.id")]
pre_no_match_df<-merge(no_match,pop_code,by="sample.id")
test<-rbind(post,pre_no_match_df)
no_match_df<-subset(test,select=-c(sample.id))
no_match_df$GSI_Time<-"Post"

#### Make a dataframe for GSI
GSI_df<-rbind(match_df,no_match_df)

#### Also filter out this individual from SBLR: 33651_94
GSI_df<-GSI_df[!GSI_df$Indiv=="33651_94",]

#### Write the GSI dataframe output to a file
write.csv(GSI_df,"pop_code_GSI.csv")

## NOW, go and recaclulate pairwise FST with just the concordant indviduals
# Compare all Pre-dam removal "populations"
# Filter out Omy5 loci 
AD_BD<-read.csv("Pre_AD_BD_FST.weir.fst.csv")
AD_BD <- AD_BD %>% 
  filter(CHROM != "NC_035081.1") %>% 
  droplevels()

AD_ID<-read.csv("Pre_AD_ID_FST.weir.fst.csv")
AD_ID <- AD_ID %>% 
  filter(CHROM != "NC_035081.1") %>% 
  droplevels()

BD_ID<-read.csv("Pre_BD_ID_FST.weir.fst.csv")
BD_ID <- BD_ID %>% 
  filter(CHROM != "NC_035081.1") %>% 
  droplevels()

AD_SBLR<-read.csv("Pre_AD_SBLR_FST.weir.fst.csv")
AD_SBLR <- AD_SBLR %>% 
  filter(CHROM != "NC_035081.1") %>% 
  droplevels()

BD_SBLR<-read.csv("Pre_BD_SBLR_FST.weir.fst.csv")
BD_SBLR <- BD_SBLR %>% 
  filter(CHROM != "NC_035081.1") %>% 
  droplevels()

ID_SBLR<-read.csv("Pre_ID_SBLR_FST.weir.fst.csv")
ID_SBLR <- ID_SBLR %>% 
  filter(CHROM != "NC_035081.1") %>% 
  droplevels()

#### Retain the top 600 SNPs, filter out the duplicates
AD_ID_high<-AD_ID %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)
AD_BD_high<-AD_BD %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)
BD_ID_high<-BD_ID %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)
AD_SBLR_high<-AD_SBLR %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)
BD_SBLR_high<-BD_SBLR %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)
ID_SBLR_high<-ID_SBLR %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)

#### Make the SNP files for output
high<-rbind.data.frame(AD_ID_high,BD_ID_high,AD_BD_high,AD_SBLR_high,BD_SBLR_high,ID_SBLR_high)
high$SNPs<-paste0(high$CHROM,sep="_",high$POS)

#### Read in the Omy5 and Omy28 loci
Omy5_Omy28<-read.csv("~/Omy5_Omy28_SNPs.csv")

#### Read in the Omy5 and Omy28 Candidate SNPs to retain even after filtering out the Omy5 loci
high<-high[c(1:2,4)]
test<-rbind(high,Omy5_Omy28)
test <- test %>% 
  distinct(SNPs,.keep_all=TRUE)
high<-test[c(3)]

#### Write the output to the GSI file, < 600 SNPs
write.csv(high,"SNPs_Location_GSI.csv")

##############################################a##################################################################################
#### Just keep the SNPs from the AD/ID/BD Pairwise FST comparisons
################################################################################################################################

#### Other SNP comparisons I tried involved slicing more of the top SNPs
AD_ID_high<-AD_ID %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)
AD_BD_high<-AD_BD %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)
BD_ID_high<-BD_ID %>% arrange(desc(WEIR_AND_COCKERHAM_FST)) %>%  slice(1:100)

#### Make the SNP files for output
high<-rbind.data.frame(AD_ID_high,BD_ID_high,AD_BD_high)
high$SNPs<-paste0(high$CHROM,sep="_",high$POS)

#### Read in the Omy5 and Omy28 Candidate SNPs to retain even after filtering out the Omy5 loci
high<-high[c(1:2,4)]
test<-rbind(high,Omy5_Omy28)
test <- test %>% 
  distinct(SNPs,.keep_all=TRUE)
high<-test[c(3)]

#### Write the output to the GSI file, < 300 SNPs
write.csv(high,"SNPs_Location_GSI_1.csv")
