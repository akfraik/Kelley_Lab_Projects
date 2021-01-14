###########################################################################################################################
## 3. Principal Components Analysis
# To taken an unsupervised approach to examining clustering patterns and investigate different filtering parameters
###########################################################################################################################

#### Set your working directory
setwd("/Users/workingdirectory")

#### Load all the necessary R packages and install if necessary
library(vcfR)
library(ggplot2)
library(dplyr)
library(stringr)
library(adegenet)
library(adespatial)
library(reconproGS)
library(gridExtra)
library(grid)

#### Read in the files that were filtered in VCFtools
# At the breaks of my script "2_VCF_to_PLINK.sh" I also used this same script to examine different filtering
# parameters to determine which I wanted to use for downstream analyses

#### However, in this script I am going to focus on after you have already done these exploratory PCAs
# and you are playing with the PCA to look at clustering for your final data set
# So, we are going to run this locally and download the VCF file to our local working directory

#### Load in VCF file from local directory and convert it into the genind format from Adegenet
vcf<-read.vcfR("~/all_81_vcf.vcf")
test<-vcfR2genind(vcf)
testing<-test$tab
f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)
}
test<-apply(testing,2,f1)

### You can use this to subset for chromosomes of interest- in this case, I was interested in looking at clustering on 
# chromosomes 5 and 28 respectively (Omy5 and Omy28)
#Omy5<-select(test,contains("NC_035081_1"))
#Omy28<-select(test,contains("NC_035104_1"))

#### Convert whatever final data set you are using into a dataframe for PCA
test<-as.data.frame(test)

#### Now actually conducte the PCA and extract the information you need
# including the coordinates, eigen values, etc.
pca.here<-dudi.pca(test,center=TRUE,scale=FALSE,scannf=FALSE,nf=4)
table<-pca.here$tab
li<-pca.here$li
l1<-pca.here$l1
co<-pca.here$co
c1<-pca.here$c1

#### Read in the metadata information 
pop_code<-read.csv("~/METADATA_File",header=TRUE)

#### Merge the sample ID values from the metadata with the genotype file in the same format "1_Sample ID"
# this is what VCF does to sample ids
pop_code$sample.id<-paste("1",pop_code$sample.id,sep="_")

#### I was interested in looking at clustering by
# Life-history cohort = point shape
# Run-timing = point size
# Population = color

#### Recode life-history cohort values (GSI) as factors and assign them a "point value" for plotting within the metadata file
pop_code$GSI<-as.factor(pop_code$GSI)
points_gsi<-c(10,8,22,9,12,11,19)
names(points_gsi)<-levels(pop_code$GSI)
points<-data.frame(GSI=names(points_gsi),points=points_gsi)
pop_code<-merge(pop_code,points,by="GSI")

#### Code all populations with colors and assign them "a "point value"color" for plotting within the metadata file
colors_location<-c("brown1","blue2","darkgreen","limegreen")
names(colors_location)<-levels(pop_code$Location)
colors_l<-data.frame(Location=names(colors_location),colors_location=colors_location)
pop_code<-merge(pop_code,colors_l,by="Location")
colors_location<-as.character(pop_code$colors_location)
names(colors_location)<-as.character(pop_code$Location)
table$sample.id<-rownames(table)
df<-merge(table,pop_code,by="sample.id")

#### Reorder file to have individuals organized chronologically: for pre THEN post-dam removal
# for facetting plots in ggplot
neworder<-c("Pre","Post")
test<-arrange(transform(df,Time=factor(Time,levels=neworder)),Time)

#### Reorder file to sort sampling locations (referred to as collections in GSI analyses) from upstream to downstream 
# In this study, segements of the river were split into populations based on relative anadromous barrier location 
# Then, within these populations, we had distinctive sampling locations that we want to be listed in the key
# from upriver -> downriver
newerorder<-c("Chicago_Camp","Near_Delabarre_Creek","Wilder","Headwaters","Hayes","Lost_River","Elkhorn","Long_Creek","Geyser","Whiskey_Bend",
              "Cat_Creek","Boulder","Glines","Altaire","Griff_Creek","Hughes_Creek","Campground_Creek","Madison_Creek",
              "South_Branch_Little_River","Little_River","Indian_Creek","Aldwell","Elwha_River")
pca_df<-arrange(transform(test,site.x=factor(site.x,levels=newerorder),site.x))
pop_code_colors<-read.csv("/Users/alexandra.fraik/Desktop/pop_code_colors.csv")
pca_df<-merge(pca_df,pop_code_colors,by="site.x",all.x=TRUE)

#### Remove underscores for collections to make the legend more readable
pca_df$site.x<-gsub('_',' ',pca_df$site.x)
color<-as.character(pca_df$color)
names(color)<-as.character(pca_df$site.x)

## Rename the levels for other factors to make the legend more readable
levels(pca_df$Time)<-c("Prior to Dam Removal","Post Dam Removal")
levels(pca_df$run.timing)<-c("Summer","Unknown","Winter")

################################################################################################################################################      
############################### Generate PCA for Adgenet output ##############################################################
################################################################################################################################################
## Make a dataframe of your PCA output to only contain the values you are interested in plotting
tab<-data.frame(sample.id=pca_df$sample.id,accession.x=pca_df$accession.x,site.x=pca_df$site.x,year=pca_df$year,Hatchery=pca_df$Nvh,
                EV1=li$Axis1,EV2=li$Axis2,EV3=li$Axis3,EV4=li$Axis4,color=pca_df$color,colors_location=pca_df$colors_location,run.timing=pca_df$run.timing,
                Time=pca_df$Time,Species=pca_df$Species.1,Location=pca_df$Location,GSI=pca_df$GSI,points_gsi=pca_df$points,
                stringsAsFactors=FALSE)
head(tab)
write.csv(tab,"PCA_EV_Adegenet.csv")
write.csv(c1,"PC_Loading_Scores_Adegenet.csv")

################################################################################################################################################      
############################### Generate PCAs by Sampling location ##############################################################
################################################################################################################################################
## colored location by population, points by life-history cohort, size of point by life-history form
pdf(file="Pop_ggplot_RAD_Location.pdf",width=24,height=20)
par(oma=c(1,2,1,10))
plot<-ggplot(tab,aes(x=EV1,y=EV2),group=Location)
plot<-plot+geom_point(aes(fill=Location,shape=GSI,color=Location,size=run.timing))
plot<-plot+labs(x="PC1",y="PC2")
plot<-plot+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plot<-plot+scale_color_manual(name="Sampling\nSites",values=colors_location)
plot<-plot+scale_fill_manual(name="Sampling\nSites",values=colors_location)
plot<-plot+scale_size_manual(name="Run-\nTiming",values=c(12,8,8))
plot<-plot+scale_shape_manual(name="Life-history\nCohort",values=c(points_gsi))
plot<-plot+facet_grid(~Time,switch="x")
plot<-plot+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="none")
plota<-ggplot(tab,aes(x=EV3,y=EV4),group=Location)
plota<-plota+geom_point(aes(fill=Location,shape=GSI,color=Location,size=run.timing))
plota<-plota+labs(x="PC3",y="PC4")
plota<-plota+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plota<-plota+scale_color_manual(name="Sampling\nSites",values=colors_location)
plota<-plota+scale_fill_manual(name="Sampling\nSites",values=colors_location)
plota<-plota+scale_size_manual(name="Run-\nTiming",values=c(12,8,8))
plota<-plota+scale_shape_manual(name="Life-history\nCohort",values=c(points_gsi),labels=c("Adult 2015","Adult 2016","Adult 2017","Adult Pre 2015","Juvenile 2016","Juvenile 2017","Pre Dam Removal"))
plota<-plota+facet_grid(~Time,switch="x")
plota<-plota+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(2.75)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(3)))
plota<-plota+theme(legend.position="bottom")
plota<-plota+guides(color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=2,byrow=FALSE,override.aes=list(size=5)),size=guide_legend(ncol=1,byrow=FALSE))
grid.arrange(plot,plota,nrow=2,ncol=1)
dev.off()

################################################################################################################################################      
############################### Generate PCAs for technical artifacts ##############################################################
################################################################################################################################################
## 1) So, now we would like to determine whether or not accession ID (DNA extraction batch, sampling ID, and frequently how individuals were fitlered)
## Let's try to generate colors with orange as the low, palegreen as mid, midnight blue as high
tab$accession.x<-as.factor(tab$accession.x)
color_pallete_function<-colorRampPalette(colors=c("orange","palegreen","midnightblue"),space="Lab")
num_colors<-nlevels(tab$accession.x)
accession_colors<-color_pallete_function(num_colors)
names(accession_colors)<-levels(tab$accession.x)

## Plot
pdf(file="PCA_Accession_RAD.pdf",width=24,height=20)
par(oma=c(1,2,1,10))
plot<-ggplot(tab,aes(x=EV1,y=EV2),group=accession.x,shape=19)
plot<-plot+geom_point(aes(color=accession.x),size=6,alpha=0.7)
plot<-plot+labs(x="PC1",y="PC2")
plot<-plot+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plot<-plot+scale_color_manual(name="Accession ID",values=accession_colors[tab$accession.x])
plot<-plot+facet_grid(~Time,switch="x")
plot<-plot+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="none")
plota<-ggplot(tab,aes(x=EV3,y=EV4),group=accession.x,shape=19)
plota<-plota+geom_point(aes(color=accession.x),size=6,alpha=0.7)
plota<-plota+labs(x="PC3",y="PC4")
plota<-plota+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plota<-plota+scale_color_manual(name="Accession ID",values=accession_colors[tab$accession.x])
plota<-plota+facet_grid(~Time,switch="x")
plota<-plota+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(2)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plota<-plota+theme(legend.position="bottom")
plota<-plota+guides(color=guide_legend(ncol=13,byrow=FALSE,override.aes=list(size=5)))
grid.arrange(plot,plota,nrow=2,ncol=1)
dev.off()

## 2) We also would like to determine whether or not read depth impacted clustering patterns
# To do this, let's read in a file with read depth information from Krista
seq_data<-read.csv("/Users/alexandra.fraik/Documents/Info/Alex_Elw_paper1_INDV_seqstats.csv")
seq_data$sample.id<-paste("1",seq_data$sample.id,sep="_")
seq_data<-seq_data[c("sample.id","readsMapped","reads_properpairs")]
tab_merged<-merge(tab,seq_data,by="sample.id")

## Now that we have a merged file, let's try to bin relative mapping rates into ~40 bins and let's start with reads proper pairs
tab<-tab_merged %>%
  mutate(MyBins=cut(reads_properpairs,
                    breaks=unique(quantile(reads_properpairs,probs=seq.int(0,1,by=1/11))),include.lowest=TRUE))
tab$MyBins<-as.factor(tab$MyBins)
color_pallete_function<-colorRampPalette(colors=c("orange","palegreen","midnightblue"),space="Lab")
num_colors<-nlevels(tab$MyBins)
bins_colors<-color_pallete_function(num_colors)
names(bins_colors)<-levels(tab$MyBins)

## Now plot it
pdf(file="PCA_Reads_Proper_Pairs_RAD.pdf",width=24,height=20)
par(oma=c(1,2,1,10))
plot<-ggplot(tab,aes(x=EV1,y=EV2),group=MyBins,shape=19)
plot<-plot+geom_point(aes(color=MyBins),size=6,alpha=0.7)
plot<-plot+labs(x="PC1",y="PC2")
plot<-plot+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plot<-plot+scale_color_manual(name="Number Reads\nProperly Paired",values=bins_colors[tab$MyBins])
plot<-plot+facet_grid(~Time,switch="x")
plot<-plot+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="none")
plota<-ggplot(tab,aes(x=EV3,y=EV4),group=MyBins,shape=19)
plota<-plota+geom_point(aes(color=MyBins),size=6,alpha=0.7)
plota<-plota+labs(x="PC3",y="PC4")
plota<-plota+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plota<-plota+scale_color_manual(name="Number Reads\nProperly Paired",values=bins_colors[tab$MyBins])
plota<-plota+facet_grid(~Time,switch="x")
plota<-plota+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(2)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plota<-plota+theme(legend.position="bottom")
plota<-plota+guides(color=guide_legend(ncol=5,byrow=FALSE,override.aes=list(size=5)))
grid.arrange(plot,plota,nrow=2,ncol=1)
dev.off()

## 3) Finally,we will determine whether or not proportion of missing data impacted clustering patterns
# Regarding counts of missing data: Each count represents the number of individuals with missing data at each locus. 
# The last column, "mean" can be thought of as the average number of individuals with missing data per locus.
# Regarding percentage missing data: This percentage is relative to the population and locus, not to the entire data set. 
# The last column, "mean" represents the average percent of the population with missing data per locus.
missing<-read.csv("/Users/alexandra.fraik/Documents/Info/all_81_missingness.csv")
colnames(missing)[2]<-"sample.id"
missing<-missing[c("sample.id","N_MISS")]
tab_merged<-merge(tab,missing,by="sample.id")

## Now that we have a merged file, let's try to bin relative mapping rates into ~40 bins and let's start with reads proper pairs
tab<-tab_merged %>%
  mutate(MissBins=cut(N_MISS,
                    breaks=unique(quantile(N_MISS,probs=seq.int(0,1,by=1/20))),include.lowest=TRUE))
tab$MissBins<-as.factor(tab$MissBins)
color_pallete_function<-colorRampPalette(colors=c("orange","palegreen","midnightblue"),space="Lab")
num_colors<-nlevels(tab$MissBins)
bins_colors<-color_pallete_function(num_colors)
names(bins_colors)<-levels(tab$MissBins)

## Now plot it
pdf(file="PCA_MISS_RAD.pdf",width=24,height=20)
par(oma=c(1,2,1,10))
plot<-ggplot(tab,aes(x=EV1,y=EV2),group=MissBins,shape=19)
plot<-plot+geom_point(aes(color=MissBins),size=6,alpha=0.7)
plot<-plot+labs(x="PC1",y="PC2")
plot<-plot+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plot<-plot+scale_color_manual(name="Missing Data\nPer Individual",values=bins_colors[tab$MissBins])
plot<-plot+facet_grid(~Time,switch="x")
plot<-plot+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="none")
plota<-ggplot(tab,aes(x=EV3,y=EV4),group=MissBins,shape=19)
plota<-plota+geom_point(aes(color=MissBins),size=6,alpha=0.7)
plota<-plota+labs(x="PC3",y="PC4")
plota<-plota+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plota<-plota+scale_color_manual(name="Missing Data\nPer Individual",values=bins_colors[tab$MissBins])
plota<-plota+facet_grid(~Time,switch="x")
plota<-plota+theme(strip.text=element_text(face="bold",size=30),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(2)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plota<-plota+theme(legend.position="bottom")
plota<-plota+guides(color=guide_legend(ncol=5,byrow=FALSE,override.aes=list(size=5)))
grid.arrange(plot,plota,nrow=2,ncol=1)
dev.off()

###########################################################################################################################################################################
########## Now let's organize the PCA allele plot by loadings
###########################################################################################################################################################################
#### We first define the chromosome order
chrorder<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","Mitochondria","Unknown")

#### Load in a file with every SNP and every chromosome that it is on for annotation
snp_position<-read.csv("~/Chromosome_Known_MAP.csv")

#### PC1
# Take the top 100 SNPs that explain the highest proportion of variance per PC
pca.c1<-pca.here$c1^2
pca.c1$SNP<-rownames(pca.c1)
pca.c1$SNP<-gsub("\\..*","",pca.c1$SNP)
pca.c1<-merge(pca.c1,snp_position,by="SNP")
pca.c1$Chromosome<-as.factor(pca.c1$Chromosome)
PC1<-pca.c1 %>% arrange(desc(CS1))%>% slice(1:713)
PC1$Chromosome[PC1$Chromosome == "30"] <- "Unknown"
PC1$PC<-"PC1"
PC1<-arrange(transform(PC1,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC1<-data.frame(SNP=PC1$SNP,Loading_Score=PC1$CS1,Chromosome=PC1$Chromosome,PC=PC1$PC)

#### PC2
# Take the top 100 SNPs that explain the highest proportion of variance per PC
PC2<-pca.c1 %>% arrange(desc(CS2))%>% slice(1:713)
PC2$Chromosome[PC2$Chromosome == "30"] <- "Unknown"
PC2$PC<-"PC2"
PC2<-arrange(transform(PC2,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC2<-data.frame(SNP=PC2$SNP,Loading_Score=PC2$CS2,Chromosome=PC2$Chromosome,PC=PC2$PC)

#### PC3
# Take the top 100 SNPs that explain the highest proportion of variance per PC
PC3<-pca.c1%>%arrange(desc(CS3))%>% slice(1:713)
PC3$Chromosome[PC3$Chromosome == "30"] <- "Unknown"
PC3$PC<-"PC3"
PC3<-arrange(transform(PC3,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC3<-data.frame(SNP=PC3$SNP,Loading_Score=PC3$CS3,Chromosome=PC3$Chromosome,PC=PC3$PC)

#### PC4
# Take the top 100 SNPs that explain the highest proportion of variance per PC
PC4<-pca.c1%>%arrange(desc(CS4))%>% slice(1:713)
PC4$Chromosome[PC4$Chromosome == "30"] <- "Unknown"
PC4$PC<-"PC4"
PC4<-arrange(transform(PC4,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC4<-data.frame(SNP=PC4$SNP,Loading_Score=PC4$CS4,Chromosome=PC4$Chromosome,PC=PC4$PC)

#### Make PCA dataframe and write an output file of the top 1% of SNPs loading onto each PC
PCA_df<-rbind.data.frame(PC1,PC2,PC3,PC4)
write.csv(PCA_df,"PCA_Loadings.csv")

#### Make a plot of each of the PCs top 1% of SNPs loading on the plot
pdf("PCA_Stacked_Loading_Scores.pdf",width=12,height=6)
plot<-ggplot(data=PCA_df,aes(x=Chromosome,y=Loading_Score,fill=Chromosome))
plot<-plot+geom_bar(stat="identity", color="black",position=position_dodge())
plot<-plot+labs(x="Chromsome Number",y="Loadings")
plot<-plot+theme(axis.text.y=element_text(color="black",size=35),axis.title.y=element_text(size=45,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=35),axis.title.x=element_text(size=45,face="bold"))
plot<-plot+theme_minimal()
plot<-plot+facet_grid(PC~.,switch="x")
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="none")
plot<-plot+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
plot
dev.off()
