###########################################################################################################################
## 5. Population structure analyses (part 2)
# Discriminant Analysis of Principal Components to perform an indivdual based clustering population structure analysis
###########################################################################################################################

#### Set your working directory
setwd("/Users/workingdirectory")

#### Load all the necessary R packages and install if necessary
library(adegenet)
library(vcfR)
library(ggplot2)
library(dplyr)
library(grid)
library(stringr)
library(gridExtra)

#### Read in the files that were filtered in VCFtools
# At the breaks of my script "2_VCF_to_PLINK.sh" I also used this same script to examine different filtering
# parameters to determine which I wanted to use for downstream analyses

#### If this looks similar to the PCA script, it should, don't panic
# Similar set of tools and analysis really except we're discretizing variance in genetic data among individuals
# to try to explain group membership/population structure

#### Similar to the 4b_FastSTRUCTURE.R script, you can definitely loop through or use "apply" to just produce DAPC plots and results
# in ggplot in R for all your iterations, but I am going to just show an example of how I made the pre-dam removal plot for K = 4

#### In part, I am doing this because you have to still take the steps to pick out how many clusters you want to examine
# Read in your VCF file - in my case, pre-dam removal
pre<-read.vcfR("~/Elw_Pre.vcf")

#### SHOULD NOT be necessary at this point, but if some non-biallelic alleles got through the filtering cracks, this will catch them
pre<-pre[is.biallelic(pre)]

#### Convert it into a genind object
genind<-vcfR2genind(pre)

#### Now, find the number of clusters
# max.n.clust = how many clusters (at most) you want the program to consider
# how many pcas TOPS (CANNOT BE MORE THAN THE NUMBER OF INDVIDUALS IN YOUR DATA SET)
genind_clusters<-find.clusters(genind,max.n.clust=20,n.pca=500,n.clust=4)

#### Finally, after playing and exploring, you will either interactively or hard code the number of clusters you want 
# DAPC to use (n.clust = n.da)
dapc<-dapc(genind,genind_clusters$grp,n.pca=500,n.da=4)

#### Extract necessary information from DAPC for plotting
# pca loading scores
pca.here<-as.data.frame(dapc$pca.loadings)
pca.here<-pca.here[1:5]
# Including indivdual loading scores
indiv_loadings<-as.data.frame(dapc$loadings)
tab<-dapc$tab
tab_df<-tab[1:4] 
# group membership
group<-dapc$grp
# individual coordintes
coords<-dapc$ind.coord
# posterior probability of membership
posterior<-dapc$posterior
dapc_df<-data.frame(tab_df,group,coords,posterior)
dapc_df$sample.id<-rownames(dapc_df)

#### Smoosh all this info into a file, just in case, or in case you want to rerun this script in the future but don't want to rerun
# the first steps of this file, you can just start from here and read this file back in
write.csv(dapc_df,"dapc_df_all.csv")

#### Load in that dang necessary meatadata file aganin
pop_code<-read.csv("~/METADATA_File",header=TRUE)

#### Merge the sample ID values from the metadata with the genotype file in the same format "1_Sample ID"
pop_code$sample.id<-paste("1",pop_code$sample.id,sep="_")
colnames(pop_code)[1]<-"sample.id"
dapc_df<-merge(dapc_df,pop_code,by="sample.id")

#### Recode all DAPC population "group" points as factors and assign them a "point value" for plotting within the metadata file
# DAPC group assignment = point shape
# Life-history form = point size

points_group<-c(10,8,22,9,12,11,19)
dapc_df$group<-as.factor(dapc_df$group)
names(points_group)<-levels(dapc_df$group)
points<-data.frame(group=names(points_group),points_group=points_group)
dapc_df<-merge(dapc_df,points,by="group")

#### Code all populations with colors and assign them "a "point value"color" for plotting within the metadata file
# Population = color
colors_location<-c("brown1","blue2","darkgreen","limegreen")
names(colors_location)<-levels(dapc_df$Location)
colors_l<-data.frame(Location=names(colors_location),colors_location=colors_location)
dapc_df<-merge(dapc_df,colors_l,by="Location")

#### Reorder file to have individuals organized chronologically: for pre THEN post-dam removal
# for facetting plots in ggplot
neworder<-c("Pre","Post")
dapc_df<-arrange(transform(dapc_df,Time=factor(Time,levels=neworder)),Time)

#### Reorder file to sort sampling locations (referred to as collections in GSI analyses) from upstream to downstream 
# In this study, segements of the river were split into populations based on relative anadromous barrier location 
# Then, within these populations, we had distinctive sampling locations that we want to be listed in the key
# from upriver -> downriver
newerorder<-c("Chicago_Camp","Near_Delabarre_Creek","Wilder","Headwaters","Hayes","Lost_River","Elkhorn","Long_Creek","Geyser","Whiskey_Bend",
              "Cat_Creek","Boulder","Glines","Altaire","Griff_Creek","Hughes_Creek","Campground_Creek","Madison_Creek",
              "South_Branch_Little_River","Little_River","Indian_Creek","Aldwell","Elwha_River")
df<-arrange(transform(dapc_df,site.x=factor(site.x,levels=newerorder),site.x))
pop_code_colors<-read.csv("/Users/alexandra.fraik/Desktop/pop_code_colors.csv")

#### Remove underscores for collections to make the legend more readable
df<-merge(df,pop_code_colors,by="site.x",all.x=TRUE)
df$site.x<-gsub('_',' ',df$site.x)
color<-as.character(df$color)
names(color)<-as.character(df$site.x)

#### Rename the levels for other factors to make the legend more readable
levels(df$Time)<-c("Post Dam Removal","Prior to Dam Removal")
levels(df$run.timing)<-c("Summer","Unknown","Winter")

######################################################################################
############ DAPC Population Groups
######################################################################################
#### Make a plot for group memebership (sampling site VS group membership)
# Pre-dam removal, K = 4
pdf(file="DAPC_Group_Figure_K4_Pre.pdf",width=14,height=12)
plotc<-ggplot(df,aes(x=site.x,y=group),group=group)
plotc<-plotc+geom_point(aes(fill=Location,color=Location,shape=group),size=8)
plotc<-plotc+labs(x="Population",y="DAPC Group")
plotc<-plotc+theme_bw()
plotc<-plotc+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
plotc<-plotc+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plotc<-plotc+theme(axis.text.x=element_text(color="black",size=25,angle=90,hjust=1),axis.title.x=element_text(size=30,face="bold"))
plotc<-plotc+scale_color_manual(name="Population",values=colors_location)
plotc<-plotc+scale_fill_manual(name="Population",values=colors_location)
plotc<-plotc+scale_shape_manual(name="DAPC\nGroup",values=c(points_group))
plotc<-plotc+facet_grid(switch="x",scales="free_x",space="free_x")
plotc<-plotc+theme(strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plotc<-plotc+theme(legend.text=element_text(size=rel(2)))
plotc<-plotc+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plotc<-plotc+guides(color=guide_legend(ncol=2,byrow=FALSE,override.aes=list(size=4)),shape=guide_legend(ncol=2,override.aes=list(size=4)))
plotc<-plotc+theme(legend.position="bottom")
plotc
dev.off()

######################################################################################
############ DAPC Population Groups
######################################################################################
#### Make a plot of DAPC Axis 1 VS DAPC Axis 2
# Pre-dam removal, K = 4
pdf(file="DAPC_K4_Pre.pdf",width=22,height=12)
plotc<-plotc<-ggplot(df,aes(x=LD1,y=LD2),group=group)
plotc<-plotc+geom_point(aes(fill=Location,color=Location,shape=group,size=Species.1))
plotc<-plotc+labs(x="DAPC Axis 1",y="DAPC Axis 2")
plotc<-plotc+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=30,face="bold"))
plotc<-plotc+theme(axis.text.x=element_text(color="black",size=25),axis.title.x=element_text(size=30,face="bold"))
plotc<-plotc+scale_color_manual(name="Population",values=colors_location)
plotc<-plotc+scale_fill_manual(name="Population",values=colors_location)
plotc<-plotc+scale_shape_manual(name="DAPC\nGroup",values=c(points_group))
plotc<-plotc+scale_size_manual(name="Phenotypic\nForm",values=c(8,12,8))
plotc<-plotc+facet_grid(switch="x",scales="free_x",space="free_x")
plotc<-plotc+theme(strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plotc<-plotc+theme(legend.text=element_text(size=rel(2.75)))
plotc<-plotc+theme(legend.title=element_text(face="bold",size=rel(3)))
plotc<-plotc+guides(color=guide_legend(ncol=2,override.aes=list(size=5)),shape=guide_legend(ncol=2,override.aes=list(size=5)),size=guide_legend(ncol=1,byrow=FALSE))
plotc<-plotc+theme(legend.position="bottom")
plotc
dev.off()

######################################################################################
############ DAPC Population Groups
######################################################################################
#### We first define the chromosome order
chrorder<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","Mitochondria","Unknown")

#### Load in a file with every SNP and every chromosome that it is on for annotation
snp_position<-read.csv("~/Chromosome_Known_MAP.csv")

#### PC1
# Take the top 100 SNPs that explain the highest proportion of variance per PC
pca.c1<-pca.here^2
pca.c1$SNP<-rownames(pca.c1)
pca.c1$SNP<-gsub("\\..*","",pca.c1$SNP)
colnames(pca.c1)<-gsub("\\-*","",colnames(pca.c1))
pca.c1<-merge(pca.c1,snp_position,by="SNP")
pca.c1$Chromosome<-as.factor(pca.c1$Chromosome)

#### PC1
PC1<-pca.c1%>%arrange(desc(PCApa.1))%>% slice(1:713)
PC1$Chromosome[PC1$Chromosome == "30"] <- "Unknown"
PC1$PC<-"PC1"
PC1<-arrange(transform(PC1,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC1<-data.frame(SNP=PC1$SNP,Loading_Score=PC1$PCApa.1,Chromosome=PC1$Chromosome,PC=PC1$PC)

#### PC2
PC2<-pca.c1%>%arrange(desc(PCApa.2))%>% slice(1:713)
PC2$Chromosome[PC2$Chromosome == "30"] <- "Unknown"
PC2$PC<-"PC2"
PC2<-arrange(transform(PC2,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC2<-data.frame(SNP=PC2$SNP,Loading_Score=PC2$PCApa.2,Chromosome=PC2$Chromosome,PC=PC2$PC)

#### PC3
PC3<-pca.c1%>%arrange(desc(PCApa.3))%>% slice(1:713)
PC3$Chromosome[PC3$Chromosome == "30"] <- "Unknown"
PC3$PC<-"PC3"
PC3<-arrange(transform(PC3,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC3<-data.frame(SNP=PC3$SNP,Loading_Score=PC3$PCApa.3,Chromosome=PC3$Chromosome,PC=PC3$PC)

#### PC4
PC4<-pca.c1%>%arrange(desc(PCApa.4))%>% slice(1:713)
PC4$Chromosome[PC4$Chromosome == "40"] <- "Unknown"
PC4$PC<-"PC4"
PC4<-arrange(transform(PC4,Chromosome=factor(Chromosome,levels=chrorder)),Chromosome)
PC4<-data.frame(SNP=PC4$SNP,Loading_Score=PC4$PCApa.4,Chromosome=PC4$Chromosome,PC=PC4$PC)

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
