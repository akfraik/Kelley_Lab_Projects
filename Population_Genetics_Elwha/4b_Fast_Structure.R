###########################################################################################################################
## 4b. Population structure analyses
# Produce fast Structure figures
###########################################################################################################################
#### I really hate the distruct plots from FastSTRUCTURE, so I'm going to make prettier ones in ggplot

#### Set your working directory
setwd("/Users/workingdirectory")

#### Load all the necessary R packages and install if necessary
library(ggplot2)
library(dplyr)
library(vcfR)
library(stringr)
library(gridExtra)
library(tidyr)
library(ggdendro)
library(reshape2)

#### Load necessary meatadata files
pop_code<-read.csv("~/METADATA_File",header=TRUE)

#### Subset the metadata file for what you need and color the populations
colors_location<-c("brown1","darkgreen","limegreen","blue2")
pop_code<-pop_code[c("sample.id","Location","site.x","Time","Species.1")]

#### You can definitely loop through or use "apply" to just produce FastSTRUCTURE plots in ggplot in R
# for all your iterations, but I am going to just show an example of how I made the pre-dam removal plot for K = 4

#### 
# Start with K = 4 which is the number of model components that FastSTRUCTURE's "ChooseK.py" script selected
samples<-read.csv("Pre_bed_nosex.csv",header=FALSE)
colnames(samples)[1]<-"Nothing"
colnames(samples)[2]<-"sample.id"
tbl_4<-read.table("Pre.4.meanQ")
here<-cbind(samples,tbl_4)

#### Merge your metadata file with the output from FastSTRUCTURE
Pre<-merge(here,pop_code,by="sample.id")

#### Gather your variables so that you can use like barplots to plot the quantiles of different ancestry per individual
here<-Pre %>% gather(Variable,Ancestry,-site.x,-sample.id,-Nothing,-Location,-Time,-Species.1)
here$Variable<-as.factor(here$Variable)
here$Genotype<-as.factor(here$Ancestry)
here <- here %>% drop_na(Variable)

#### I was interested in looking at clustering by population abd sampling location
# population = color
# order = upriver sampling locations -> down river samplig locations
names(colors_location)<-levels(as.factor(here$Variable))
colors_l<-data.frame(Variable=names(colors_location),colors_location=colors_location)
df<-merge(here,colors_l,by="Variable")

#### Defined the sampling location order
newerorder<-c("Chicago_Camp","Near_Delabarre_Creek","Wilder","Headwaters","Hayes","Lost_River","Elkhorn","Long_Creek","Geyser","Whiskey_Bend",
              "Cat_Creek","Boulder","Glines","Altaire","Griff_Creek","Hughes_Creek","Campground_Creek","Madison_Creek",
              "Little_River","Indian_Creek","Aldwell","South_Branch_Little_River","Elwha_River")

#### Define the populations order
neworder<-c("AD","SBLR","ID","BD")

#### Ordered upriver to downriver
Pre<-arrange(transform(df,Location=factor(Location,levels=neworder)),site.x=factor(site.x,levels=newerorder))
Pre<-arrange(transform(Pre,site.x=factor(site.x,levels=newerorder)),site.x)

#### Plot K = 4 in ggplot pre-dam removal
# Organize individuals by sampling location and split plot by population
pdf("K_4_Pre.pdf",height=15,width=45)
plota<-ggplot(Pre,aes(x=sample.id,y=Ancestry),group=Variable)+theme(legend.position="none")
plota<-plota+geom_bar(aes(y=Ancestry,x=sample.id,fill=Variable),stat="identity")
plota<-plota+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=40,face="bold"))
plota<-plota+theme(axis.text.x=element_blank(),axis.title.x=element_text(size=40,face="bold"))
plota<-plota+ylab(label="Ancestry")+xlab(label="Sample Locations")
plota<-plota+scale_fill_manual(values=colors_location)
plota<-plota+scale_color_manual(values=colors_location)
plota<-plota+facet_grid(~site.x+Location,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(1.5)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(1.5)))
plota<-plota+theme(legend.position="bottom")
plota
dev.off()
