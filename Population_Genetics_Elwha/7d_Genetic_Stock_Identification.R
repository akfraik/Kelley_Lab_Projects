###########################################################################################################################
## 7d. Genetic Stock Identification
# Plot the RUBIAS output
###########################################################################################################################

#### Set your working directory
setwd("/Users/workingdirectory")

#### Load all the necessary R packages and install if necessary
library(ggplot2)
library(dplyr)
library(vcfR)
library(stringr)
library(gridExtra)
library(grid)
library(tidyr)
library(ggdendro)
library(reshape2)

#### Read in and subset the metadata file for what you need
pop_code<-read.csv("/Users/alexandra.fraik/Documents/Elw_All/Pre_Post_Steelhead/GSI/GSI_07_20/GSI_DAPC_FST/pop_code_GSI.csv")
colnames(pop_code)[1]<-"sample.id"
pop_code<-pop_code[c("sample.id","Time","year","broodYearEst","lifestage.x","Species.1","run.timing","source.Lat","source.Long","GSI")]
pop_code<-pop_code[pop_code$Time=="Post",]

#### Let's read in the GSI files
location<-read.csv("Rep_Indiv_Ests_Location.csv")
df<-merge(location,pop_code,by="sample.id")
colors_location<-c("brown1","blue2","darkgreen","limegreen")
names(colors_location)<-levels(df$repunit)
colors_l<-data.frame(repunit=names(colors_location),colors_location=colors_location)
pop_code<-merge(df,colors_l,by="repunit")
colors_location<-as.character(pop_code$colors_location)
names(colors_location)<-as.character(pop_code$repunit)
newerorder<-c("AD","ID","SBLR","ID")

#### Subset adults by return year
pre<-pop_code[pop_code$GSI=="Adult_Pre_2015",]
fifteen<-pop_code[pop_code$GSI=="Adult_2015",]
sixteen<-pop_code[pop_code$GSI=="Adult_2016",]
seventeen<-pop_code[pop_code$GSI=="Adult_2017",]
df.list<-list(pre,fifteen,sixteen,seventeen)

#### Loop through function to reorder individuals from increasing to decreasing reporting unit ancestry proportions
neworder<-c("AD","ID","SBLR","BD")
res<-lapply(df.list, function(x) {
  test <- x %>% 
    group_by(sample.id) %>% 
    top_n(1, rep_pofz) %>%
    arrange(repunit=factor(repunit,levels=neworder)) %>%
    ungroup() %>%
    mutate(order = row_number())
  test$index<-as.numeric(row.names(test))
  order<-test[c("sample.id","index")]
  test<-merge(x,order,by="sample.id")
  x<-test[with(test,order(index)),]
  x } )
names(res)<-c("pre","fifteen","sixteen","seventeen")
list2env(res,envir=.GlobalEnv)
seventeen<-arrange(transform(seventeen,repunit=factor(repunit,levels=newerorder)),repunit)

#### GSI stock output by capture date
pdf("GSI_Year_Run_Timing_Adult.pdf",width=30,height=18)
plota<-ggplot(pre,aes(x=index,y=rep_pofz),group=repunit)
plota<-plota+geom_bar(stat="identity",aes(y=rep_pofz,x=index,fill=repunit,color=factor(run.timing)),size=1.5)
plota<-plota+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plota<-plota+scale_color_manual(name="Reference Sampling Sites",values=c("NA"))
plota<-plota+theme_bw()
plota<-plota+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_blank())
plota<-plota+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plota<-plota+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plota<-plota+theme(legend.position="none")
plota<-plota+ylab(label="")
plota<-plota+scale_x_continuous(expand=c(0, 0))
plotb<-ggplot(fifteen,aes(x=index,y=rep_pofz),group=repunit)
plotb<-plotb+geom_bar(stat="identity",aes(y=rep_pofz,x=index,fill=repunit,color=factor(run.timing)),size=1.5)
plotb<-plotb+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plotb<-plotb+scale_color_manual(name="Reference Sampling Sites",values=c("black","NA"))
plotb<-plotb+theme_bw()
plotb<-plotb+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_blank())
plotb<-plotb+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plotb<-plotb+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plotb<-plotb+theme(legend.position="none")
plotb<-plotb+ylab(label="")
plotb<-plotb+scale_x_continuous(expand=c(0,0))
plotc<-ggplot(sixteen,aes(x=index,y=rep_pofz),group=repunit)
plotc<-plotc+geom_bar(stat="identity",aes(y=rep_pofz,x=index,fill=repunit,color=factor(run.timing)),size=1.5)
plotc<-plotc+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plotc<-plotc+scale_color_manual(name="Reference Sampling Sites",values=c("NA"))
plotc<-plotc+theme_bw()
plotc<-plotc+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_blank())
plotc<-plotc+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plotc<-plotc+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plotc<-plotc+theme(legend.position="none")
plotc<-plotc+ylab(label="")
plotc<-plotc+scale_x_continuous(expand=c(0,0))
plotd<-ggplot(seventeen,aes(x=index,y=rep_pofz),group=repunit)
plotd<-plotd+geom_bar(stat="identity",aes(y=rep_pofz,x=index,fill=repunit,color=factor(run.timing)),size=1.5)
plotd<-plotd+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plotd<-plotd+scale_color_manual(name="Reference Sampling Sites",values=c("black","NA"))
plotd<-plotd+theme_bw()
plotd<-plotd+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_blank())
plotd<-plotd+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plotd<-plotd+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plotd<-plotd+theme(legend.text=element_text(size=rel(1.5)))
plotd<-plotd+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plotd<-plotd+theme(legend.position="bottom")
plotd<-plotd+ylab(label="")
plotd<-plotd+scale_x_continuous(expand=c(0,0))
grid.arrange(plota,plotb,plotc,plotd,left=textGrob("Inferred Genetic Ancestry",gp=gpar(fontsize=25,fontface='bold'),rot=90),nrow=4,ncol=1)
dev.off()

#### Let's try maybe binding these seqeuences in a dataframe
df<-rbind(pre,fifteen,sixteen,seventeen)
neworder<-c("Adult_Pre_2015","Adult_2015","Adult_2016","Adult_2017")
df<-arrange(transform(df,GSI=factor(GSI,levels=neworder,labels=c("Adults Prior to 2015","Adults 2015","Adults 2016","Adults 2017"))))
newerorder<-c("AD","ID","SBLR","BD")
df<-arrange(transform(df,repunit=factor(repunit,levels=newerorder)),repunit)

#### Then plot an alternative figure
pdf("GSI_Year_Run_Timing_Adult_2.pdf",width=30,height=14)
plotd<-ggplot(df,aes(x=index,y=rep_pofz),group=repunit)
plotd<-plotd+geom_bar(aes(y=rep_pofz,x=index,fill=repunit),stat="identity")
plotd<-plotd+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plotd<-plotd+ylab(label="Inferred Genetic Ancestry")
plotd<-plotd+theme_bw()
plotd<-plotd+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_text(size=20,face="bold"))
plotd<-plotd+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plotd<-plotd+facet_grid(GSI~.,scales="fixed")
plotd<-plotd+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plotd<-plotd+theme(legend.text=element_text(size=rel(1.5)))
plotd<-plotd+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plotd<-plotd+theme(legend.position="bottom")
plotd<-plotd+scale_x_continuous(expand=c(0, 0))
plotd
dev.off()

#### subset juveniles by year
j_sixteen<-pop_code[pop_code$GSI=="Juvenile_2016",]
j_seventeen<-pop_code[pop_code$GSI=="Juvenile_2017",]
df.list<-list(j_sixteen,j_seventeen)

#### Loop through function
neworder<-c("AD","ID","SBLR","BD")
res<-lapply(df.list, function(x) {
  test <- x %>% 
    group_by(sample.id) %>% 
    top_n(1, rep_pofz) %>%
    arrange(repunit=factor(repunit,levels=neworder)) %>%
    ungroup() %>%
    mutate(order = row_number())
  test$index<-as.numeric(row.names(test))
  order<-test[c("sample.id","index")]
  test<-merge(x,order,by="sample.id")
  x<-test[with(test,order(index)),]
  x } )
names(res)<-c("j_sixteen","j_seventeen")
list2env(res,envir=.GlobalEnv)
j_seventeen<-arrange(transform(j_seventeen,repunit=factor(repunit,levels=newerorder)),repunit)

#### GSI stock output by capture date
pdf("GSI_Year_Run_Timing_Juvenile.pdf",width=30,height=14)
plotc<-ggplot(j_sixteen,aes(x=index,y=rep_pofz),group=repunit)
plotc<-plotc+geom_bar(stat="identity",aes(y=rep_pofz,x=index,fill=repunit))
plotc<-plotc+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plotc<-plotc+theme_bw()
plotc<-plotc+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_blank())
plotc<-plotc+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plotc<-plotc+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plotc<-plotc+theme(legend.position="none")
plotc<-plotc+ylab(label="")
plotc<-plotc+scale_x_continuous(expand=c(0,0))
plotd<-ggplot(j_seventeen,aes(x=index,y=rep_pofz),group=repunit)
plotd<-plotd+geom_bar(stat="identity",aes(y=rep_pofz,x=index,fill=repunit))
plotd<-plotd+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plotd<-plotd+theme_bw()
plotd<-plotd+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_blank())
plotd<-plotd+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plotd<-plotd+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plotd<-plotd+theme(legend.text=element_text(size=rel(1.5)))
plotd<-plotd+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plotd<-plotd+theme(legend.position="bottom")
plotd<-plotd+ylab(label="")
plotd<-plotd+scale_x_continuous(expand=c(0,0))
grid.arrange(plotc,plotd,left=textGrob("Inferred Genetic Ancestry",gp=gpar(fontsize=25,fontface='bold'),rot=90),nrow=2,ncol=1)
dev.off()

#### Let's try maybe binding these seqeuences in a dataframe
df<-rbind(j_sixteen,j_seventeen)
neworder<-c("Juvenile_2016","Juvenile_2017")
df<-arrange(transform(df,GSI=factor(GSI,levels=neworder,labels=c("Juveniles 2016","Juveniles 2017"))))
newerorder<-c("AD","ID","SBLR","BD")
df<-arrange(transform(df,repunit=factor(repunit,levels=newerorder)),repunit)

#### Then plot an alternative figure
pdf("GSI_Year_Run_Timing_Juvenile_2.pdf",width=30,height=14)
plotd<-ggplot(df,aes(x=index,y=rep_pofz),group=repunit)
plotd<-plotd+geom_bar(aes(y=rep_pofz,x=index,fill=repunit),stat="identity")
plotd<-plotd+scale_fill_manual(name="Reference Sampling Sites",values=colors_location)
plotd<-plotd+ylab(label="Inferred Genetic Ancestry")
plotd<-plotd+theme_bw()
plotd<-plotd+theme(axis.text.y=element_text(color="black",size=20),axis.title.y=element_text(size=20,face="bold"))
plotd<-plotd+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plotd<-plotd+facet_grid(GSI~.,scales="fixed")
plotd<-plotd+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plotd<-plotd+theme(legend.text=element_text(size=rel(1.5)))
plotd<-plotd+theme(legend.title=element_text(face="bold",size=rel(2.5)))
plotd<-plotd+theme(legend.position="bottom")
plotd<-plotd+scale_x_continuous(expand=c(0, 0))
plotd
dev.off()

#### Let's read in the GSI files for collections
location<-read.csv("Rep_Mix_Ests_Location.csv")
colors_location<-c("brown1","blue2","darkgreen","limegreen")
location<-location[c("mixture_collection","repunit","indiv_repprop")]
names(colors_location)<-levels(location$repunit)
colors_l<-data.frame(repunit=names(colors_location),colors_location=colors_location)
pop_code<-merge(location,colors_l,by="repunit")
colors_location<-as.character(pop_code$colors_location)
names(colors_location)<-as.character(pop_code$repunit)
pop_code<-pop_code[!pop_code$mixture_collection=="Pre_Dam_Removal",]
levels(pop_code$mixture_collection)[levels(pop_code$mixture_collection)=="Adult_Pre_2015"]<-"Adults_Pre 2015"
pop_code <- pop_code %>% separate(mixture_collection, c("lifestage.x","year"),"_")
pop_code<-as.data.frame(pop_code)
neworder<-c("Pre 2015","2015","2016","2017")
df<-arrange(transform(pop_code,year=factor(year,levels=neworder,labels=c("Prior to 2015","2015","2016","2017"))))
newerorder<-c("AD","ID","SBLR","BD")
df<-arrange(transform(df,repunit=factor(repunit,levels=newerorder)),repunit)
df$lifestage.x<-as.factor(df$lifestage.x)
levels(df$lifestage.x)[levels(df$lifestage.x)=="Adults"]<-"Adult"

#### GSI stock output by capture date
names<-c('Pre_2015'="Prior to 2015",'2015'="2015",'2016'="2016",'2017'="2017")
pdf("GSI_Year_Run_Timing_Mixtures.pdf",width=15,height=10)
plota<-ggplot(df,aes(fill=repunit,y=indiv_repprop,x=year))
plota<-plota+geom_bar(position="fill",stat="identity")
plota<-plota+scale_fill_manual(values=colors_location,name="Reporting Unit")
plota<-plota+scale_color_manual(values=colors_location,name="Reporting Unit")
plota<-plota+xlab(label="Life-history Cohorts")
plota<-plota+ylab(label="Inferred Reporting Unit Proportion")
plota<-plota+theme_bw()
plota<-plota+theme(axis.text.y=element_text(color="black",size=15),axis.title.y=element_text(size=20,face="bold"))
plota<-plota+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())
plota<-plota+facet_grid(lifestage.x~year,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text=element_text(face="bold",size=20),strip.background=element_rect(fill="white",color="white"))
plota<-plota+theme(legend.text=element_text(size=rel(1.5)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(2)))
plota
dev.off()
