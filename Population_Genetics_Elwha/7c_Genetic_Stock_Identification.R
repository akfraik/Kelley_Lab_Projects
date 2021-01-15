###########################################################################################################################
## 7c. Genetic Stock Identification
# Run Rubias
###########################################################################################################################

#### Set your working directory
setwd("/Users/workingdirectory")

#### Load all the necessary R packages and install if necessary
library(rubias)
library(tidyverse)
library(dplyr)

#### Name columns for loci
prefix<-"locus"
test<-rep(seq(1:265),each=2)
names<-c("a","b")
testa<-rep(names,times=265)
loci.names<-paste(prefix,test,names,sep="_")

#### Read in reference data file
pop_code<-read.csv("/METADATA_FILE")
colnames(metadata)[1]<-"sample.id"

#### Subset metadata files
metadata<-metadata[c("sample.id","Location","site.x","GSI")]
write.csv(metadata,"Post_ID_File.csv")

#### Read in GSI Reference file
ref<-read.csv("GSI_Ref.csv",header=FALSE)
ref<-ref[c(-1,-3,-4,-5,-6)]

#### Subset just the loci from the data frame
loci<-ref[c(2:531)]

#### Subset just the loci from the data frame
colnames(loci)<-loci.names
here<-as.data.frame(ref[c(1)])
colnames(here)[1]<-"sample.id"
metadata_test<-merge(here,metadata,by="sample.id")
sample_type<-as.character(rep("reference",times=479))

#### Generate the ultimate reference data frame
df_ref<-data.frame(sample_type=sample_type,repunit=metadata_test$Location,collection=metadata_test$site.x,indiv=metadata_test$sample.id,loci)
df_ref$sample_type<-as.character(df_ref$sample_type)
df_ref$repunit<-as.character(df_ref$repunit)
df_ref$collection<-as.character(df_ref$collection)
df_ref$indiv<-as.character(df_ref$indiv)

#### Assess the self-assignments from the reference
sa<-self_assign(reference=df_ref,gen_start_col=5)
head(sa,n=100)
# Log likelihood is the log probability of the fish's genotype given it is from the collection using the leave-out approach
max_log_likelihood <- sa %>%
  group_by(collection, repunit) %>%
  filter(log_likelihood == max(log_likelihood)) 
# The posterior probability is that of assigning the fish to the collection given an equal prior on every collection in the reference
test<-as.data.frame(max_log_likelihood)
# Z-scores can be used to assess the distribution of the z-score statistic for fish from known reference populations
test<-test[c(1:8)]

#### Reorder the populations by sampling site from upriver -> downriver
newerorder<-c("Chicago_Camp","Wilder","Hayes","Elkhorn","Geyser","Cat_Creek","Altaire","Campground_Creek",
              "Madison_Creek","Little_River","Indian_Creek","Aldwell","South_Branch_Little_River","Elwha_River")

#### Now let's summarize the leave-one-out approach results by the individual mixture inferred collection 
# The scaled_likelihood is the posterior prob of assigning the fish to the inferred_collection given an equal prior on every collection in the reference
# The repu scaled likelihood is the summation of the probabilities of assigning a fish to a collection per repunit
sa_to_indiv<-sa %>%
  select(indiv,collection,inferred_collection,repunit,inferred_repunit,z_score,log_likelihood,scaled_likelihood)

#### Write the output from the self-assign LOO simulations to file for individuals
write.csv(sa_to_indiv,"Assessment_GSI_References_Individuals.csv")

#### Now let's summarize the leave-one-out approach by the individual mixture inferred collection 
# The scaled likelihood is the likehood of assigning a fish back to the inferred collection
# The repu scaled likelihood is the summation of the likelihoods of assigning a fish back to the inferred repunit
sa_to_collec <- sa_to_indiv %>% 
  group_by(collection, inferred_collection) %>% 
  summarise(collec_mean_scaled_like = mean(scaled_likelihood)) %>% 
  spread(inferred_collection,collec_mean_scaled_like) %>% 
  arrange(factor(collection,levels=newerorder),collection) %>% 
  select(as.vector(newerorder))

#### Write the output from the self-assign LOO simulations to file for collections
write.csv(sa_to_collec,"Assessment_GSI_References_Collections.csv")

#### Now let's summarize the leave-one-out approach by the individual mixture inferred collection 
sa_to_repu <- sa %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))

#### Write the output from the self-assign LOO simulations to file for reporting units
write.csv(sa_to_repu,"Assessment_GSI_References_Repunit.csv")

#### Sumarize by collection by taking the average of the scaled likelihoods per inferred rep unit
# Likelihood of membership in inferred collection)
sa_to_collec <- sa %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood)) %>%
  group_by(collection, repunit, inferred_repunit) %>%
  summarise(repu_mean_scaled_like = mean(repu_scaled_like)) %>% 
  spread(inferred_repunit,repu_mean_scaled_like) %>% 
  arrange(factor(collection,levels=newerorder),collection)
  
#### Write the output from the self-assign LOO simulations to file for inferred reporting units
write.csv(sa_to_collec,"Assessment_GSI_References_Inferred_Repunits.csv")

#### Filter out potentially problem individuals
problem_individuals<- sa %>% 
  filter(is.na(scaled_likelihood)) %>%
  distinct(indiv,.keep_all = FALSE)
problem_individuals<-data.frame(problem_individuals)
problem_individuals<-as.vector(problem_individuals$indiv)
df_ref<-filter(df_ref,!(indiv %in% problem_individuals))

#### Let's also use the LOO approach to simulate different mixtures, now, call those repunits and then sum up over reporting units
# Resampling_unit = "gene_copies"
omy_sims<-assess_reference_loo(reference=df_ref,gen_start_col=5,reps=500,mixsize=50,return_indiv_posteriors = TRUE)
mix_props<-omy_sims$mixing_proportions
tmp <- mix_props %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))

#### Set the colors for the repunits
colors_location<-c("brown1","blue2","darkgreen","limegreen")
names(colors_location)<-c("AD","BD","ID","SBLR")
names<-c('AD'="Above the Dams",'BD'="Below the Dams",'ID'="In Between the Dams",'SBLR'="South Branch of the Little River")
neworder<-c("AD","ID","SBLR","BD")
tmp<-arrange(transform(tmp,repunit=factor(repunit,levels=neworder)),repunit)

#### So, really, to determine these fish, I should just determine how many SDs from the mean a fish has to be (for its PofZ) 
# I am going to compare the distributions for the number of fish that have a PofZ for their max posterior membership
sa_to_repu_z_scores <- sa %>%
  group_by(indiv, collection, inferred_collection, repunit, inferred_repunit, z_score) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))
write.csv(sa_to_repu_z_scores,"Individual_Posterior_Locations_Pre_Dam_Removal.csv")

#### PofZ is the mean posterior membership over the MCMC of the posterior probability that the individual originated from the COLLECTION, not repunit
indiv_props<-omy_sims$indiv_posteriors

#### Plot the mixutre and simulated proportions per repunit
# Plot out the simulated LOO refernce individuals' assignment to reporting unit and compare to real assignment
# Simulated_repunit is the reporting unit the individual was simulated from 
# Simulated_collection is the collection the simulated genotype came from 
pdf(file="Simulations_Reference_Repunits.pdf",width=8,height=8)
plota<-ggplot(tmp,aes(x=true_repprop,y=reprop_posterior_mean,colour=repunit))+geom_point() 
plota<-plota+scale_color_manual(name="Reference Reporting Units",values=colors_location)
plota<-plota+geom_abline(intercept=0,slope=1)
plota<-plota+facet_wrap(~repunit,labeller=as_labeller(names))
plota<-plota+labs(x="True Posterior Mixture Proprtion",y="Simulated Posterior Mixture Proprtion")
plota<-plota+scale_x_continuous(limits=c(0,1))
plota<-plota+scale_y_continuous(limits=c(0,1))
plota<-plota+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plota<-plota+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(1)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(1.5)))
plota<-plota+theme(legend.position="none")
plota
dev.off()

#### We can see what the distribution of posteriors to the correct reporting unit is for fish from the different simulated collections
# We’ll do that with a boxplot, coloring by repunit, first aggregating over reporting units
repu_pofzs <- indiv_props %>%
  filter(repunit == simulated_repunit) %>%
  group_by(iter, indiv, simulated_collection, repunit) %>%  
  summarise(repu_PofZ = sum(PofZ)) %>%
  ungroup() %>%
  arrange(repunit, simulated_collection) %>%
  mutate(simulated_collection = factor(simulated_collection, levels = unique(simulated_collection)))
df<-arrange(transform(repu_pofzs,simulated_collection=factor(simulated_collection,levels=newerorder),simulated_collection))
repu_pofzs<- df %>% drop_na(repu_PofZ)

#### Also get the number of simulated individuals from each collection
num_simmed <- indiv_props %>%
  group_by(iter, indiv) %>%
  slice(1) %>%
  ungroup() %>%
  count(simulated_collection)

#### Plot the mixutre and simulated proportions per individual per collection
# Note, the last few steps make simulated collection a factor so that collections within the same repunit are grouped together in the plot.
pdf(file="Simulations_Reference_Collections.pdf",width=8,height=8)
plota<-ggplot(repu_pofzs, aes(x = simulated_collection, y = repu_PofZ))
plota<-plota+geom_boxplot(aes(colour = repunit))
plota<-plota+labs(x="Simulated Collection",y="Reporting Unit Posterior Probability")
plota<-plota+geom_text(data = num_simmed, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 3)
plota<-plota+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5))
plota<-plota+ylim(c(0, 1.05))
plota<-plota+scale_x_discrete(labels=newerorder)
plota<-plota+scale_color_manual(name="Reference Reporting Units",values=colors_location)
plota<-plota+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plota<-plota+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(1)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(1.5)))
plota<-plota+theme(legend.position="none")
plota
dev.off()

#### Read in GSI mixture data file
mix<-read.csv("GSI_Mix.csv",header=FALSE)
mix<-mix[c(-1,-3,-4,-5,-6)]

#### Subset just the loci from the data frame
loci<-mix[c(2:531)]

#### Subset just the loci from the data frame
colnames(loci)<-loci.names
here<-as.data.frame(mix[c(1)])
colnames(here)[1]<-"sample.id"
metadata_test<-merge(here,metadata,by="sample.id")
sample_type<-as.character(rep("mixture",times=558))
repunit<-as.character(rep("here",times=558))

#### Generate the ultimate reference data frame
df_mix<-data.frame(sample_type=sample_type,repunit=repunit,collection=metadata_test$GSI,indiv=metadata_test$sample.id,loci)
df_mix$sample_type<-as.character(df_mix$sample_type)
df_mix$repunit<-as.character(df_mix$repunit)
df_mix$collection<-as.character(df_mix$collection)
df_mix$indiv<-as.character(df_mix$indiv)
df_mix[df_mix=="here"]<-NA

#### For the first round, let's try this
# Method "PB" runs the full Bayesian model and looks for 
mix_est<-infer_mixture(reference=df_ref,mixture=df_mix,gen_start_col=5,method="PB")
lapply(mix_est, head)
location_post<-as.data.frame(mix_est$indiv_posteriors)
location_post<-location_post[c(1:7)]
write.csv(location_post,"Individual_Posterior_Location.csv")

#### Aggregating collections into reporting units for mixing proportions adding mixing proportions over collections in the repunit
rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  

#### Aggregating collections into reporting units for mixing proportions adding mixing proportions over collections in the repunit
rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  

#### Now let's try seeing if our distribution of samples are significantly different from the distribution of Z-scores for each reportig units
# Just a reminder, the Z scores for the reference samples is in here
# We need to subset based on collection, not repunit since the log-likelihood is the likelihood of collection
# But we also only want ones where the collection and the inferred collection match
ref_z_scores <- sa_to_repu_z_scores %>% 
  filter(collection == inferred_collection)

#### To see if individuals are "non-elwha" by checking to see if their Z-score falls out of the distribution of Z-scores for the reference 
# collection I think the best way to do this is to compare the maximum a-posteriori repunit of each individual in the dataset (so this should be in map_rows) 
# to the Z scores that make up the distribution for the AD/ID/BD get the maximum-a-posteriori population for each individual
mix_z_scores <- mix_est$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

#### Now since I am too lazy and it is too late at night for this b**** I am going to write it all out for each population....
Little_River_mix<- mix_z_scores %>% filter(collection == "Little_River")
Little_River_ref <- ref_z_scores %>% filter(collection == "Little_River")
Little_River_mix$Alien<-ifelse(sapply(Little_River_mix$z_score, function(x)
  any(min(Little_River_ref$z_score) >= x | max(Little_River_ref$z_score) <= x)),"YES", "NA")

#### Too many lines of code, so pretend I already ran it

#### There are X fish that are outside of the bounds, but they are so close to not being outside of the bounds that it just seems unlikely
df<-rbind(Little_River_mix,South_Branch_Little_River_mix,Indian_Creek_mix,Altaire_mix,Geyser_mix,Elkhorn_mix,Hayes_mix,
          Wilder_mix,Chicago_Camp_mix,Elwha_River_mix,Madison_Creek_mix,Aldwell_mix,Cat_Creek_mix,Campground_Creek_mix)

#### Filter out all of the fish that were "alien" or fell outside the bounds
alien_fish <- df %>% filter(Alien == "YES") %>% 
  filter(PofZ < 0.90)

#### Create alien fish vector and write to output file
alien_fish<-data.frame(mixutre_collection=alien_fish$mixture_collection,indiv=alien_fish$indiv,repunit=alien_fish$repunit,
                       collection=alien_fish$collection,PofZ=alien_fish$PofZ,log_likelihood=alien_fish$log_likelihood)
write.csv(alien_fish,"Alien_GSI_Fish.csv")
alien_fish<-as.vector(alien_fish$indiv)

#### Now let's go back and remove those "alien" fish from the mixture collection 
mix<-read.csv("GSI_Mix.csv",header=FALSE)
mix<-filter(mix,!(V2 %in% alien_fish))
value<-length(mix$V2)
mix<-mix[c(-1,-3,-4,-5,-6)]

#### Subset just the loci from the data frame
loci<-mix[c(2:531)]

#### Subset just the loci from the data frame
colnames(loci)<-loci.names
here<-as.data.frame(mix[c(1)])
colnames(here)[1]<-"sample.id"
metadata_test<-merge(here,metadata,by="sample.id")

#### THE TWO LINES BELOW NEED TO BE CHANGED TO MATCH THE NUMBER OF FISH RETAINED AFTER FILTERING OUT ALIEN FISH
sample_type<-as.character(rep("mixture",times=value))
repunit<-as.character(rep("here",times=value))

#### Generate the ultimate reference data frame
df_mix<-data.frame(sample_type=sample_type,repunit=repunit,collection=metadata_test$GSI,indiv=metadata_test$sample.id,loci)
df_mix$sample_type<-as.character(df_mix$sample_type)
df_mix$repunit<-as.character(df_mix$repunit)
df_mix$collection<-as.character(df_mix$collection)
df_mix$indiv<-as.character(df_mix$indiv)
df_mix[df_mix=="here"]<-NA

#### For the first round, let's try this
# Method "PB" runs the full Bayesian model and looks for 
mix_est<-infer_mixture(reference=df_ref,mixture=df_mix,gen_start_col=5,method="PB")
lapply(mix_est, head)
location_post<-as.data.frame(mix_est$indiv_posteriors)
location_post<-location_post[c(1:7)]

#### Creates bootstrap corrected CIs
mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi)) %>%
  left_join(mix_est$bootstrapped_proportions) %>%
  ungroup() %>%
  filter(mixture_collection == "Adult_2015") %>%
  arrange(desc(repprop)) %>%
  slice(1:10)

#### Aggregating collections into reporting units for mixing proportions adding mixing proportions over collections in the repunit
rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  
write.csv(rep_mix_ests,"Rep_Mix_Ests_Location.csv")

#### For individuals posteriors
rep_indiv_ests <- mix_est$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))
colnames(rep_indiv_ests)[2]<-"sample.id"
write.csv(rep_indiv_ests,"Rep_Indiv_Ests_Location.csv")


############################################### Juvenile 2016 ###############################################
#### Now we are going to create the posterior density 
# Now we are going to read in these files find the 3 most abundant:
# Relative density
top4 <- rep_mix_ests %>%
  filter(mixture_collection == "Juvenile_2016") %>% 
  arrange(desc(repprop)) %>%
  slice(1:4)

#### Check how many MCMC sweeps were done:
nsweeps <- max(mix_est$mix_prop_traces$sweep)

#### Keep only the collection of interest (Juvenile_2016 in this case) then discard the first 200 sweeps 
# as burn-in and then aggregate over reporting units and then keep only the 
# top 3 reporting units (as we only use 3) from above
trace_subset_Juvenile_2016 <- mix_est$mix_prop_traces %>%
  filter(mixture_collection == "Juvenile_2016", sweep > 200) %>%
  group_by(sweep, repunit) %>%
  summarise(repprop = sum(pi)) %>% 
  filter(repunit %in% top4$repunit)

#### Calculate pi = mixing proportion of each mixture collection
top4_cis <- trace_subset_Juvenile_2016 %>%
  group_by(repunit) %>%
  summarise(loCI = quantile(repprop, probs = 0.025),
            hiCI = quantile(repprop, probs = 0.975),
            Mean_Prob = mean(repprop))

#### Geerate Pi CIs and write to output file
top4_cis<-as.data.frame(top4_cis)
write.csv(top4_cis,"Juvenile_2016_CIs.csv")

#### Now let's plot those
pdf(file="Juvenile_2016_Trace.pdf")
plota<-ggplot(trace_subset_Juvenile_2016, aes(x = repprop, color = repunit))+geom_density()                                                                     
plota<-plota+scale_color_manual(name="Representative Populations",values=colors_location)
plota
dev.off()

#### Figure out if any of these individuals are not from any population
# get the maximum-a-posteriori population for each individual
map_rows <- mix_est$indiv_posteriors %>%
  filter(mixture_collection == "Juvenile_2016") %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()
map_rows_juvenile_2016<-as.data.frame(map_rows)

#### If everything is kosher, then we expect that the z-scores we see will be roughly normally distributed
normo <- tibble(z_score = rnorm(1e06))
# Plot the Z score figure
pdf(file="Z_Score_Juvenile_2016.pdf")
ggplot(map_rows_juvenile_2016, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")
dev.off()

#### Let's conduct a KS test to see if the distribution of values is significantly different from the normal distribution
output<-apply(map_rows[7],1,function(x) {
  ks<-ks.test(x,normo %>% pull(z_score),alterantive=two.sided)
  c(statistic=ks$statistic, p.value=ks$p.value)
  setNames(c(ks$statistic, ks$p.value), c("statistic", "p.value"))
})

#### Finally create a dataframe for this life-history cohort
juvenile_2016_df<-as.data.frame(t(output))
juvenile_2016_df$Individual<-map_rows$indiv
juvenile_2016_df$Reporting_Unit<-"Juvenile_2016"

##### Now repeat the above chunk of code for each life-history cohort or "collection" post-dam removal
# Then just a reminder, the Z scores for the reference samples is in hereif it changes drastically among 
hist(map_rows$PofZ,breaks=200,freq=FALSE)

#### let's also makea a figure with all of the reference populations now we can plot those:
here<-rbind(trace_subset_Adult_Pre_2015,trace_subset_Adult_2015,trace_subset_Adult_2016,trace_subset_Adult_2017,trace_subset_Juvenile_2016,trace_subset_Juvenile_2017)
names(colors)<-levels(as.factor(here$repunit))
colors_l<-data.frame(repunit=names(colors_location),colors=colors_location)
trace_subset<-merge(here,colors_l,by="repunit")

#### Make whole trace subset plot
pdf(file="Trace_Subset_All.pdf")
plota<-ggplot(trace_subset,aes(x=repprop,color=repunit))
plota<-plota+geom_density()                                                                     
plota<-plota+scale_color_manual(name="Representative Populations",values=colors_location)
plota
dev.off()

#### Now let's try seeing if our distribution of samples are significantly different from  the normal distribution
# Now let's bind these data frames & filter these guys out of the big data set and give this thing a try
test<-rbind(juvenile_2016_df,juvenile_2017_df,adult_pre_2015_df,adult_2015_df,adult_2016_df,adult_2017_df)
alien_fish <- test %>% filter(p.value < 0.01)
dim(alien_fish)
alien_fish_df <- mix_est$indiv_posteriors %>%
  group_by(indiv) %>%
  filter(indiv %in% alien_fish$Individual)
write.csv(test,"Non_Normal_Distributed_GSI_Fish.csv")
