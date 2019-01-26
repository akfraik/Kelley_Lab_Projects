###########################################################################################################################
## 11c. SNP2GO (R Script)
###########################################################################################################################

## First step set the working directory and load the required package
  #install.packages("hash")
  #source("http://bioconductor.org/biocLite.R")
  #biocLite("goProfiles")
  #biocLite("GenomicRanges")
  #install.packages(pkgs="SNP2GO_1.0.5.tar.gz", type="source") 
  library(SNP2GO)
  
##Read in the VCF files
  snps<-read.delim("/Bed_Files/Annotated_Outliers_Bed_Variable.vcf",header=FALSE,comment.char="#")
  snps[,2]<-as.numeric(snps[,2])
  snps<-GRanges(seqnames=snps[,1],ranges=IRanges(snps[,2],snps[,2]))

## The manual recommends creating two lists: the candidates and non-candidates but I would like to take two different approaches
 #First use the non-candidate set as the entire set of SNPs, including our candidates
 #The second would be to create a set of non-candidates from the entire set sans the candidates
  #This means we need to end up generating two sadditional VCF files
  #Defining the candidates SNPs
 #Input files for candidate SNps to come from Gene_Annotations_Pre_VCF.sh script
  cand<-read.delim("/Bed_Files/Annotated_Outliers_Bed_Variable.vcf",header=FALSE,comment.char="#")
  cand[,2]<-as.numeric(cand[,2])
  cand<-GRanges(seqnames=cand[,1],ranges=IRanges(cand[,2],cand[,2]))

## Defining the non-candidates SNPs that are sans candidates
  noncand<-read.delim("/Bed_Files/Sans_Outliers_Bed_Variable.vcf",header=FALSE,comment.char="#")
  noncand[,2]<-as.numeric(noncand[,2])
  noncand<-GRanges(seqnames=noncand[,1],ranges=IRanges(noncand[,2],noncand[,2]))

## Case 1: Using a GFF file with candidates and full set 
 # Make sure you set the path to your directory containing your reference genome
  reference<-snp2go(gtf="/Sarcophilus_harrisii.DEVIL7.0.88.gtf.gz",candidateSNPs=cand,noncandidateSNPs=snps,FDR=0.05,runs=10000,extension=50)

## Get the first of the enriched GO terms then print all the regions with at least one enriched GO term
  gff.significant.terms<-reference$enriched$GO
  first.term<-gff.significant.terms[1]
  print(reference$regions)
  print(reference$regions[unlist(as.list(reference$go2ranges[["regions"]][first.term]))])

## Although version 2 seems more complicated, it allows to get the regions associated with more than one term. 
 # In the following example, all regions associated with the first ten enriched GO terms (gff.significant.terms[1:10]) are printed:
  print(reference$regions[unlist(as.list(reference$go2ranges[["regions"]][gff.significant.terms[1:10]]))])


## Print candidate and non-candidate SNPs associated with 1st enriched GO terms: 
 #Like for GO regions, there are also two possibilities to do that:
  print(cand[unlist(as.list(reference$go2ranges[["candidates"]][first.term]))])
  print(noncand[unlist(as.list(reference$go2ranges[["noncandidates"]][first.term]))])

## Get the number of informative candidates of the GTF analysis and print all GO terms associated to at least one gene in the
 # Store the results in tab-seperated files, GFF analysis:
  list<-reference$informative.candidate.snps
  print(reference$goterms)
  write.table(file="snp2go_pre_gtf.tsv",reference$enriched,sep="\t",row.names=F)
