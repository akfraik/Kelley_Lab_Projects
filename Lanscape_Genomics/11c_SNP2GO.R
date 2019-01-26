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
  snps<-read.delim("/project/devils_alexf/try3/data/info/Rapture/rap_pre.vcf",header=FALSE,comment.char="#")
  snps[,2]<-as.numeric(snps[,2])
  snps<-GRanges(seqnames=snps[,1],ranges=IRanges(snps[,2],snps[,2]))

#The manual recommends creating two lists: the candidates and non-candidates but I would like to take two different approaches
#First use the non-candidate set as the entire set of SNPs, including our candidates
#The second would be to create a set of non-candidates from the entire set sans the candidates
#This means we need to end up generating two sadditional VCF files

#Defining the candidates SNPs
#Input files for candidate SNps to come from Gene_Annotations_Pre_VCF.sh script
        cand <- read.delim("/project/devils_alexf/try3/data/info/Rapture/rap_pre_snp.vcf",header=FALSE,comment.char="#")
        cand[,2] <- as.numeric(cand[,2])
        cand <- GRanges(seqnames=cand[,1],ranges=IRanges(cand[,2],cand[,2]))

#Defining the non-candidates SNPs that are sans candidates
#       noncand <- read.delim("/project/devils_alexf/try3/data/Rapture/VCFtools/rap_pre_minus_SNPs.vcf",header=FALSE,comment.char="#")
#       noncand[,2] <- as.numeric(noncand[,2])
#       noncand <- GRanges(seqnames=noncand[,1],ranges=IRanges(noncand[,2],noncand[,2]))

        # Case 1: Using a GFF file with candidates and full set 
#       x <- snp2go(gff="/project/devils_alexf/try3/data/Rapture/Bed_Files/Sarcophilus_harrisii.DEVIL.gff3.gz",
 #            candidateSNPs=cand,
  #           noncandidateSNPs=snps)
