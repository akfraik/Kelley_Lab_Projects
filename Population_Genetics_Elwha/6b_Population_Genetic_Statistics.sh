###########################################################################################################################
## 6b. Population genetics analyses
# Calculate Pairwise FST
###########################################################################################################################

#### Setup
module load bio/vcftools
module load bio/plink/1.90b6.6

#### Set up directory paths to where files are stored
BASE="/INSERT PATH TO BASE DIRECTORY"
INDIR="/INSERT PATH TO VCF FILES"
INFODIR="/INSERT PATH TO INFO"
OUTDIR="/INSERT PATH TO OUTPUT FILE LOCATION"
METADATA="${INFODIR}/INSERT PATH TO METADATA"

#### Now let's do the same for calculating FST for Pre/Post dams
# Run this script after running the 6a_Population_Genetics_Statistics.sh script

#### Compare AD to BD Pre
  vcftools --gzvcf "${INDIR}/all_81_vcf.vcf.gz" --weir-fst-pop "${INFODIR}/Pre_FST.Keep" \
   --weir-fst-pop "${INFODIR}/Post_FST.Keep" --out "${OUTDIR}/Pre_Post_All_FST"

#### Compare AD to ID Pre
  vcftools --gzvcf "${INDIR}/Elw_Pre.vcf.gz" --weir-fst-pop "${INFODIR}/AD.Keep" \
   --weir-fst-pop "${INFODIR}/ID.Keep" --out "${OUTDIR}/Pre_AD_ID_FST"

#### Compare BD to ID Pre
  vcftools --gzvcf "${INDIR}/Elw_Pre.vcf.gz" --weir-fst-pop "${INFODIR}/BD.Keep" \
    --weir-fst-pop "${INFODIR}/ID.Keep" --out "${OUTDIR}/Pre_BD_ID_FST"

#### Compare BD to SBLR Pre
  vcftools --gzvcf "${INDIR}/Elw_Pre.vcf.gz" --weir-fst-pop "${INFODIR}/BD.Keep" \
    --weir-fst-pop "${INFODIR}/SBLR.Keep" --out "${OUTDIR}/Pre_BD_SBLR_FST"

#### Compare AD to SBLR Pre
  vcftools --gzvcf "${INDIR}/Elw_Pre.vcf.gz" --weir-fst-pop "${INFODIR}/AD.Keep" \
    --weir-fst-pop "${INFODIR}/SBLR.Keep" --out "${OUTDIR}/Pre_AD_SBLR_FST"

#### Compare ID to SBLR Pre
  vcftools --gzvcf "${INDIR}/Elw_Pre.vcf.gz" --weir-fst-pop "${INFODIR}/ID.Keep" \
    --weir-fst-pop "${INFODIR}/SBLR.Keep" --out "${OUTDIR}/Pre_ID_SBLR_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/BD.Keep" \
    --weir-fst-pop "${INFODIR}/ID.Keep" --out "${OUTDIR}/Post_BD_ID_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_Pre_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Adult_2015.Keep" --out "${OUTDIR}/Post_Adult_Pre_2015_2015_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_Pre_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Adult_2016.Keep" --out "${OUTDIR}/Post_Adult_Pre_2015_2016_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_Pre_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Adult_2017.Keep" --out "${OUTDIR}/Post_Adult_Pre_2015_2017_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_Pre_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2016.Keep" --out "${OUTDIR}/Post_Adult_Pre_2015_Juvenile_2016_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_Pre_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2017.Keep" --out "${OUTDIR}/Post_Adult_Pre_2015_Juvenile_2017_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Adult_2016.Keep" --out "${OUTDIR}/Post_Adult_2015_2016_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Adult_2017.Keep" --out "${OUTDIR}/Post_Adult_2015_2017_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2016.Keep" --out "${OUTDIR}/Post_Adult_2015_Juvenile_2016_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2015.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2017.Keep" --out "${OUTDIR}/Post_Adult_2015_Juvenile_2017_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2016.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Adult_2017.Keep" --out "${OUTDIR}/Post_Adult_2016_2017_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2016.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2016.Keep" --out "${OUTDIR}/Post_Adult_2016_Juvenile_2016_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2016.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2017.Keep" --out "${OUTDIR}/Post_Adult_2016_Juvenile_2017_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2017.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2016.Keep" --out "${OUTDIR}/Post_Adult_2017_Juvenile_2016_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Adult_2017.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2017.Keep" --out "${OUTDIR}/Post_Adult_2017_Juvenile_2017_FST"

#### Compare BD to ID Post
  vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --weir-fst-pop "${INFODIR}/Post_Juvenile_2016.Keep" \
    --weir-fst-pop "${INFODIR}/Post_Juvenile_2017.Keep" --out "${OUTDIR}/Post_Juvnile_2016_2017_FST"
