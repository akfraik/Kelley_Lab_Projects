###########################################################################################################################
## 4b. Popualtion Genetics Statistics
###########################################################################################################################

############################# First portion of this script calculates heterozygosity #####################################
## Filtering options that decided "SNP" list (same as in plink files)
 # Directories where files are stored
INDIR="/Population_Genetics"
OUTDIR="/Population_Genetics/FST"
POPULATIONS="1 2 3 4 5 6 7"   # Replace numbers with your population names
KEEPDIR="/Plink_Input"

#Load vcf tools
  module load vcftools

## Use converted VCF pre-dftd file to start off this script from the last script
  gzip "${INDIR}/pre_vcf.recode.vcf" 
  gzip "${INDIR}/post_vcf.recode.vcf" 
  
## Compare population 1 to 2:7 for pre-disease, then post-disease
for i in 2 3 4 5 6 7;
do

 # $KEEPDIR/1.before - population 1 "keep" file for pre-disease  to the ith population
  vcftools --gzvcf "${INDIR}/pre_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/1.before" \
    --weir-fst-pop "$KEEPDIR/${i}.before" --out "$OUTDIR/1_${i}_pre_FST"
  
  # $KEEPDIR/1.after - population 1 "keep" file for post-disease to the ith population
   vcftools --gzvcf "${INDIR}/post_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/1.after" \
    --weir-fst-pop "$KEEPDIR/${i}.after" --out "$OUTDIR/1_${i}_post_FST"
    
 done
 
## Compare population 2 to 3:7 for pre-disease, then post-disease
for i in 3 4 5 6 7;
do

 # $KEEPDIR/2.before - population 2 "keep" file for pre-disease  to the ith population
  vcftools --gzvcf "${INDIR}/pre_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/2.before" \
    --weir-fst-pop "$KEEPDIR/${i}.before" --out "$OUTDIR/2_${i}_pre_FST"
  
  # $KEEPDIR/2.after - population 2 "keep" file for post-disease to the ith population
   vcftools --gzvcf "${INDIR}/post_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/2.after" \
    --weir-fst-pop "$KEEPDIR/${i}.after" --out "$OUTDIR/2_${i}_post_FST"
    
 done
 
## Compare population 3 to 4:7 for pre-disease, then post-disease
for i in 4 5 6 7;
do

 # $KEEPDIR/3.before - population 3 "keep" file for pre-disease  to the ith population
  vcftools --gzvcf "${INDIR}/pre_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/3.before" \
    --weir-fst-pop "$KEEPDIR/${i}.before" --out "$OUTDIR/3_${i}_pre_FST"
  
  # $KEEPDIR/3.after - population 3 "keep" file for post-disease to the ith population
   vcftools --gzvcf "${INDIR}/post_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/3.after" \
    --weir-fst-pop "$KEEPDIR/${i}.after" --out "$OUTDIR/3_${i}_post_FST"
    
 done
 
## Compare population 4 to 5:7 for pre-disease, then post-disease
for i in 5 6 7;
do

 # $KEEPDIR/4.before - population 4 "keep" file for pre-disease  to the ith population
  vcftools --gzvcf "${INDIR}/pre_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/4.before" \
    --weir-fst-pop "$KEEPDIR/${i}.before" --out "$OUTDIR/4_${i}_pre_FST"
  
  # $KEEPDIR/4.after - population 4 "keep" file for post-disease to the ith population
   vcftools --gzvcf "${INDIR}/post_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/4.after" \
    --weir-fst-pop "$KEEPDIR/${i}.after" --out "$OUTDIR/4_${i}_post_FST"
    
 done
 
## Compare population 5 to 6:7 for pre-disease, then post-disease
for i in 6 7;
do

 # $KEEPDIR/5.before - population 5 "keep" file for pre-disease  to the ith population
  vcftools --gzvcf "${INDIR}/pre_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/5.before" \
    --weir-fst-pop "$KEEPDIR/${i}.before" --out "$OUTDIR/5_${i}_pre_FST"
  
  # $KEEPDIR/5.after - population 5 "keep" file for post-disease to the ith population
   vcftools --gzvcf "${INDIR}/post_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/5.after" \
    --weir-fst-pop "$KEEPDIR/${i}.after" --out "$OUTDIR/5_${i}_post_FST"
    
 done
 
## Compare population 6 to 7 for pre-disease, then post-disease
  vcftools --gzvcf "${INDIR}/pre_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/6.before" \
    --weir-fst-pop "$KEEPDIR/7.before" --out "$OUTDIR/6_7_pre_FST"
   vcftools --gzvcf "${INDIR}/post_vcf.recode.vcf.gz" --weir-fst-pop "$KEEPDIR/6.after" \
    --weir-fst-pop "$KEEPDIR/7.after" --out "$OUTDIR/6_7_post_FST"
