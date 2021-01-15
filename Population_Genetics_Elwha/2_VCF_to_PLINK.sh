###########################################################################################################################
## 2. Filter and convert VCF File
# Filter the STACKS output file, subset to useable VCF files for population genetics analyses and convert to plink
###########################################################################################################################

####################################################################################################################################
# In my scripts, if there are breaks like this that means run the code in chunks the first time through to check on
# intermediate files or to explore filtering parameters
####################################################################################################################################

## First step, create one master VCF file
# Setup
module load bio/vcftools
module load bio/plink/1.90b6.6

## Set up directory paths to where files are stored
BASE="/INSERT PATH TO BASE DIRECTORY"
INDIR="/INSERT PATH TO STACKS SNP CALLS"
INFODIR="/INSERT PATH TO INFO"
OUTDIR="/Filtered VCF Files home"
MIN_GP=1

#### ZIP VCF file containing the SNP calls from STACKS
#gzip "${INDIR}/populations.snps.vcf"

##### Retain only biallelic sites - don't need to use remove indels from the STACKS file as they do it internally
vcftools --gzvcf "${INDIR}/populations.snps.vcf.gz" --min-alleles 2 --max-alleles 2 \
        --recode-INFO-all --recode --out "${INDIR}/Elwha_test_STACKS"

####################################################################################################################################
#### Explore how much missing data per individual you want to allow
# MIND = missing data per individual - ran the code below to produce this list of individuals
# Now filter based on MIND to create files for PCAs to explore missing data patterns and create filtered files for downstream analyses
for j in 6 7 8 9;
do

## Now let's try looping through like we did for plink 
## Where in PLINK this is the maximum missing data rate per SNP, in VCFtools this is 
## The maximum missing data ALLOWED - so the metrics should be the inverse
MAXMISSINGIND=0.${j}   # Maximum missing data rate per individual

## First convert the vcf to bed file
plink --vcf "${INDIR}/Elwha_tmp_STACKS.recode.vcf" --vcf-min-gp 1 \
       --biallelic-only strict --make-bed --out "${OUTDIR}/Elwha_bed" \
       --allow-extra-chr --noweb --const-fid 1  

## This is a temporary file to remove individuals with REALLY high missing data
plink --bfile "${OUTDIR}/Elwha_bed" --allow-extra-chr --recode tab \
       --mind $MAXMISSINGIND --out "${OUTDIR}/all_MIND_${j}" 

## Create file of individuals with MIND (missing data per individual) at ${k}0% 
## To use to remove the individuals in subsequent fitlering scripts
# First make a "nosex" file or file with all the individual IDs in them that can be used by PLINK
awk '{print $1"\t"$2}' "${OUTDIR}/all_MIND_${j}.ped" > "${OUTDIR}/all_MIND_${j}.nosex"
sed -e 's/.bwamemalignPE_10kbunmap.sorted//g' "${OUTDIR}/all_MIND_${j}.nosex" | sort -k2 > "${INFODIR}/All_${j}"
cp "${INFODIR}/Subset_38.Exclude" "${INFODIR}/Subset_38_copy.Exclude"

# Then make one that can be used by VCFtools
sed -e 's/.bwamemalignPE_10kbunmap.sorted//g' "${INFODIR}/All_2710.keep" | sort -k2 > "${INFODIR}/All_2710"
comm "${INFODIR}/All_2710" "${INFODIR}/All_${j}" -2 -3 > "${INFODIR}/Subset_${j}.exclude"
cat "${INFODIR}/Subset_${j}.exclude" >> "${INFODIR}/Subset_38_copy.Exclude"

# Create a final set of individuals to "exclude" based on varying levels of ${j} = MIND
awk -F'\t' '!seen[$2]++' "${INFODIR}/Subset_38_copy.Exclude" > "${INFODIR}/Subset_${j}.exclude"

done
####################################################################################################################################

#### Remove individuals with >60 % MIND + bad individuals from the data set
# MIND = missing data per individual - ran the code below to produce this list of individuals
## Using the output from exploring misisng data, using the commented code above, produced an output file "${INFODIR}/Exclude_List.exclude" 
## Excludes all of those individuals that had > 60% missing data per individual 
vcftools --gzvcf "${INDIR}/Elwha_test_STACKS.recode.vcf.gz" --remove "${INFODIR}/Exclude_List.exclude" \
	--recode-INFO-all --recode --out "${INDIR}/Elwha_tmp_STACKS"

####################################################################################################################################
#### Explore how much missing data per site you want to allow
# MISS = missing data per site (SNP) - ran the code below to produce this list of SNPS
# MSD = maximum missing sequencing depth per site
# Now filter based on MISS and to create files for PCAs to explore missing data patterns and create filtered files for downstream analyses

#### Now filter based on MSD to create files for PCAs
for k in 5 6 7 8;
do

#### Loop over files for MISS parameter
for i in 1 2 3 4;
do

#### Filtering parameters
MAXMEANDEPTH=${k}   # Minimum sequencing depth per SNP
MAXMISSINGSNP=0.${i}   # Maximum missing data allowed per SNP

#### First things first, filter VCF for minimum mean depth per site [$k]
vcftools --vcf "${INDIR}/Elwha_tmp_STACKS.recode.vcf" --recode --recode-INFO-all \
       --min-meanDP $MAXMEANDEPTH --out "${OUTDIR}/subset_dp_${k}"

#### Now zip up these files
gzip "${OUTDIR}/subset_dp_${k}.recode.vcf"

#### Now convert the vcf to a bed file
plink --vcf "${OUTDIR}/subset_dp_${k}.recode.vcf.gz" --vcf-min-gp 1 \
       --biallelic-only strict --make-bed --out "${OUTDIR}/Elwha_${k}_bed" \
       --allow-extra-chr --noweb --const-fid 1  

#### Now because the sample IDs and the variant IDs are different we are going to standardize them between VCF and PLINK
awk '{print $1"\t""test.test""\t"$3"\t"$4"\t"$5"\t"$6}' "${OUTDIR}/Elwha_${k}_bed.bim" | sed -e 's/test//g' > "${OUTDIR}/test.bim"
mv "${OUTDIR}/test.bim" "${OUTDIR}/Elwha_${k}_bed.bim"

#### Now convert the bed to a ped file 
# Set missing variant ids in a way that will be comparable between plink/vcftools
plink --bfile "${OUTDIR}/Elwha_${k}_bed"  --allow-extra-chr --recode tab --set-missing-var-ids @_# \
	--exclude "${INFODIR}/Paralogs_All" --out "${OUTDIR}/Elwha_${k}"

#### We need to fix all of the nomenclature for the names by removing ".bwamemalignPE_10kbunmap.sorted" 
# from the PED files and the nosex files
sed -e 's/.bwamemalignPE_10kbunmap.sorted//g' "${OUTDIR}/Elwha_${k}.ped" | sed -e 's/_CAT//g' > "${OUTDIR}/Elwha_${k}_test.ped"
awk '{print $1"\t"$2}' "${OUTDIR}/Elwha_${k}_test.ped" > "${OUTDIR}/Elwha_${k}_test.nosex"

#### Make the map and log files for PLINK
mv "${OUTDIR}/Elwha_${k}.map" "${OUTDIR}/Elwha_${k}_test.map"
mv "${OUTDIR}/Elwha_${k}.log" "${OUTDIR}/Elwha_${k}_test.log"

#### Now figure out what removing high levels of MISS looks like after removing MSD
plink --file "${OUTDIR}/Elwha_${k}_test" --allow-extra-chr --recode tab --maf 0.01 \
	--geno $MAXMISSINGSNP --out "${OUTDIR}/all_MISS_${k}${i}"

#### Create popfile for information regarding individuals
awk '{print $1"\t"$2}' "${OUTDIR}/all_MISS_${k}${i}.ped" > "${OUTDIR}/all_MISS_${k}${i}.nosex"

#### Close both for loops
done
done

####################################################################################################################################

#### After exploring for missing data by MIND, MISS, MSD, then filter and create a file you like
# Remove sites with sequencing depth < 8 and > 10 % missing data per site (therefore 81)
# Use the output from the previous steps that contains the combination of filtering parameters that you prefer
# Then convert the PLINK files to proper BED files for downstream population structure analyses
plink --noweb --file "${OUTDIR}/all_MISS_81" --make-bed --allow-extra-chr --keep "${INFODIR}/All.Keep" \
  --out "${OUTDIR}/all_81_bed"

#### For loop to create special bed files for PCA
for ext in ".bed" ".fam" ".nosex" ".bim"; do
   cp "${OUTDIR}/all_81_bed$ext" \
       "${OUTDIR}/all_81_bed_ChrNum$ext"
done

#### Create the proper BED files for Admixture/Fast Structure
sed 's/.\+_GL//g' "${OUTDIR}/all_81_bed.bim" | \
   sed 's/_random\t/\t/g' > \
   "${OUTDIR}/all_81_bed_ChrNum.bim"

#### Print a line to remind yourself what your selected filtering parameters were
echo "removed individuals with > 60% missing data per individual and sites with sequencing depth < 8 and > 10 % missing data per site"

#### Convert BED to VCF to make a file with ALL the samples you want to keep
plink --noweb --bfile "${OUTDIR}/all_81_bed" --recode vcf --allow-extra-chr --keep "${INFODIR}/All.Keep" \
	--out "${OUTDIR}/all_81_vcf"

#### Then rezip your vcf file
gzip "${OUTDIR}/all_81_vcf.vcf"

#### Zip up original STACKS files after use
gzip "${INDIR}/Elwha_tmp_STACKS.recode.vcf"
gzip "${INDIR}/Elwha_test_STACKS.recode.vcf"
