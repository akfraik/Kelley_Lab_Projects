###########################################################################################################################
## 6a. Population genetics analyses
# Calculate Pi, Ajk, Tajima's D
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

################################################################################################
#### Start by looping over the metadata file to creae a keep file for each time period
# Loop over every unique time period during dam removal: pre/post
################################################################################################
cut -d',' -f21 $METADATA | sort | uniq > "${INFODIR}/All_Times.txt"

#### Loop over time periods
while read time;
do
echo $time
awk -F $'\t' -v i="$time" '$21 == i' $METADATA | sed -e 's/,/\t/g' | awk '{print $1}' | awk '{print 1"\t"$1}' > "${INFODIR}/${time}.Keep"

#### Loop over time periods to convert PLINK to VCF files
plink --noweb --bfile "${INDIR}/all_81_bed" --recode vcf --keep "${INFODIR}/${time}.Keep" \
        --allow-extra-chr --out "${INDIR}/Elw_${time}"

#### Zip up the VCF files per dam removal time period
gzip "${INDIR}/Elw_${time}.vcf"

#### Calculate Tajima's D for Pre and Post-dam Removal
vcftools --gzvcf "${INDIR}/Elw_${time}.vcf.gz" --TajimaD 10000 --out "${OUTDIR}/${time}_D"

#### Calculate Relatedness (Ajk) among individuals Pre and Post-dam Removal
vcftools --gzvcf "${INDIR}/Elw_${time}.vcf.gz" --relatedness --out "${OUTDIR}/${time}"

#### Calculate Nucleotide Diversity (Pi) Pre and Post-dam Removal
vcftools --gzvcf "${INDIR}/Elw_${time}.vcf.gz" --site-pi --out "${OUTDIR}/${time}_Pi"

done < "${INFODIR}/All_Times.txt"

################################################################################################
#### Start by looping over the metadata file to creae a keep file for each life-history form
# Loop over Steelhead and Unknown fish
################################################################################################
cut -d',' -f20 $METADATA | sort | uniq | sed -e 's/\t/,/g' > "${INFODIR}/All_Species.txt"

#### Loop over every unique life-history form
while read species;
do
echo $species
awk -F $'\t' -v i="$species" '$20 == i' $METADATA | sed -e 's/,/\t/g' | awk '{print $1}' | sed -e 's/^/1_/g' > "${INFODIR}/${species}.Keep"

#### Loop through time periods and subset by species
vcftools --gzvcf "${INDIR}/Elw_Pre.vcf.gz" --keep "${INFODIR}/${species}.Keep" --recode --recode-INFO-all --out "${INDIR}/Elw_Pre_${species}"

#### Zip up the VCF files per dam removal time period
gzip "${INDIR}/Elw_Pre_${species}.recode.vcf"

#### Calculate Tajima's D among life-history forms
vcftools --gzvcf "${INDIR}/Elw_Pre_${species}.recode.vcf.gz" --TajimaD 10000 --out "${OUTDIR}/Elw_Pre_${species}_D"

#### Calculate Relatedness (Ajk) among life-history forms
vcftools --gzvcf "${INDIR}/Elw_Pre_${species}.recode.vcf.gz" --relatedness "${OUTDIR}/Elw_Pre_${species}"

#### Calculate Nucleotide Diversity (Pi) among life-history forms
vcftools --gzvcf "${INDIR}/Elw_Pre_${species}.recode.vcf.gz" --site-pi --out "${OUTDIR}/Elw_Pre_${species}_Pi"

done < "${INFODIR}/All_Species.txt"

################################################################################################
#### Start by looping over the metadata file to creae a keep file for each population
# Loop over AD/ID/BD/SBLR
################################################################################################
cut -d',' -f22 $METADATA | sort | uniq | sed -e 's/\t/,/g' > ${INFODIR}/All_Locations.txt"

#### Loop over every unique population ID
while read locations;
do
echo $locations
awk -F $'\t' -v i="$locations" '$22 == i' $METADATA | sed -e 's/,/\t/g' | awk '{print $1}' | sed -e 's/^/1_/g' > "${INFODIR}/${locations}.Keep"

#### Loop through time periods and subset by AD/BD/ID 
vcftools --gzvcf "${INDIR}/Elw_Pre.vcf.gz" --keep "${INFODIR}/${locations}.Keep" --recode --recode-INFO-all --out "${INDIR}/Elw_Pre_${locations}"

#### Zip up the VCF files AD/ID/BD
gzip "${INDIR}/Elw_Pre_${locations}.recode.vcf"

#### Calculate Tajima's D among life-history forms
vcftools --gzvcf "${INDIR}/Elw_Pre_${locations}.recode.vcf.gz" --TajimaD 10000 --out "${OUTDIR}/Pre_${locations}_D"

#### Calculate Relatedness (Ajk) among life-history forms
vcftools --gzvcf "${INDIR}/Elw_Pre_${locations}.recode.vcf.gz" --relatedness "${OUTDIR}/Pre_${locations}"

#### Calculate Nucleotide Diversity (Pi) among life-history forms
vcftools --gzvcf "${INDIR}/Elw_Pre_${locations}.recode.vcf.gz" --site-pi --out "${OUTDIR}/Pre_${locations}_Pi"

done < "${INFODIR}/All_Locations.txt"

################################################################################################
#### Start by looping over the metadata file to creae a keep file for each life-history cohort
# Loop over Juveniles/Adults for each natal year
################################################################################################
cut -d',' -f26 $METADATA | sort | uniq | sed -e 's/\t/,/g' > ${INFODIR}/All_Cohorts.txt"

#### Loop over every unique cohort
while read cohorts;
do
echo $cohorts
awk -F $'\t' -v i="$cohorts" '$26 == i' $METADATA | sed -e 's/,/\t/g' | awk '{print $1}' | sed -e 's/^/1_/g' > "${INFODIR}/Post_${cohorts}.Keep"

#### Loop through time periods and subset by AD/BD/ID 
vcftools --gzvcf "${INDIR}/Elw_Post.vcf.gz" --keep "${INFODIR}/Post_${cohorts}.Keep" --recode --recode-INFO-all --out "${INDIR}/Elw_Post_${cohorts}"

#### Zip up the VCF files Post_Cohorts
gzip "${INDIR}/Elw_Post_${cohorts}.recode.vcf"

#### Calculate Tajima's D among life-history cohorts
vcftools --gzvcf "${INDIR}/Elw_Post_${cohorts}.recode.vcf.gz" --TajimaD 10000 --out "${OUTDIR}/Post_${cohorts}_D"

#### Calculate Relatedness among life-history cohorts
vcftools --gzvcf "${INDIR}/Elw_Post_${cohorts}.recode.vcf.gz" --relatedness "${OUTDIR}/Post_${cohorts}"

#### Calculate Nucleotide Diversity (Pi) among life-history cohorts
vcftools --gzvcf "${INDIR}/Elw_Post_${cohorts}.recode.vcf.gz" --site-pi --out "${OUTDIR}/Post_${cohorts}_Pi"

done < "${INFODIR}/All_Cohorts.txt"
