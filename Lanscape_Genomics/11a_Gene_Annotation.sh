###########################################################################################################################
## 11a. Find the gene annotations by using Bedtools
###########################################################################################################################

## First step, is to take our converted PED (pre and post-disease) and convert them to a BED file
INDIR="${PROJDIR}/Plink_Input"
OUTDIR="${PROJDIR}/Bed_Files"

module load vcftools

if [[ $PBS_O_WORKDIR ]]; then

## Load plink
   source /usr/share/modules/init/bash
   module load plink/1.90_extrachr
   set -e
   set -u
   set -o pipefail

## For loop with pre and post-disease status animals
  for i in pre post;
  do
  
## Location of input directory and SNP file
  mkdir /Minotaur/${i}
  SNPs_List="/Minotaur/${i}/SNPs_List_${i}"
  
## Creating Plink files with just the desired SNPs
  plink --allow-extra-chr --file "${INDIR}/beagle_rap_pre" --recode --out "${OUTDIR}/outliers_${i}_disease" --extract $SNPs_List \
   || { echo "final plink failed"; exit 1; }
  awk '{print $1 "\t" $4-1 "\t" $4 "\t" $1 "-" $4};' "${OUTDIR}/outliers_${i}_disease.map" > "${OUTDIR}/outliers_${i}.bed"

## Look for overlap using the "Window" option
 # bedtools window [OPTIONS] [-a|-abam] -b <BED/GFF/VCF>
  bedtools window -a "${OUTDIR}/outliers_${i}.bed" -b "${OUTDIR}/biomart_reduced.bed" > "${OUTDIR}/Annotated_Outliers_Bed_${i}"

## VCF tools to create pre-DFTD VCF files for SNP2GO from bedtools overlap tool
   vcftools --vcf "${INDIR}/ped_all" --keep "${INDIR}/focal_${i}_individuals.keep" --recode --recode-INFO-all 
    --snps "${OUTDIR}/Annotated_Outliers_Bed_${i}" --out "${OUTDIR}/Annotated_Outliers_Bed_${i}.vcf"

## Make scripts for SNP2GO
   mkdir SNP2GO
   mkdir SNP2GO/${i}
   sed -e 's/Variable/${i}/g' "$SCRIPTDIR/11b_SNP2GO_${i}.R" > "$SCRIPTDIR/11b_SNP2GO_${i}.R"

## now give permission to run this script within this for loop
   chmod u+x "$SCRIPTDIR/11b_SNP2GO_${i}.R" 
   "./$SCRIPTDIR/11b_SNP2GO_${i}.R"

done

else

# SLURM equivalent to qsub to create copies of the scripts and create output and error files
# Makes a copy of the script every time it runs. The only requirement is that you don't try to
# run this script with source; relies on the fact that the name of this script will be contained in $0 variable, 
# so it can just copy itself to a new directory. The "basename" command gets rid of the directory name part of a file. 
# The nice thing about this is that it leaves an unambiguous record of what you did every time the script runs.
    d=`date "+%Y%m%d-%H%M"`
    sfile="${SCRIPTDIR}/$(basename $0)-${d}.sh"
    rm -f $sfile
    cp $0 $sfile
    qsub $sfile
 fi
