###########################################################################################################################
# 1. Convert VCF to PED file
# 1a. Convert the Beagle imputed VCF output file a to Plink file
###########################################################################################################################

## First step, create one master VCF file

# Minimum genotype probability to allow (otherwise set genotype to missing)
  MIN_GP=1

# Directories where files are stored
INDIR="/PATH_TO_INPUT_FILES"
OUTDIR="/Plink_Files"
# Replace numbers with your population names
POPULATIONS="1 2 3 4 5 6 7"
Keep_File="$OUTDIR/focal_individuals.keep"

# A little trick in SLURM so that you create output and input
# Load plink
  if [[ $PBS_O_WORKDIR ]]; then=
  
  source /usr/share/modules/init/bash
  module load bcftools/2016-04-08
  module load plink/1.90b3
    
  set -u
  set -e
  set -o pipefail
   
## bcftools to merge population files into one VCF
bcftools merge --output-type "z" --output "${OUTDIR}/all.gp.vcf.gz" \
        "${INDIR}/1.gp.vcf.gz" "${INDIR}/2.gp.vcf.gz" "${INDIR}/3.gp.vcf.gz" "${INDIR}/4.gp.vcf.gz" \
        "${INDIR}/5.gp.vcf.gz" "${INDIR}/6.gp.vcf.gz" "${INDIR}/7.gp.vcf.gz" \
    || { echo "merging failed"; exit 1; }
    
## Then convert the master, merged VCF into a PED files
# Filtering and quality control steps provided in subsequent scripts for filtering out variants
   plink --vcf "${OUTDIR}/all.gp.vcf.gz" --vcf-min-gp $MIN_GP --allow-extra-chr \
        --const-fid 1 --recode --out "${OUTDIR}/ped_all" --keep $Keep_File \
    || { echo "creating plink file failed"; exit 1; }

## Now create filtered plink files for pre-disease and post-diseased populations
# SNP filtering parameters
MAXMISSINGIND=0.8   # Maximum missing data rate per ind.; applied before the other filters
MINMAF=0.05         # Minimum minor allele frequency across all individuals
MAXMISSINGSNP=0.6   # Maximum missing data rate per SNP across all individuals

## Looking for lines that don't have the pattern and we want to remove all SNPs just on the X chromosome
    grep -v "chrX" $mapfile | cut -f 2 > "${INDIR}/noX.txt" \
        || { echo "Combining SNPs filter files failed"; exit 1; }

#Creating "pre-disease" and "post-disease" files
     plink --allow-extra-chr --file "${OUTDIR}/ped_all" --extract "${INDIR}/noX.txt" --mind $MAXMISSINGIND \
     --maf $MINMAF --geno $MAXMISSINGSNP --recode --out "${OUTDIR}/ped_tmp" --keep $keep_file \
        || { echo "plink failed"; exit 1; }

 ## The --keep and --extract parameters in this command are optional for generating a pre-disease PED file for analyses
 # --keep can be used for pulling out specific indivduals (in this case, the pre-disease individuals)
 # --extract can be used for pulling out specific SNPs
      plink --allow-extra-chr --file "${OUTDIR}/ped_tmp" --recode \
      --out "${OUTDIR}/ped_pre_disease" --keep "${OUTDIR}/focal_pre_disease.keep" \
        || { echo "final pre-disease plink failed"; exit 1; }

 ## Post-disease PED file for analyses
      plink --allow-extra-chr --file "${OUTDIR}/ped_tmp" --recode \
      --out "${OUTDIR}/ped_post_disease" --keep "${OUTDIR}/focal_post_disease.keep" \
        || { echo "final pre-disease plink failed"; exit 1; }

  rm "${OUTDIR}/*tmp"

else

# SLURM equivalent to qsub to create copies of the scripts and create output and error files
# Makes a copy of the script every time it runs. The only requirement is that you don't try to
# run this script with source; relies on the fact that the name of this script will be contained in $0 variable, 
# so it can just copy itself to a new directory. The "basename" command gets rid of the directory name part of a file. 
# The nice thing about this is that it leaves an unambiguous record of what you did every time the script runs.
    sfile="${INDIR}/script_copies/$(date '+%Y-%m-%d')-$(basename $0)"
    cp $0 $sfile
    cd "${INDIR}/working"
    qsub $sfile
    echo "${INDIR}/working"
  
