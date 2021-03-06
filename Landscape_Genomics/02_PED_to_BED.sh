###########################################################################################################################
# 2. Converting PED to BED files
# Input required for Fast Structure and Admixture 
###########################################################################################################################

## First step, is to convert our PED files to BED files removing SNPs that are putatively in strong LD
# So we are going to go back to an intermediate file in the previous script

INDIR="/Plink_Files""
OUTDIR="/Structure_Analyses"
# Replace numbers with your population names
POPULATIONS="1 2 3 4 5 6 7"
Keep_File="$OUTDIR/focal_individuals.keep"

# SNP filtering parameters
MAXMISSINGIND=0.8   # Maximum missing data rate per ind.; applied before the other filters
MINMAF=0.05         # Minimum minor allele frequency across all individuals
MAXMISSINGSNP=0.1   # Maximum missing data rate per SNP across all individuals

# A little trick in SLURM so that you create output and input
  if [[ $PBS_O_WORKDIR ]]; then

# Load plink
   source /usr/share/modules/init/bash
   module load plink/1.90_extrachr

   set -e
   set -u
   set -o pipefail
   
 ## Creating pre-disease PED files with more stringent filtering for structural analyses
   plink --allow-extra-chr --file "$INDIR/ped_all.ped" --mind $MAXMISSINGIND --maf $MINMAF --geno $MAXMISSINGSNP \
        --recode --out "${OUTDIR}/ped_pre_dftd_fast_struc_tmp" --keep $keep_file \
        || { echo "ptmp pre-disease plink failed"; exit 1; }

        plink --allow-extra-chr --file "${OUTDIR}/ped_pre_dftd_fast_struc_tmp" \
        --recode --out "${OUTDIR}/ped_pre_dftd_fast_struc" --keep $keep_file \
        || { echo "final plink failed"; exit 1; }

 ## Creating post-disease PED files with more stringent filtering for structural analyses
   plink --allow-extra-chr --file "$INDIR/ped_all.ped" --mind $MAXMISSINGIND --maf $MINMAF --geno $MAXMISSINGSNP \
        --recode --out "${OUTDIR}/ped_post_dftd_fast_struc_tmp" --keep $keep_file \
        || { echo "tmp post-disease plink failed"; exit 1; }

        plink --allow-extra-chr --file "${OUTDIR}/ped_post_dftd_fast_struc_tmp" \
        --recode --out "${OUTDIR}/ped_post_dftd_fast_struc" --keep $keep_file \
        || { echo "final plink failed"; exit 1; }
  
## Creating BED files for pre-disease
 # Use all the above to create regular bed file pre-disease
  plink --noweb --file "${OUTDIR}/ped_pre_dftd_fast_struc" --make-bed --allow-extra-chr --out "${OUTDIR}/bed_pre_dftd_fast_struc"

## For loop to create additional requried files
   for ext in ".bed" ".fam" ".nosex" ".bim"; do
    cp "${OUTDIR}/bed_pre_dftd_fast_struc$ext" "${OUTDIR}/bed_pre_dftd_fast_struc_ChrNum$ext"
   done
   
   sed 's/.\+_GL//g' "${OUTDIR}/bed_pre_dftd_fast_struc.bim" | sed 's/_random\t/\t/g' > "${OUTDIR}/bed_pre_dftd_fast_struc_ChrNum.bim"

## Creating BED files for post-disease
 # Use all the above to create regular bed file post-disease
  plink --noweb --file "${OUTDIR}/ped_post_dftd_fast_struc" --make-bed --allow-extra-chr --out "${OUTDIR}/bed_post_dftd_fast_struc"

## For loop to create additional requried files
   for ext in ".bed" ".fam" ".nosex" ".bim"; do
    cp "${OUTDIR}/bed_post_dftd_fast_struc$ext" "${OUTDIR}/bed_post_dftd_fast_struc_ChrNum$ext"
   done
   
   sed 's/.\+_GL//g' "${OUTDIR}/bed_post_dftd_fast_struc.bim" | sed 's/_random\t/\t/g' > "${OUTDIR}/bed_post_dftd_fast_struc_ChrNum.bim"

  rm "${OUTDIR}/*tmp"  
  
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
