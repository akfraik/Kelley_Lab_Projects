###########################################################################################################################
## 4b. Run XTX
  # Bayenv 2 -FST outlier test analog
###########################################################################################################################

## First step, is to take our converted BED files (pre and post-disease)
INDIR="/XTX"
POPULATIONS="1 2 3 4 5 6 7"     # Replace numbers with your population names
SCRIPTDIR="/Scripts"            # Specify the directory you have been storing all of your scripts (including this one)
KEEPDIR="/Plink_Input"

## For loop with pre and post-disease status animals
  for i in pre post;
  do

## Now set the output directory
  OUTDIR="/XTX/${i}"
  cd $OUTDIR

## Load plink
  source /usr/share/modules/init/bash
  module load plink/1.90_extrachr  
  set -u
  set -e
  set -o pipefail

## Before even iterating through SNPs to calculate the XTX or Bayes Factors, we need to generate a matrix
 module load bayenv2
 bayenv2 -i "{INDIR}/SNPs_file_${i}.txt" -p 6 -k 100000 > "${INDIR}/Cov_Matrix_${i}.txt"
|| { echo "failed on matrix generation; exit 1; }
 
 ## Number of SNPs = NSNPS=(10184), i.e., lines (alleles) in SNPs file / 2
    for s in $(seq 2 2 $((10184 * 2)))
    do
 
    ## Create temporary files in which the first two lines (frequency of allele 1 and 2 for snp -n respectively)
      # Then delete the temporary files
      head -n $s "${INDIR}/SNPs_file_${i}" | tail -n 2 > "${OUTDIR}/tmp_${i}.snp"

     ## Then run the actual Bayenv2 command including the "X" flag to calculate the XTX value for each SNP across pops
      # DUM_ENVIRON = "dummy" environmental file that contains dummy variables for each population (p = 6 in my case)
      bayenv2 -i "${OUTDIR}/tmp_${i}.snp" -m "${INDIR}/Cov_Matrix_${i}.txt" -e "${INDIR}/DUM_ENVIRON" -p 6 -k 100000 -n 2 -t -c -X -f \
       || { echo "failed on snp ${s}"; exit 1; }
    done
    
rm "${OUTDIR}/tmp_${i}.snp"

## This done is indicative of the pre and post-DFTD status for loop
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
