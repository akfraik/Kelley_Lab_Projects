###########################################################################################################################
## 7a. XTX Input file generation
  # Bayenv 2 -FST outlier test analog
###########################################################################################################################

mkdir XTX

## First step, is to take our converted BED files (pre and post-disease)
INDIR="/Plink_Files"
OUTDIR="/XTX"
POPULATIONS="1 2 3 4 5 6 7"     # Replace numbers with your population names
SCRIPTDIR="/Scripts"            # Specify the directory you have been storing all of your scripts (including this one) in
KEEPDIR="/Plink_Input"

## For loop with pre and post-disease status animals
  for i in pre post;
  do

## Make directories for pre and post-disease status animals
  mkdir XTX/${i}
  POPFILE="/Structure_Analyses/Devils_Locations_${i}.tsv"

## Python code for transposing a column
transpose="""
import sys
if len(sys.argv) > 1:
    delimiter = sys.argv[1]
else:
    delimiter = None
NF = 0
data = []
with sys.stdin as handle:
    for line in handle:
        if delimiter is None:
            fields = line.strip().split()
        else:
            fields = line.strip().split(delimiter)
        if NF == 0:
            NF = len(fields)
        data.append(fields)
with sys.stdout as out:
    for i in range(NF):
        out.write('\t'.join(line[i] for line in data) + '\n')
"""

mkdir -p $SCRIPTDIR
if [[ $PBS_O_WORKDIR ]]; then

## Load plink
  source /usr/share/modules/init/bash
  module load plink/1.90_extrachr  
  set -u
  set -e
  set -o pipefail
  
## Loop over populations to get allele count data
   for pop in $POPULATIONS; 
   do
     keep_file="${KEEPDIR}/${pop}.${i}.keep"
     count_file="${KEEPDIR}/${pop}.${i}.frq.count"
     allele_file="${KEEPDIR}/${pop}.${i}.allele"
     popfile="${KEEPDIR}/${i}.popfile.tmp"
      p=$(echo $pop)

 ## Make a list of individuals from this population
    awk '{print $4 "\t" $3}' $POPFILE > $popfile
    grep -P "${pop}" $popfile | sed 's/'${pop}'/1/g' > $keep_file \
     || { echo "Listing individuals for ${pop} failed"; exit 1; }

  ## Use plink to get the allele counts. By setting --keep-allele-order
   # we make sure that the alleles are the same for every population.
     plink --file $inprefix --keep $keep_file --allow-extra-chr --freq counts --keep-allele-order --out $count_file \
      || { echo "Getting counts for ${pop} failed"; exit 1; }

  ## Extract the snps so we know the order of SNPs later
    awk '{print $2};' "$count_file" > "${OUTDIR}/snps.txt" \
     || { echo "Extracting snps for ${pop} failed"; exit 1; }

    ## Now, write a 3rd for loop to produce files format for allele counts for SNPs for population 1
     file_line_number=$(cat "${OUTDIR}/$count_file" | wc -l)
        for ((k=2; k<=file_line_number; k++));
         do
            sed -n ''$k'p' "${OUTDIR}/$count_file" | awk -F " " '{print $5}' >> "${OUTDIR}/$allele_file"
            sed -n ''$k'p' "${OUTDIR}/$count_file" | awk -F " " '{print $6}' >> "${OUTDIR}/$allele_file"
          ## Now we tell the inner loop to run this many times before the outer loop
           # End of 3rd level loop
         done
       
    ## Now we tell the inner loop to run this many times before the outer loop
     # End of 2nd level loop.
    done 
    
## Dump the count columns into a temporary file, with tabs separating the entries. The temporary file will need to be
 # transposed before using it with bayenv.
  paste -d"\t" "${OUTDIR}/1.${i}.allele" "${OUTDIR}/2.${i}.allele" "${OUTDIR}/3.${i}.allele" "${OUTDIR}/4.${i}.allele" 
   "${OUTDIR}/5.${i}.allele" "${OUTDIR}/6.${i}.allele" >> "${OUTDIR}/SNPs_file_${i}" \
  || { echo "extracting output for population failed"; exit 1; }
    
## This "done" complete the for loop for each disease status (pre and post, the initial for loop)        
 # End of final, outer loop
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
