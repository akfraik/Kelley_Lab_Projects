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

outer=1

## For loop with pre and post-disease status animals
  for i in pre post;
  do

## Make directories for pre and post-disease status animals
  mkdir XTX/${i}
  POPFILE="/Structure_Analyses/Devils_Locations_${i}.tsv"

## Now pass the command for the beginning of the outer loop
  echo "Pass $outer in outer loop."
  echo "---------------------"
  inner=1           # Reset inner loop counter.

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
  set -o pipefail
  
## Loop over populations to get allele count data
   for pop in $POPULATIONS; 
   do
     keep_file="${KEEPDIR}/${pop}.${i}"
     count_file="${KEEPDIR}/${pop}.{i}"
     allele_file="${KEEPDIR}/${pop}.allele.${i}"
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
    awk '{print $2};' "${count_file}.frq.count" > "${OUTDIR}/snps.txt" \
     || { echo "Extracting snps for ${pop} failed"; exit 1; }

    ## Now, write a 3rd for loop to produce files format for allele counts for SNPs for population 1
     file_line_number=$(cat "${OUTDIR}/Freycinet.frq.count" | wc -l)
        for ((k=2; k<=file_line_number; k++));
         do
            sed -n ''$k'p' "${OUTDIR}/${pop}.frq.count" | awk -F " " '{print $5}' >> "${OUTDIR}/${pop}.allele.before"
            sed -n ''$k'p' "${OUTDIR}/${pop}.frq.count" | awk -F " " '{print $6}' >> "${OUTDIR}/${pop}.allele.before"
          ## Now we tell the inner loop to run this many times before the outer loop
            echo "Pass $inner_inner in inner_inner loop."
            let "inner_inner+=1"  # Increment inner, inner loop counter, or the 3rd level of loops
          # End of 3rd level loop
         done
       
    ## Now we tell the inner loop to run this many times before the outer loop
      echo "Pass $inner in inner loop."
      let "inner+=1"  # Increment inner loop counter.
     # End of 2nd level loop.
    done 
    
## This "done" complete the for loop for each disease status (pre and post, the initial for loop)        
 # End of outer loop.
 let "outer+=1"    # Increment outer loop counter. 
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
