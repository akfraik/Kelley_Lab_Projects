###########################################################################################################################
## 7b. Genetic Stock Identification
# Create the input files for RUBIAS
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

#### Loop through SNPs to retain for Genetic Stock Identification Analyses
for SNP in GSI GSI_1;
do

#### Change outdirectory for output
OUTDIR=${BASE}/${SNP}

#### Take the output from R and make it useable here
sed -e 's/"//g' "${INFODIR}/SNPs_Location_${SNP}.csv" | sed -e 's/,/\t/g' | awk '{print $1}' | sed 1d > "${INFODIR}/${SNP}_Set"

#### Make a file that has less missing data
plink --noweb --bfile "${INDIR}/all_81_bed" --make-bed --allow-extra-chr --extract "${INFODIR}/${SNP}_Set" --out "${OUTDIR}/GSI"

##### Create the Pre/Post Dam removal Keep files
cut -d',' -f27 "${INFODIR}/pop_code_subset.csv" | sort | uniq > "${INFODIR}/All_Times.txt"

#### Use pre and post-dam removal keep files for generating GSI files
for time in Pre Post;
do 

## Now we are going to make a file for the csv file for GSI
## Now convert to the right kind of file
plink --noweb --bfile "${OUTDIR}/GSI" --allele1234 --allow-extra-chr --recode tab --keep "${INFODIR}/${time}.Keep" --out "${OUTDIR}/GSI_${time}"

#### Make "keep" file 
awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1]' "${INFODIR}/all_81.keep" $METADATA > "${INFODIR}/Elwha_All"
awk -F $'\t' -v i="$time" '$21 == i' "${INFODIR}/Elwha_All" | sed -e 's/^M//g' > "${INFODIR}/Elwha_${time}"

##### Now we want to edit these files so that we paste the second column with population info
# Replace column of 1s with site.x
awk 'FNR==NR{a[NR]=$2;next}{$1=a[FNR]}1' "${INFODIR}/Elwha_${time}" "${OUTDIR}/GSI_${time}.ped" > "${OUTDIR}/test.ped"
awk 'FNR==NR{a[NR]=$2;next}{$1=a[FNR]}1' "${INFODIR}/Elwha_${time}" "${OUTDIR}/GSI_${time}.nosex" > "${OUTDIR}/test.nosex"
mv "${OUTDIR}/test.ped" "${OUTDIR}/GSI_${time}.ped"
mv "${OUTDIR}/test.nosex" "${OUTDIR}/GSI_${time}.nosex"

#### Replace spaces with commas
sed -e 's/ /,/g' "${OUTDIR}/GSI_${time}.ped" > "${OUTDIR}/GSI_${time}.csv"

done

#### Now we rename everything appropriately
mv "${OUTDIR}/GSI_Pre.csv" "${OUTDIR}/GSI_Ref.csv"
mv "${OUTDIR}/GSI_Post.csv" "${OUTDIR}/GSI_Mix.csv"

done
