###########################################################################################################################
## 6c. Population genetics analyses
# Create the input files for Adegenet
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

#### Make files ready for adegenet keep file for Pre and Post-dam removal
for time in Pre Post;
do

#### Make Elwha files
awk -F $'\t' -v i="$time" '$21 == i' $METADATA | sed -e 's/^M//g' > "${INFODIR}/Elwha_${time}"

#### Convert PED to allele files for Adgenet 
plink --noweb --bfile "${INDIR}/all_81_bed" --allele1234 --allow-extra-chr --recode tab --keep "${INFODIR}/${time}.Keep" --out "${OUTDIR}/Adegenet_${time}"

##### Now we want to edit these files so that we paste the second column with population info
# Replace column 1 of ped and nosex files with column 2 from the METADATA file
# Replace 2nd column in f1 with 3rd column from f2
awk 'FNR==NR{a[NR]=$22;next}{$1=a[FNR]}1' "${INFODIR}/Elwha_${time}" "${OUTDIR}/Adegenet_${time}.ped" > "${OUTDIR}/${time}.ped"
awk 'FNR==NR{a[NR]=$22;next}{$1=a[FNR]}1' "${INFODIR}/Elwha_${time}" "${OUTDIR}/Adegenet_${time}.nosex" > "${OUTDIR}/${time}.nosex"
mv "${OUTDIR}/${time}.ped" "${OUTDIR}/Adegenet_${time}.ped"
mv "${OUTDIR}/${time}.nosex" "${OUTDIR}/Adegenet_${time}.nosex"

#### Replace spaces with commas
sed -e 's/ /,/g' "${OUTDIR}/Adegenet_${time}.ped" > "${OUTDIR}/Adegenet_${time}.csv"

#### Finish loop
done
