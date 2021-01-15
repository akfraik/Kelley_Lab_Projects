###########################################################################################################################
## 4a. Population structure analyses
# Produce input files for fast Structure runs
###########################################################################################################################

#### Setup
module load bio/vcftools
module load bio/plink/1.90b6.6
module load bio/faststructure/1.0 

#### Set up directory paths to where files are stored
BASE="/INSERT PATH TO BASE DIRECTORY"
INDIR="/INSERT PATH TO VCF FILES"
INFODIR="/INSERT PATH TO INFO"
OUTDIR="/INSERT PATH TO OUTPUT FILE LOCATION"
METADATA="${INFODIR}/INSERT PATH TO METADATA"

#### Loop over every unique time period during dam removal: pre/post
# then make "keep" files for individuals from each time period
while read time;
do
echo $time
awk -F $'\t' -v i="$time" '$21 == i' $METADATA | sed -e 's/,/\t/g' | awk '{print $1}' | awk '{print 1"\t"$1}' > "${INFODIR}/${time}.Keep"

done <"${INFODIR}/All_Times.txt"

#### Loop through pre and post dam removal individuals
for i in Pre Post;
do

#### We want to only look at all individuals first prior to dam removal
# then post-dam removal to evaluate if/how population structure has changed
plink --bfile "${INDIR}/all_81_bed" --make-bed --keep "${INFODIR}/${i}.Keep" \
	--out "${OUTDIR}/Pre_bed" --allow-extra-chr --noweb

#### For loop to create special bed files for FastSTRUCTURE
for ext in ".bed" ".fam" ".nosex" ".bim"; do
   cp "${OUTDIR}/${i}_bed$ext" \
       "${OUTDIR}/${i}_bed_ChrNum$ext"
done

#### Still formating bed files for FastSTRUCTURE
sed 's/.\+_GL//g' "${OUTDIR}/${i}_bed.bim" | \
   sed 's/_random\t/\t/g' > \
   "${OUTDIR}/${i}_bed_ChrNum.bim"

#### Iterate through number of "K" values or population sizes you want to test
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
do 

#### Create fast structure output
structure.py -K $K --input "${OUTDIR}/${i}_bed" --output "${OUTDIR}/${i}" --prior simple --full --seed 100


#### Finish looping over possible K values
done

#### Choose K
chooseK.py --input "${OUTDIR}/${i}"

#### Finish looping over time
done
