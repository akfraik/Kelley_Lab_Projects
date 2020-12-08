###########################################################################################################################
# 1. GATK SNP calling pipeline
# 1a. BAM clean up and renaming RG
###########################################################################################################################

### First change directories to the directory with your mapped files that have been converted to BAMs
cd $MAPDIR

## Start with output from the output from the script 2_Map_and_Trim.sh in the Devil_Trancriptome_Uninfected directory
module load java
# $LANE = lane information
# $ID = Illumina library prep ID
# $LIBRARY = Library prep number
# ${i} = sample number or sample ID depending on how you coded it in previous script (I coded it as sample number)

## $Sample_List= List of samples (space separated) in order of the samples you're feeding into program
# we create an array prior to the for loop to generate matching of sample number to sample ID for read group info
arr=($Sample_List)

for i in 1..20;

do

RG_ID="echo ${arr[${i}]}"

java -jar picard.jar AddOrReplaceReadGroups I=$MAPDIR/${i}.bam O=$MAPDIR/${i}_rg_added_sorted.bam SO=coordinate \
RGID=$RG_ID RGLB=$LIBRARY RGPL=Illumina RGPU=$LANE RGSM=${i}

java -jar picard.jar MarkDuplicates I=$MAPDIR/${i}_rg_added_sorted.bam O=$MAPDIR/${i}_dedupped.bam \
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$MAPDIR/${i}_output.metric

# Remove temporary files and zip up the others to reduce space
rm $MAPDIR/${i}_rg_added_sorted.bam
bgzip $MAPDIR/${i}.bam

done

###########################################################################################################################
# 1b. Split and Trim
###########################################################################################################################

# $REFDIR = reference annotation directory
# ref.fa = reference file in the reference directory
module load java

# Samples to run
for i in 1..20;

do

# Command for split N Trim right from the Broad Institute
java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R $REFDIR/ref.fa -I $MAPDIR/${i}_dedupped.bam -o $MAPDIR/${i}_split.bam \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

# Echo the sample that completes 
echo ${i}_split.bam

# Remove temporary files and zip up the others to reduce space
rm $MAPDIR/${i}_dedupped.bam

done

###########################################################################################################################
# 1c. Base Quality Recalibration Score
###########################################################################################################################

# $REFDIR = reference annotation directory
module load java

# Samples to run
for i in 1..20;

do

# Command for BQRS
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REFDIR/ref.fa -I $MAPDIR/${i}_split.bam 
-o $MAPDIR/${i}_recal_data.table --run_without_dbsnp_potentially_ruining_quality

# Recreate recalibrate BAM files with improved scores
java -jar GenomeAnalysisTK.jar -T PrintReads -R $REFDIR/ref.fa -I $INDIR/${i}_split.bam 
-BQSR $MAPDIR/${i}_recal_data.table -o $MAPDIR/${i}_recal.bam

#Echo the sample that completes 
echo $MAPDIR/${i}_recal.bam

samtools index $MAPDIR/${i}_recal.bam

done

###########################################################################################################################
# 1d. Variant Calling
###########################################################################################################################

# $REFDIR = reference annotation directory
module load java

#Create a directory or designate a directory for SNPs to go in
# Samples to run
for i in 1..20;

do

#Sample 1-1 was run as the test
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $REFDIR/ref.fa -I $INDIR/${i}_recal.bam -dontUseSoftClippedBases \
-stand_call_conf 20.0 -o $SNPDIR/${i}_SNPs.vcf

#Echo the sample that completes 
echo ${i}_SNPs.vcf

gzip $SNPDIR/${i}_SNPs.vcf

done

###########################################################################################################################
# 1e. Variant Filtering
###########################################################################################################################

# $REFDIR = reference annotation directory
module load java

#Create a directory or designate a directory for SNPs to go in
# Samples to run
for i in 1..20;

do

# Unzip prior to filtering
gunzip $SNPDIR/${i}_SNPs.vcf

# Filtering commands
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R $REFDIR/ref.fa -V $SNPDIR/${i}_SNPs.vcf -window 35 -cluster 3 \
-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o $SNPDIR/${i}_SNPs_filtered.vcf

# Echo the sample that completes 
echo ${i}_SNPs_filtered.vcf

# Rezip files
gzip $INDIR/*.vcf.gz

done
