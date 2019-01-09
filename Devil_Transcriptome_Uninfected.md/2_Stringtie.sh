###########################################################################################################################
# 2. Stringtie to generate gene and transcript count matrices
# 2a. Convert SAM to BAM files and order
###########################################################################################################################

## We use samtools to create sorted and indexed BAM files for our samples
for i in 1 .. 20;

do

# Convert SAM to BAM files
samtools view -bS $MAPDIR/${i}.sam > $MAPDIR/${i}.bam

# Sort our BAM files
samtools sort -@ 8 $MAPDIR/${i}.bam -o $MAPDIR/${i}_sorted.bam 

# Index our sorted BAM files
samtools index $MAPDIR/${i}_sorted.bam 

# "echo" or print the information for the sample that completes then "rm" or delete the temporary files
echo ${i}.bam
rm $MAPDIR/${i}.bam 
rm $MAPDIR/${i}.sam

done

###########################################################################################################################
# 2b. Stringtie Alignment
###########################################################################################################################

## We used Stringtie to assemble and quantify expressed genes and transcripts 
# BAM files -> gtf files
# $REFANN = reference annotation (gtf file)
# $LANE = lane number for samples
# $COUNTDIR = output directory of your choosing for count files
for i in 1 .. 20;

do

stringtie -p 4 -e -B -G $REF_ANN -o $COUNTDIR/${i}.gtf -l ${i}_$LANE $MAPDIR/${i}.bam

done

###########################################################################################################################
# 2c. Stringtie assembly
###########################################################################################################################

## Estimate transcript abundances using published reference genome as merged reference 
# Also creates table counts for Ballgown (R)
for i in 1 .. 20;

do

stringtie -e -B -p 4 -G $REF_ANN -o $COUNTDIR/${i}.gtf $COUNTDIR/${i}.bam

done

###########################################################################################################################
# 2d. Stringtie to EdgeR
###########################################################################################################################

## Take the output from stringtie and convert it to gene and transcript count matrices for input for EdgeR
# Make sure you have a working version of python locally installed
module loady python

# Create a "merged.txt file" which is a text file that contains the location of your gtf files in the first column 
# and the name of the files pre .gtf in the second column
# Be sure the python script is in the same directory as the script you are running
prepDE.py -i $INDIR/merged.txt -g $COUNTDIR/gene_count_matrix.csv -t $OUTDIR/transcript_count_matrix.csv
