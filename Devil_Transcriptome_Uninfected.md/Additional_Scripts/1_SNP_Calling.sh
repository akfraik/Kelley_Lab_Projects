###########################################################################################################################
# 1. GATK SNP calling pipeline
# 1a. BAM clean up and renaming RG
###########################################################################################################################

### First change directories to the directory with your mapped files that have been converted to BAMs
cd $MAPDIR

## Start with output from the output from the script 2_Map_and_Trim.sh in the Devil_Trancriptome_Uninfected directory
module load java

for i in 1..20;

do

