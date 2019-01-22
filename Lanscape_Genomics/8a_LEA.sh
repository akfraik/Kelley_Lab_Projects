###########################################################################################################################
## 8a. Landscape Ecological Association (LEA)
###########################################################################################################################

# Navigate into the directory with your input files
cd Plink_Files
SCRIPTDIR=/Scripts

## Create a pre and post-DFTD loop to run the next script in
 # Like with Outflank, replace variable with pre or post
for i in pre post;
do

## Create directories for your output files
 # Now, we want to run the R script (in the directory) in this loop and 
mkdir LEA
mkdir LEA/${i}
sed -e 's/Variable/${i}/g' "$SCRIPTDIR/8b_LEA_${i}.R" > "$SCRIPTDIR/8b_LEA_${i}.R"

 ## now give permission to run this script within this for loop
chmod u+x "$SCRIPTDIR/8b_LEA_${i}.R"
"./$SCRIPTDIR/8b_LEA_${i}.R"

done

else

# SLURM equivalent to qsub to create copies of the scripts and create output and error files
# Makes a copy of the script every time it runs. The only requirement is that you don't try to
# run this script with source; relies on the fact that the name of this script will be contained in $0 variable, 
# so it can just copy itself to a new directory. The "basename" command gets rid of the directory name part of a file. 
# The nice thing about this is that it leaves an unambiguous record of what you did every time the script runs.
    sfile="${INDIR}/script_copies/$(date '+%Y-%m-%d')-$(basename $0)"
    cp $0 $sfile
    cd "${INDIR}/working"
    qsub $sfile
    echo "${INDIR}/working"
