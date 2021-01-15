###########################################################################################################################
## 6a. Outflank Input file generation
###########################################################################################################################

## First step, is to take our converted BED files (pre and post-disease)
INDIR="/Plink_Files"
OUTDIR="/Outflank"
POPULATIONS="1 2 3 4 5 6 7"     # Replace numbers with your population names
SCRIPTDIR="/Scripts"            # Specify the directory you have been storing all of your scripts (including this one) in

for i in pre post;

POPFILE="/Structure_Analyses/Devils_Locations_${i}.tsv"

do

## Make copies of all of the files you are manipulating so you don't tamper the original file
 # Put them in your output directory
  sort -k3 $POPFILE
   cp "${INDIR}/ped_${i}_disease.ped" "${OUTDIR}/ped_${i}_disease.ped"
   cp "${INDIR}/ped_${i}_disease.map" "${OUTDIR}/ped_${i}_disease.map"
   cp "${INDIR}/ped_${i}_disease.nosex" "${OUTDIR}/ped_${i}_disease.nosex"
  sort -k2 "${OUTDIR}/ped_${i}_disease.ped" > "${OUTDIR}/ped_${i}_disease.ped"

  mv "${OUTDIR}/ped_${i}_disease.ped" "${OUTDIR}/ped_${i}_disease.ped"

if [[ $PBS_O_WORKDIR ]]; then

## Load plink
    source /usr/share/modules/init/bash
    module load plink/1.90_extrachr  
  set -o pipefail

## Filter the plink file and change the format
    plink --file "${OUTDIR}/ped_${i}_disease.ped" --recode A tabx \
      --out "${OUTDIR}/outflank_${i}" --allow-extra-chr \
     || { echo "plink conversion to raw format failed"; exit 1; }

## Modify the plink file so that the first column contains population names
 # Start by making a file with population names and individual names
    rm -f "${OUTDIR}/population_table_${i}.tsv.tmp"
 
 # Replace 1:7 with population names   
 for pop in "1 2 3 4 5 6 7"; 
  do
     p=$(echo $pop)
    rep -P "\t${pop}\t" $POPFILE | awk '{print "'$p'\tsample_" $3}' >> "${OUTDIR}/population_table_${i}.tsv.tmp"
 done

 ## Sort the output and filter out the unique entries
    sort -k 2 "${OUTDIR}/population_table_${i}.tsv.tmp" | uniq > "${OUTDIR}/population_table_${i}.tsv" || \
  { echo "conversion to of population_table_${i} from tmp to tsv failed"; exit 1; }

## Remove the temporary file
    rm -f "${OUTDIR}/population_table_${i}.tsv.tmp"
    
## Then join this to the raw file; note that the result has the population IDs in the second column 
 # they really should be in the first column
    tail -n +2 "${OUTDIR}/outflank_${i}" | join -1 2 -2 2 "${OUTDIR}/population_table_${i}.tsv" - | \
      sed 's/ /\t/g' | cut -f 1,2,3- > "${OUTDIR}/outflank_${i}.populations.raw" || \
    { echo "joining population information failed"; exit 1; }
    
 ## Then chop up the file with various Unix tools by taking the 2nd column of the "raw" file (which contains the population IDs) 
  # and pipes that into a tail command, which removes the first row 
  # (-n +2 means take everything starting from the second row). 
  # The first row is just a header.
    cut -f 2 "${OUTDIR}/outflank_${i}.populations.raw" > "${OUTDIR}/outflank_${i}.ped.popnames" || \
    { echo "Making popnames for ${transect} failed"; exit 1; }

## This command takes the first row, and removes the first 6 columns 
 # (-f "7-" means 7th column to the end; -d " " means the input is space-delimited). The header row contains the locus names.
 # The sed commands replace spaces with linebreaks and remove some junk added to the ends of the IDs.
    head -n 1 "${outfile}.raw" | cut -f "7-" | sed 's/\t/\n/g' | sed 's/_.$//g' | sed 's/_random//g' | \
        sed 's/chr.*_GL//g' > "${OUTDIR}/outflank_${i}.locusnames" || \
    { echo "Making locusnames for ${transect} failed"; exit 1; }

    # This command first removes the header row, then removes the first
    # six columns, replaces spaces with tabs, and makes missing data "9".
    tail -n +2 "${outfile}.raw" | cut -f "7-" | sed 's/ /\t/g' | \
        sed 's/NA/9/g' > "${OUTDIR}/outflank_${i}.snpmat" || \
    { echo "Making snpmat for ${transect} failed"; exit 1; }
    
else

 ## Now, we want to run the R script (in the directory) in this loop and 
  mkdir Ouftlank_Output
  mkdir Ouftlank_Output/${i}
  sed -e 's/Variable/${i}/g' "$SCRIPTDIR/6b_Outflank.R" > "$SCRIPTDIR/6b_Outflank_${i}.R"

 ## now give permission to run this script within this for loop
  chmod u+x "$SCRIPTDIR/6b_Outflank_${i}.R"
  "./$SCRIPTDIR/6b_Outflank_${i}.R"

done

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
