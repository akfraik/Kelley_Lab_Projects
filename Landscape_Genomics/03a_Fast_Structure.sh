###########################################################################################################################
## 3a. Run Fast Structure
  # Ran two programs quantifying underying demographic structure 
###########################################################################################################################

## First step, is to take our converted BED files (pre and post-disease)
INDIR="/Structure_Analyses"
OUTDIR="/Structure_Analyses/Fast_Structure"
POPULATIONS="1 2 3 4 5 6 7"     # Replace numbers with your population names

### Pre-disease
## Fast Structure; run K 1:20 to see what the K value (putative population size) could be pre-disease
  for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
        do python structure.py -K $K --input "${INDIR}/bed_pre_dftd_fast_struc" --output "${OUTDIR}/pre_dftd_fast_struc_output" \
           --prior simple --full --seed 100
  done

## Choosing K - spits out the most likely K value based on the model complexity that maximizes marginal likelihood 
 # and the model components used to explain structure in data (pre-disease)
        python chooseK.py --input "${OUTDIR}/pre_dftd_fast_struc_output"

#Fast Structure Figures
  for K in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
       do python distruct.py -K $K --input "${OUTDIR}/pre_dftd_fast_struc_output" --output "${OUTDIR}/pre_dftd_fast_struc_distruct".$K 
  done
  
### Post-disease
## Fast Structure; run K 1:20 to see what the K value (putative population size) could be post-disease
  for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
        do python structure.py -K $K --input "${INDIR}/bed_post_dftd_admix" --output "${OUTDIR}/post_dftd_fast_struc_output" \
           --prior simple --full --seed 100
  done

## Choosing K - spits out the most likely K value based on the model complexity that maximizes marginal likelihood 
 # and the model components used to explain structure in data (post-disease)
        python chooseK.py --input "${OUTDIR}/post_dftd_fast_struc_output"

#Fast Structure Figures
  for K in 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20;
       do python distruct.py -K $K --input "${OUTDIR}/post_dftd_fast_struc_output" --output "${OUTDIR}/post_dftd_fast_struc_distruct".$K 
  done
