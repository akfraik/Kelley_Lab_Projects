###########################################################################################################################
## 3b. Run Admixture
  # Ran two programs quantifying underying demographic structure 
  # SCRIPT NEEDS WORK ********
###########################################################################################################################
#!/usr/bin/env Rscript

## I would run this script from the inside of the directory your input files are in 
 # cd "/Structure_Analyses"

## Admixture can be run in R so we are going to use lapply
 # Usage: admixture <input file> <K>

module load admixture

admixture /project/devils_alexf/try3/data/Rapture/Structural_Analyses_Input/Admixture_Input/rap_bed_pre_admix_ChrNum.bed 2 --cv

# This needs to be fixed next time I run this - I think I need to rewrite this in lapply format...
  for (k in 1:20) {
    pdf("/Structure_Analyses/Admixture/bed_pre_dftd$K.pdf")
  tbl.k<-read.table("/project/devils_alexf/try2/results/Rapture/Admixture/rap_bed_pre_ChrNum.2.Q")
  *****<-barplot(t(as.matrix(tbl2)),col=rainbow(2),xlab="Individuals",ylab="Ancestry",border=NA)
  dev.off()
    print(k)
}

