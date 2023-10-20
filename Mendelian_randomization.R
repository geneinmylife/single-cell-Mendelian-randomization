#library(TwoSampleMR)

###exposure input file
The exposure input file (single_cell_file) need to be a tab separated table with the following columns.
1. `CHR` -- chromosome identifier
2. `POS` -- base pair position
3. `RSID` -- name of the SNP
4. `REF_ALLELE` -- reference allele (usually the major allele)
5. `EFFECT_ALLELE` -- alternative allele (usually the minor allele)
6. `BETA` -- beta effect size of the marginal estimate
7. `SE` -- standard error of the effect size of the marginal estimate.
10. `MAF` -- allele frequency of the alternative allele (alternative allele)
11. `SAMPLE` -- number of observations used for the estimate.

#single_cell <- read_exposure_data("single_cell_file")

###Format your exposure data
#exposure_dat <- format_data(single_cell, type = "exposure")

###outcome input file (local files)
#outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = file,
    sep = ",",
    snp_col = "rsid",
    beta_col = "effect",
    se_col = "SE",
    effect_allele_col = "a1",
    other_allele_col = "a2",
    eaf_col = "a1_freq",
    pval_col = "p-value",
    units_col = "Units",
    gene_col = "Gene",
    samplesize_col = "n"
)

#harmdat <- harmonise_data(exposure_dat, outcome_dat, action=2)

###run MR analysis
#mr_results <- mr(harmdat,method_list=c("mr_wald_ratio"))
