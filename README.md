# Z-lab

#library(TwoSampleMR)
#library("readxl")
#exposure_dat <- read_exposure_data(
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
#mr_results <- mr(harmdat,method_list=c("mr_wald_ratio"))
