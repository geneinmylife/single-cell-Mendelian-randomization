The single-cell-Mendelian-randomization contains TwoSampleMR package, a Mendelian randomization (MR) method that efficiently identifies causal effects between single cell gene expression and complex traits.

#Installation
To install the latest version of TwoSampleMR, perform as normal:

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

To update the package just run the remotes::install_github("MRCIEU/TwoSampleMR") command again.

#Run MR analysis    
Firstly, to run the single-cell MR analysis, you need to input a exposure file with correct format:
The exposure input file need to be a tab separated table with the following columns.
1. `CHR` -- chromosome identifier
2. `POS` -- base pair position
3. `RSID` -- name of the SNP
4. `REF_ALLELE` -- reference allele (usually the major allele)
5. `EFFECT_ALLELE` -- alternative allele (usually the minor allele)
6. `BETA` -- beta effect size of the marginal estimate
7. `SE` -- standard error of the effect size of the marginal estimate.
10. `MAF` -- allele frequency of the alternative allele (alternative allele)
11. `SAMPLE` -- number of observations used for the estimate.

Secondly, preparing a phenotype input file (outcome):
The phenotype input files need to be organized in the following format：
“rsid,effect,SE,a1,a2,a1_freq,p-value,Units,Gene,n”

To extract the exposure SNPs from this data, we could use the following command:

outcome_dat <- read_outcome_data(
    snps = bmi_exp_dat$SNP,
    filename = "gwas_summary.csv",
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

Thirdly, it is important to harmonise the effects. This means that the effect of a SNP on the exposure and the effect of that SNP on the outcome must each correspond to the same allele.

dat <- harmonise_data(exposure_dat, outcome_dat)

Fourthly, once the exposure and outcome data are harmonised, we have effects and standard errors for each instrument SNP available for the exposure and outcome traits. We can use this information to perform Mendelian randomisation. To do this, simply run:
MR_result <- mr(dat)

#MR methods
The list of available MR methods can be obtained:

mr_method_list()
