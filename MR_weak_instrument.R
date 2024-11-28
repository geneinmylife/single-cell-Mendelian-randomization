library(TwoSampleMR)#Wald ratio,Weighted median,Weighted mode
library(data.table)
library(vcfR)
library(stringr)
library(dplyr)
library(MRPRESSO)#mr_presso
library(mr.raps)#mr raps
library(MendelianRandomization)
library(readxl)
data_res_2 <- read_xlsx('cancer_1/data/check/Main_result.xlsx')
cancer_list <- unique(data_res_2$outcome)
data_res <- read_xlsx('cancer_1/data/check/NC_IV_final.xlsx')
fil_list <- c('ieu-a-1126',
              'ieu-b-85',
              'ieu-a-966',
              'ieu-a-1120',
              'thyroid_cancer',
              'cervix_cancer')




IV_data_list <- list.files('cancer_1/data/exposure_use_iv_exposure_associations/')

data_res$celltype <- str_split_fixed(data_res$exposure,'_',2)[,2]


weak_IV_list <- list()
for(i in 1:46){
  IV <- read.csv(paste0('cancer_1/data/exposure_use_iv_exposure_associations/',IV_data_list[i]))
  IV$celltype <- gsub('_500kb_combined_rsID_IV.csv','',IV_data_list[i])
  IV$exposure <- paste(IV$exposure,IV$celltype,sep = '_')
  weak_IV_list[[i]] <- IV
}
weak_IV_data <- do.call(rbind,weak_IV_list)


for(j in 6){
  load(paste0('cancer_1/data/summary/',fil_list[j],'.Rdata'))
  data_res_cancer <- data_res
  weak_IV_data_2 <- subset(weak_IV_data,exposure %in% data_res_cancer$exposure)
  weak_IV_data_2_formatted <- format_data(weak_IV_data_2,
                                          type = "exposure",
                                          phenotype_col = "exposure",
                                          snp_col = "SNP",
                                          beta_col = "beta.exposure",
                                          se_col = "se.exposure",
                                          eaf_col = "eaf.exposure",
                                          effect_allele_col = "effect_allele.exposure",
                                          other_allele_col = "other_allele.exposure",
                                          pval_col = "pval.exposure",
                                          samplesize_col = "samplesize.exposure",
                                          gene_col = "gene",
                                          id_col = "id.exposure",
                                          chr_col = "chr.exposure",
                                          pos_col = "pos.exposure")
  outcome_data2 <- subset(outcome_data1,snp %in% weak_IV_data_2_formatted$SNP)
  outcome_data2_formatted <- NULL
  if(j<6){
    outcome_data2_formatted <- format_data(outcome_data2,
                                           type = "outcome",
                                           phenotype_col = "cancer",
                                           snp_col = "snp",
                                           beta_col = "beta",
                                           se_col = "se",
                                           eaf_col = "EAF",
                                           effect_allele_col = "ALT",
                                           other_allele_col = "REF",
                                           pval_col = "P",
                                           chr_col = "CHROM",
                                           pos_col = "POS")
  } else{
    outcome_data2_formatted <- format_data(outcome_data2,
                                           type = "outcome",
                                           phenotype_col = "cancer",
                                           snp_col = "snp",
                                           beta_col = "beta",
                                           se_col = "se",
                                           eaf_col = "EAF",
                                           effect_allele_col = "ALT",
                                           other_allele_col = "REF",
                                           pval_col = "P")
  }
  
  
  dat <- harmonise_data(
    exposure_dat = weak_IV_data_2_formatted, 
    outcome_dat = outcome_data2_formatted
  )
  
  results_wr <- mr(dat,method_list = "mr_wald_ratio")
  # results_wr<-mr_singlesnp(dat)
  if(nrow(results_wr)>0){
    results_wr$outcome <- cancer_list[j]
    write.csv(results_wr,paste("LHY/cancer_1/data/check/wr/",cancer_list[j],'_',"wr.txt",sep=""),row.names=FALSE)
  }
  
  cd4_data_check<-dat[which(duplicated(dat$id.exposure)),]
  cd4_data_unique<-dat[-which(dat$id.exposure%in%cd4_data_check$id.exposure),]
  if(nrow(cd4_data_unique)>0){
    cd4_data_2plus<-dat[-which(dat$id.exposure%in%cd4_data_unique$id.exposure),]
  } else{
    cd4_data_2plus<-dat
  }
  
  
  results_wm<-mr(cd4_data_2plus,method_list = "mr_weighted_median")
  results_wm$outcome <- cancer_list[j]
  write.csv(results_wm,paste("cancer_1/data/check/wm/",cancer_list[j],'_',"wm.txt",sep=""),row.names=FALSE)
  
  
  results_wm2<-mr(cd4_data_2plus,method_list = "mr_weighted_mode")
  results_wm2$outcome <- cancer_list[j]
  write.csv(results_wm2,paste("cancer_1/data/check/wm2/",cancer_list[j],'_',"wm2.txt",sep=""),row.names=FALSE)
  
  exposure_names<-unique(cd4_data_2plus$exposure)
  
  
  for (k in 1:length(exposure_names)){
    cd4_data_2plus_use<-cd4_data_2plus[which(cd4_data_2plus$exposure==exposure_names[k]),]
    f <- try(results_raps<-mr.raps.all(b_exp = cd4_data_2plus_use$beta.exposure,
                                       b_out = cd4_data_2plus_use$beta.outcome,
                                       se_exp = cd4_data_2plus_use$se.exposure,
                                       se_out = cd4_data_2plus_use$se.outcome),silent = TRUE)
    if(class(f)=='try-error'){
      next
    }
    MR_RAPS_results<-cbind(exposure_names[k],nrow(cd4_data_2plus_use),
                           results_raps[1,1],results_raps[1,2],results_raps[1,3],results_raps[1,4])
    write.table(MR_RAPS_results,
                paste("cancer_1/data/check/RAPS/",cancer_list[j],'_',"RAPS.txt",sep=""),
                col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
  }
  for (l in 1:length(exposure_names)){
    data_for_MR<-as.data.frame(cd4_data_2plus[which(cd4_data_2plus$exposure==exposure_names[l]),])
    if(nrow(data_for_MR)>=4){
      # Run MR-PRESSO global method
      MR_PRESSO_results<-mr_presso(BetaOutcome = "beta.outcome", 
                                   BetaExposure = "beta.exposure", 
                                   SdOutcome = "se.outcome", 
                                   SdExposure = "se.exposure", 
                                   OUTLIERtest = TRUE, 
                                   DISTORTIONtest = TRUE, 
                                   data = data_for_MR, 
                                   NbDistribution = 1000,  
                                   SignifThreshold = 0.05)
      MR_PRESSO_results<-cbind(exposure_names[l],nrow(data_for_MR),
                               MR_PRESSO_results$`Main MR results`[1,2],MR_PRESSO_results$`Main MR results`[1,3],MR_PRESSO_results$`Main MR results`[1,4],MR_PRESSO_results$`Main MR results`[1,6],
                               MR_PRESSO_results$`Main MR results`[2,2],MR_PRESSO_results$`Main MR results`[2,3],MR_PRESSO_results$`Main MR results`[2,4],MR_PRESSO_results$`Main MR results`[2,6])
      write.table(MR_PRESSO_results,
                  paste("cancer_1/data/check/PRESSO/",cancer_list[j],'_',"PRESSO.txt",sep=""),
                  col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
    }
  }
  
  
  
  cd4_data_2plus$outcome<-cancer_list[j]
  cd4_data_2plus_mr<-dat_to_MRInput(cd4_data_2plus)
  for (k in 1:length(cd4_data_2plus_mr)){
    data_for_MR<-cd4_data_2plus_mr[[k]]
    f <- try(MR_CML<-mr_cML(data_for_MR,n=11900),silent = TRUE)
    if(class(f)=='try-error'){
      next
    }
    MR_CML_results<-cbind(names(cd4_data_2plus_mr[k]),MR_CML$SNPs,MR_CML$Estimate,MR_CML$StdError,MR_CML$CILower,MR_CML$CIUpper,MR_CML$Pvalue)
    write.table(MR_CML_results,
                paste("cancer_1/data/check/MR_CML/",cancer_list[j],'_',"MR_CML.txt",sep=""),
                col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
  }
  
  for (l in 1:length(cd4_data_2plus_mr)){
    data_for_MR<-cd4_data_2plus_mr[[l]]
    f <- try(MR_Robust<-mr_ivw(data_for_MR,robust = TRUE))
    if(class(f)=='try-error'){
      next
    }
    MR_Robust_results<-cbind(names(cd4_data_2plus_mr[l]),MR_Robust$SNPs,MR_Robust$Estimate,MR_Robust$StdError,MR_Robust$CILower,MR_Robust$CIUpper,MR_Robust$Pvalue)
    write.table(MR_Robust_results,
                paste("cancer_1/data/check/MRRobust/",cancer_list[j],'_',"MRRobust.txt",sep=""),
                col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
  }
  for (m in 1:length(cd4_data_2plus_mr)){
    data_for_MR<-cd4_data_2plus_mr[[m]]
    f <- try(MR_DIVW<-mr_divw(data_for_MR),silent = TRUE)
    if(class(f)=='try-error'){
      next
    }
    MR_DIVW_results<-cbind(names(cd4_data_2plus_mr[m]),MR_DIVW$SNPs,MR_DIVW$Estimate,MR_DIVW$StdError,MR_DIVW$CILower,MR_DIVW$CIUpper,MR_DIVW$Pvalue)
    write.table(MR_DIVW_results,
                paste("cancer_1/data/check/MRDIVW/",cancer_list[j],'_',"MRDIVW.txt",sep=""),
                col.names=FALSE,append=TRUE,row.names = FALSE,quote=FALSE)
  }
  
}
