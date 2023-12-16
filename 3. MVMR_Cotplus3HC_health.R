# This script was created by Jasmine Khouja 03.10.22. 

# The script conducts univariable and multivariable MR exploring the effects of
# nicotine and non-nicotine constituents of tobacco smoke (measured by 3HC plus
# cotinine [Cplus3HC] and cigarettes per day [CPD]) on health outcomes 
# (chronic obstructive pulmonary disease [COPD], coronary heart disease [CHD], 
# forced expiratory volume [FEV], forced vital capacity [FVC], heart rate [HR],
# and body mass index [BMI])

# The SNPs used in this script were selected based on the finding from GSCAN (Liu et al., 2019) and the cotinine consortium.
# The script calls a proxy searching script - contact Jasmine for access.

################################################################################
##### Load packages #####
# Note: some packages already installed 
################################################################################

library(usethis)
library(TwoSampleMR)
library(ieugwasr)
library(googleAuthR)
library(tidyverse)
library(stringr)
library(plyr)
library(dplyr)
library(forestplot)
library(gtable)
library(reshape)
library(gplots)
require(ggplot2)
library(ggplot2)
library(gridExtra)
library(grid)
library(extrafont)
library(plotly)
library(data.table)
library(curl)
library(MVMR)
library(MendelianRandomization)
library(simex)
library(writexl)

################################################################################
##### Set directory and memory #####
################################################################################

wd<-"folder for working directory"
exp1_wd<-"folder which contains exposure 1 (Cplus3HC) data"
exp1_file<-"insert file name of exposure 1 data here"
exp2_wd<-"folder which contains exposure 2 (CPD) datas"
exp2_file<-"insert file name of exposure 2 data here"
out_wd<-"folder for outputs"

setwd(wd)
memory.limit(size = 80000)

################################################################################
##### Set SNP lists for exposure 1 (exp1), exposure 2 (exp2) and both #####
# if there is a header, header=TRUE, if not, header=FALSE 
# ensure that the SNP column is labelled SNP
################################################################################

exp1_SNPlist<- read.table("COTPLUS3HC_SNPlist.txt", header=FALSE)
exp1_SNPlist<-rename(exp1_SNPlist, c("SNP" = "V1"))
exp2_SNPlist<- read.table("CPD_SNPlist.txt", header=FALSE)
exp2_SNPlist<-rename(exp2_SNPlist, c("SNP" = "V1"))

################################################################################
#####Extract exposure data for MR of exp1 and health outcomes #####
################################################################################
setwd(exp1_wd)
exp1_dat <- read_exposure_data(exp1_file,
                              clump = FALSE,
                              sep = "\t",
                              snp_col = "RSID",
                              beta_col = "beta",
                              se_col = "se",
                              eaf_col = "eaf",
                              effect_allele_col = "EA",
                              other_allele_col = "NEA",
                              pval_col = "p",
                              samplesize_col = "N",
                              min_pval = 1e-200,
                              log_pval = FALSE
)


##### Creating a SNP list for exp1 whereby p<5e-6 because only 3 ##### 
##### conditionally independent SNPs ##### 
##### 16 SNPs remain after clumping ##### 

exp1_SNPlist_5e6 <- exp1_dat[exp1_dat$pval.exposure <= 0.000005, ]
exp1_SNPlist_5e6 <- clump_data(
  exp1_SNPlist_5e6,
  clump_kb = 500,
  clump_r2 = 0.1,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

exp1_dat_mr <- format_data(
  exp1_dat,
  type = "exposure",
  snps = exp1_SNPlist$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

exp1_dat_mr_5e6 <- format_data(
  exp1_dat,
  type = "exposure",
  snps = exp1_SNPlist_5e6$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd(wd)
write.csv(exp1_dat_mr,"exp1_dat_mr.csv", row.names = FALSE)
write.csv(exp1_dat_mr_5e6,"exp1_dat_mr_5e6.csv", row.names = FALSE)
write.csv(exp1_dat,"exp1_dat.csv", row.names = FALSE)

################################################################################
##### Extract exposure data for MR of exp2 and health outcomes #####
#Note - this data excludes 23andMe and UKB to avoid sample overlap.
################################################################################

setwd(exp2_wd)
memory.limit(size = 80000)

exp2_dat <- read_exposure_data(exp2_file,
                              clump = FALSE,
                              sep = "\t",
                              snp_col = "RSID",
                              beta_col = "BETA",
                              se_col = "SE",
                              eaf_col = "AF",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              pval_col = "PVALUE",
                              samplesize_col = "N",
                              min_pval = 1e-200,
                              log_pval = FALSE
)

exp2_dat_mr <- format_data(
  exp2_dat,
  type = "exposure",
  snps = exp2_SNPlist$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd(wd)
write.csv(exp2_dat_mr,"exp2_dat_mr.csv", row.names = FALSE)
write.csv(exp2_dat,"exp2_dat.csv", row.names = FALSE)

################################################################################
##### Extract outcome data for MR #####
# Set outcomes, samples and sample sizes list for loops
# Extract outcome data
# Organise outcome names
################################################################################

outcomes<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
samples<-c('current', 'ever', 'former', 'never')
samp_sizes <- c(49721, 213341, 163620, 258056) # current, ever, former, never

setwd(out_wd)

mr_dfs_exp1 <- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp1[[outcomes[i]]][[samples[j]]] <- read_outcome_data(
      paste0(outcomes[i], "_", samples[j], "_imputed.txt/",outcomes[i], "_", samples[j], "_imputed.txt"),
      snps = exp1_dat_mr$SNP,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      pval_col = "P_BOLT_LMM_INF",
      samplesize_col = samp_sizes[j],
      min_pval = 1e-200,
      log_pval = FALSE
    )
  }
}

mr_dfs_exp1_5e6 <- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]] <- read_outcome_data(
      paste0(outcomes[i], "_", samples[j], "_imputed.txt/",outcomes[i], "_", samples[j], "_imputed.txt"),
      snps = exp1_dat_mr_5e6$SNP,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      pval_col = "P_BOLT_LMM_INF",
      samplesize_col = samp_sizes[j],
      min_pval = 1e-200,
      log_pval = FALSE
    )
  }
}

mr_dfs_exp2 <- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp2[[outcomes[i]]][[samples[j]]] <- read_outcome_data(
      paste0(outcomes[i], "_", samples[j], "_imputed.txt/",outcomes[i], "_", samples[j], "_imputed.txt"),
      snps = exp2_dat_mr$SNP,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      pval_col = "P_BOLT_LMM_INF",
      samplesize_col = samp_sizes[j],
      min_pval = 1e-200,
      log_pval = FALSE
    )
  }
}

outcome_names<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp1[[outcomes[i]]][[samples[j]]]["Phenotype"]<-NA
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$Phenotype<-outcome_names[i]
  }
}

outcome_names<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]["Phenotype"]<-NA
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$Phenotype<-outcome_names[i]
  }
}

outcome_names<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp2[[outcomes[i]]][[samples[j]]]["Phenotype"]<-NA
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$Phenotype<-outcome_names[i]
  }
}


################################################################################
##### Harmonising #####
################################################################################

harmonise_mr_dfs_exp1<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]] <- harmonise_data(exp1_dat_mr, mr_dfs_exp1[[outcomes[i]]][[samples[j]]])
  }
}

harmonise_mr_dfs_exp1_5e6<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]] <- harmonise_data(exp1_dat_mr_5e6, mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]])
  }
}

harmonise_mr_dfs_exp2<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]] <- harmonise_data(exp2_dat_mr, mr_dfs_exp2[[outcomes[i]]][[samples[j]]])
  }
}

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed_exp1<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    proxy_needed_exp1[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(exp1_SNPlist$SNP, harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$SNP))
  }
}

proxy_needed_exp1_5e6<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    proxy_needed_exp1_5e6[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(exp1_SNPlist_5e6$SNP, harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$SNP))
  }
}

proxy_needed_exp2<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    proxy_needed_exp2[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(exp2_SNPlist$SNP, harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$SNP))
  }
}

################################################################################
##### No proxies needed #####
################################################################################

################################################################################
##### Harmonising negative effect alleles #####
#Note: Must ensure the beta is positive for the MR results to be correct
# This means you need to flip the beta and alleles so that the effect is positive
################################################################################

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for(k in 1:length(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]][["beta.exposure"]])){
      if(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure[k]<0){harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.outcome[k] <- -1*harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.outcome[k] }
      if(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure[k]<0){harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure[k] <- -1*harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure[k] }
    }
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for(k in 1:length(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]][["beta.exposure"]])){
      if(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure[k]<0){harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome[k] <- -1*harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome[k] }
      if(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure[k]<0){harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure[k] <- -1*harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure[k] }
    }
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for(k in 1:length(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]][["beta.exposure"]])){
      if(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure[k]<0){harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.outcome[k] <- -1*harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.outcome[k] }
      if(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure[k]<0){harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure[k] <- -1*harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure[k] }
    }
  }
}


################################################################################
##### Save dataframes #####
#Note: This is not necessary but can be useful for sharing the datasets 
# Also can avoid rerunning time-consuming SNP extraction steps above
################################################################################

setwd(wd)
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    write.csv(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]], paste0("MR_exp1_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    write.csv(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]], paste0("MR_exp1_5e6_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    write.csv(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]], paste0("MR_exp2_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
  }
}

outcomes<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
samples<-c('current', 'ever', 'former', 'never')
samp_sizes <- c(49721, 213341, 163620, 258056) # current, ever, former, never

################################################################################
################################## MR ##########################################
################################################################################

################################################################################
##### Generate results inc. F and Q stats for heterogeneity #####
################################################################################

result_exp1<-list()
ivw_exp1<-list()
ptr_exp1<-list()
egger_exp1<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure<- as.numeric(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure)
    result_exp1[[outcomes[i]]][[samples[j]]] <- mr(
      harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]],
      parameters = default_parameters(),
      method_list = subset(mr_method_list(), use_by_default)$obj
    )
    result_exp1[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(result_exp1[[outcomes[i]]][[samples[j]]])
    
    harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]<- subset(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]], mr_keep)
    ivw_exp1[[outcomes[i]]][[samples[j]]]<-TwoSampleMR::mr_ivw(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                                              harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                              harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                              harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.outcome,  
                                                              parameters = default_parameters())
    ptr_exp1[[outcomes[i]]][[samples[j]]]<-data.frame(ivw_exp1[[outcomes[i]]][[samples[j]]]["Q"])
    egger_exp1[[outcomes[i]]][[samples[j]]]<-mr_egger_regression(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure,
                                                                harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                                harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                                harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.outcome, 
                                                                parameters)
    ptr_exp1[[outcomes[i]]][[samples[j]]][2,1]<-egger_exp1[[outcomes[i]]][[samples[j]]]["Q"]
    F = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure^2/harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure^2
    mF = mean(F)
    ptr_exp1[[outcomes[i]]][[samples[j]]][3,1]<-mF
  }
}

result_exp1_5e6<-list()
ivw_exp1_5e6<-list()
ptr_exp1_5e6<-list()
egger_exp1_5e6<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure<- as.numeric(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure)
    result_exp1_5e6[[outcomes[i]]][[samples[j]]] <- mr(
      harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]],
      parameters = default_parameters(),
      method_list = subset(mr_method_list(), use_by_default)$obj
    )
    result_exp1_5e6[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(result_exp1_5e6[[outcomes[i]]][[samples[j]]])
    
    harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]<- subset(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]], mr_keep)
    ivw_exp1_5e6[[outcomes[i]]][[samples[j]]]<-TwoSampleMR::mr_ivw(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                                                  harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                                  harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                                  harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.outcome,  
                                                                  parameters = default_parameters())
    ptr_exp1_5e6[[outcomes[i]]][[samples[j]]]<-data.frame(ivw_exp1_5e6[[outcomes[i]]][[samples[j]]]["Q"])
    egger_exp1_5e6[[outcomes[i]]][[samples[j]]]<-mr_egger_regression(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure,
                                                                    harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                                    harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                                    harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.outcome, 
                                                                    parameters)
    ptr_exp1_5e6[[outcomes[i]]][[samples[j]]][2,1]<-egger_exp1_5e6[[outcomes[i]]][[samples[j]]]["Q"]
    F = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure^2/harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure^2
    mF = mean(F)
    ptr_exp1_5e6[[outcomes[i]]][[samples[j]]][3,1]<-mF
  }
}

result_exp2<-list()
ivw_exp2<-list()
ptr_exp2<-list()
egger_exp2<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure<- as.numeric(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure)
    result_exp2[[outcomes[i]]][[samples[j]]] <- mr(
      harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]],
      parameters = default_parameters(),
      method_list = subset(mr_method_list(), use_by_default)$obj
    )
    result_exp2[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(result_exp2[[outcomes[i]]][[samples[j]]])
    
    
    harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]<- subset(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]], mr_keep)
    ivw_exp2[[outcomes[i]]][[samples[j]]]<-TwoSampleMR::mr_ivw(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                                              harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                              harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                              harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.outcome,  
                                                              parameters = default_parameters())
    ptr_exp2[[outcomes[i]]][[samples[j]]]<-data.frame(ivw_exp2[[outcomes[i]]][[samples[j]]]["Q"])
    egger_exp2[[outcomes[i]]][[samples[j]]]<-mr_egger_regression(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure,
                                                                harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                                harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                                harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.outcome, 
                                                                parameters)
    ptr_exp2[[outcomes[i]]][[samples[j]]][2,1]<-egger_exp2[[outcomes[i]]][[samples[j]]]["Q"]
    F = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure^2/harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure^2
    mF = mean(F)
    ptr_exp2[[outcomes[i]]][[samples[j]]][3,1]<-mF
  }
}

################################################################################
##### Calculate I2GX #####
# Isq(y, s) where y = vector of effects and s = vector of standard errors
# Apply simulation extrapolation SIMEX corrections to MR-Egger analysis where 
# estimates < 0.9 
# This would indicate the effect estimate is biased by 10% due to measurement
# error
# (Bowden, Del Greco, et al., 2016)
# (Lederer & K?chenhoff, 2006)
# exp2 is okay 0.967 but exp1 is not - 0.807 and 0.863
################################################################################

ISQ_exp1<-data.frame(matrix(NA, nrow = length(outcomes), ncol = length(samples)))

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ISQ_exp1[i,j]<-Isq(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure, harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure)
  }
}

ISQ_exp1_5e6<-data.frame(matrix(NA, nrow = length(outcomes), ncol = length(samples)))

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ISQ_exp1_5e6[i,j]<-Isq(harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure, harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure)
  }
}
ISQ_exp2<-data.frame(matrix(NA, nrow = length(outcomes), ncol = length(samples)))

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ISQ_exp2[i,j]<-Isq(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure, harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure)
  }
}

source("SIMEX_Cotplus3HC.R")

################################################################################
##### Forest plot results #####
##### Note: doesn't include simex corrections
################################################################################

forest<-list()
forest_exp1<- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    
    forest_exp1[[outcomes[i]]][[samples[j]]] = mr_forest(mr_input(bx = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                                                 bxse = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                                 by = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                                 byse = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.outcome),
                                                        methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE)
    forest = mr_forest(mr_input(bx = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                bxse = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                by = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                byse = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.outcome)) 
    forest[[outcomes[i]]][[samples[j]]] = forest + coord_cartesian(xlim=c(-2,2))
    
    mr_funnel(mr_input(bx = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                       bxse = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure, 
                       by = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                       byse = harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.outcome))
    
  }
}

forest_exp1_5e6<- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    
    forest_exp1_5e6[[outcomes[i]]][[samples[j]]] = mr_forest(mr_input(bx = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                                                     bxse = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                                     by = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                                     byse = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.outcome),
                                                            methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE)
    forest = mr_forest(mr_input(bx = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                bxse = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                by = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                byse = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.outcome)) 
    forest[[outcomes[i]]][[samples[j]]] = forest + coord_cartesian(xlim=c(-2,2))
    
    mr_funnel(mr_input(bx = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                       bxse = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.exposure, 
                       by = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                       byse = harmonise_mr_dfs_exp1_5e6[[outcomes[i]]][[samples[j]]]$se.outcome))
    
  }
}

forest_exp2<- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    
    forest_exp2[[outcomes[i]]][[samples[j]]] = mr_forest(mr_input(bx = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                                                 bxse = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                                                 by = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                                                 byse = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.outcome),
                                                        methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE)
    forest = mr_forest(mr_input(bx = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                                bxse = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure, 
                                by = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                                byse = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.outcome)) 
    forest[[outcomes[i]]][[samples[j]]] = forest + coord_cartesian(xlim=c(-2,2))
    
    mr_funnel(mr_input(bx = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure, 
                       bxse = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure, 
                       by = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.outcome, 
                       byse = harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.outcome))
    
  }
}

################################################################################
################################   MVMR  #######################################
################################################################################

################################################################################
######Merge exposure SNP lists #####
# Load data if not already loaded
################################################################################

setwd(wd)
exp1_dat_mr<-read.csv("exp1_dat_mr.csv", header = TRUE)
exp1_dat_mr_5e6<-read.csv("exp1_dat_mr_5e6.csv", header = TRUE)
exp2_dat_mr<-read.csv("exp2_dat_mr.csv", header = TRUE)
exp1_dat<-read.csv("exp1_dat.csv", header = TRUE)
exp2_dat<-read.csv("exp2_dat_.csv", header = TRUE)

exp1_instr = exp1_dat_mr[,c("SNP","effect_allele.exposure","se.exposure","pval.exposure","beta.exposure","exposure","mr_keep.exposure","id.exposure","samplesize.exposure")]
exp1_5e6_instr = exp1_dat_mr_5e6[,c("SNP","effect_allele.exposure","se.exposure","pval.exposure","beta.exposure","exposure","mr_keep.exposure","id.exposure","samplesize.exposure")]
exp2_instr = exp2_dat_mr[,c("SNP","effect_allele.exposure","se.exposure","pval.exposure","beta.exposure","exposure","mr_keep.exposure","id.exposure","samplesize.exposure")]

INSTR<- do.call("rbind", list(exp1_instr, exp2_instr))
INSTR_5e6<- do.call("rbind", list(exp1_5e6_instr, exp2_instr))

################################################################################
##### Check for overlapping SNPs between the exposure instruments #####
################################################################################

n_occur<- data.frame(table(INSTR$SNP))
n_occur[n_occur$Freq >1, ]
INSTR[INSTR$SNP %in% n_occur$Var1[n_occur$Freq >1], ]

n_occur_5e6<- data.frame(table(INSTR_5e6$SNP))
n_occur_5e6[n_occur_5e6$Freq >1, ]
INSTR_5e6[INSTR_5e6$SNP %in% n_occur_5e6$Var1[n_occur_5e6$Freq >1], ]

################################################################################
##### Extract instruments from the exp1 data set #####
################################################################################

exp1_dat_mvmr <- format_data(
  exp1_dat,
  type = "exposure",
  snps = INSTR$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

##### Change name of GWAS and check n SNPs #####
exp1_dat_mvmr$id.exposure<- "1"
str(exp1_dat_mvmr)

################################################################################
##### Extract 5e6 instruments from the exp1 data set #####
################################################################################

exp1_dat_mvmr_5e6 <- format_data(
  exp1_dat,
  type = "exposure",
  snps = INSTR_5e6$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

##### Change name of GWAS and check n SNPs #####
exp1_dat_mvmr_5e6$id.exposure<- "1"
str(exp1_dat_mvmr_5e6)

################################################################################
##### Extract instruments from the exp2 data set #####
################################################################################

exp2_dat_mvmr <- format_data(
  exp2_dat,
  type = "exposure",
  snps = INSTR$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

##### Change name of GWAS and check n SNPs #####
exp2_dat_mvmr$id.exposure<- "2"
str(exp2_dat_mvmr)

################################################################################
##### Extract 5e6 instruments from the exp2 data set #####
################################################################################

exp2_dat_mvmr_5e6 <- format_data(
  exp2_dat,
  type = "exposure",
  snps = INSTR_5e6$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

##### Change name of GWAS and check n SNPs #####
exp2_dat_mvmr_5e6$id.exposure<- "2"
str(exp2_dat_mvmr_5e6)

################################################################################
##### Merge them to extract from outcome #####
################################################################################

exp1_dat_mvmr_1 = subset(exp1_dat_mvmr, select = c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "exposure", "mr_keep.exposure", "id.exposure", "eaf.exposure", "samplesize.exposure"))
exp2_dat_mvmr_1 = subset(exp2_dat_mvmr, select = c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "exposure", "mr_keep.exposure", "id.exposure", "eaf.exposure", "samplesize.exposure"))

exp1_dat_mvmr_5e6_1 = subset(exp1_dat_mvmr_5e6, select = c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "exposure", "mr_keep.exposure", "id.exposure", "eaf.exposure", "samplesize.exposure"))
exp2_dat_mvmr_5e6_1 = subset(exp2_dat_mvmr_5e6, select = c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "exposure", "mr_keep.exposure", "id.exposure", "eaf.exposure", "samplesize.exposure"))


##### Check structure is the same #####
str(exp1_dat_mvmr_1)
str(exp2_dat_mvmr_1)
str(exp1_dat_mvmr_5e6_1)
str(exp2_dat_mvmr_5e6_1)

################################################################################
##### Merge ####
################################################################################

exposures<- do.call("rbind", list(exp1_dat_mvmr_1, exp2_dat_mvmr_1))
exposures_5e6<- do.call("rbind", list(exp1_dat_mvmr_5e6_1, exp2_dat_mvmr_5e6_1))

##### Save dataframes #####
setwd(wd)
write.csv(exposures,"exp1_exp2.csv", row.names = FALSE)
write.csv(exposures_5e6,"exp1_exp2_5e6.csv", row.names = FALSE)

################################################################################
##### Find proxies missing from either the exp1 or exp2 dataset to add to ##### 
##### missing outcome SNPs #####
# Note: rs4886550 is not available in the ref panel 
# Proxy files are split into a and b where some SNPs are missing from the exp2 
# dataset and some are missing from the exp1 dataset
# 1a = exp2 snps missing from the exp1 data (no missing exp1 snps)
# 1b
# 2a = exp2 snps missing from the exp1 data (5e6 - essentially same as above)
# 2b = exp1 snps missing from the exp2 data
# No proxies needed
################################################################################

n_occur<- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq <2, ]
exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq <2], ]

n_occur_5e6<- data.frame(table(exposures_5e6$SNP))
n_occur_5e6[n_occur_5e6$Freq <2, ]
exposures_5e6[exposures_5e6$SNP %in% n_occur_5e6$Var1[n_occur_5e6$Freq <2], ]

proxy_needed1 <- data.frame(exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq <2], ])
proxy_needed1 <- data.frame(proxy_needed1[1])
proxy_needed2 <- data.frame(exposures_5e6[exposures_5e6$SNP %in% n_occur_5e6$Var1[n_occur_5e6$Freq <2], ])
proxy_needed2a<- subset(proxy_needed2, id.exposure==id.exposure[2])
proxy_needed2b<- subset(proxy_needed2, id.exposure==id.exposure[1])
proxy_needed2 <- data.frame(proxy_needed2[1])
proxy_needed2a <- data.frame(proxy_needed2a[1])
proxy_needed2b <- data.frame(proxy_needed2b[1])

################################################################################
##### Clumping ####
################################################################################

##### Change all p-values for exp1 to 1e-200 for clumping so that none are dropped #####
##### Save old p-values first #####
exposures$oldpvalues <-exposures$pval.exposure
exposures_5e6$oldpvalues <-exposures_5e6$pval.exposure

exposures<- exposures %>% 
  mutate(pval.exposure = if_else(exposures$SNP %in% exp1_instr$SNP, 1e-201, pval.exposure))

exposures_5e6<- exposures_5e6 %>% 
  mutate(pval.exposure = if_else(exposures_5e6$SNP %in% exp1_5e6_instr$SNP, 1e-201, pval.exposure))

##### Clump the data #####
exposures$id.exposure[exposures$id.exposure == "2"] <- "1"
exposures<- clump_data(exposures, clump_kb=500, clump_r2=0.1) 
str(exposures) 

exposures_5e6$id.exposure[exposures_5e6$id.exposure == "2"] <- "1"
exposures_5e6<- clump_data(exposures_5e6, clump_kb=500, clump_r2=0.1) 
str(exposures_5e6) 

##### Add ID's back #####
exposures$id.exposure[exposures$samplesize.exposure<6000] <- "1"
exposures$id.exposure[exposures$samplesize.exposure>6000] <- "2"

exposures_5e6$id.exposure[exposures_5e6$samplesize.exposure<6000] <- "1"
exposures_5e6$id.exposure[exposures_5e6$samplesize.exposure>6000] <- "2"

##### Revert all p-values for exp1 from 1e-200 #####
exposures$pval.exposure<-exposures$oldpvalues
exposures<-select(exposures,-c(oldpvalues))

exposures_5e6$pval.exposure<-exposures_5e6$oldpvalues
exposures_5e6<-select(exposures_5e6,-c(oldpvalues))

##### Split again to harmonise based on exposure id #####
exp1 = split(exposures, exposures$id.exposure)[['1']]
exp2 = split(exposures, exposures$id.exposure)[['2']]

exp1_5e6 = split(exposures_5e6, exposures_5e6$id.exposure)[['1']]
exp2_5e6 = split(exposures_5e6, exposures_5e6$id.exposure)[['2']]

################################################################################
##### Harmonise exp1 on exp2 #####
################################################################################

names(exp1) = gsub( "exposure", "outcome", names(exp1))
exp1_exp2 = harmonise_data(exp2, exp1)

names(exp1_5e6) = gsub( "exposure", "outcome", names(exp1_5e6))
exp1_exp2_5e6 = harmonise_data(exp2_5e6, exp1_5e6)

################################################################################
##### Keep only snps that are present across both exposures ####
#Note: they would have frequency 1 if only available in one dataset
################################################################################

n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq == 2,]
exposures<- exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq == 2],]
str(exposures)

n_occur_5e6 <- data.frame(table(exposures_5e6$SNP))
n_occur_5e6[n_occur_5e6$Freq == 2,]
exposures_5e6<- exposures_5e6[exposures_5e6$SNP %in% n_occur_5e6$Var1[n_occur_5e6$Freq == 2],]
str(exposures_5e6)

################################################################################
##### Format exposures #####
################################################################################

##### Keep only snps MrKeep= TRUE #####
exp1_exp2 = exp1_exp2[exp1_exp2$mr_keep== TRUE, ]
str(exp1_exp2) 

exp1_exp2_5e6 = exp1_exp2_5e6[exp1_exp2_5e6$mr_keep== TRUE, ]
str(exp1_exp2_5e6) 

##### Split the tables - exp2 ##### 
exp2_H<- subset(exp1_exp2, id.exposure=="2", select= c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure))
exp2_H_5e6<- subset(exp1_exp2_5e6, id.exposure=="2", select= c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure))

##### Split the tables - exp1 ##### 
exp1_H<- subset(exp1_exp2, id.outcome=="1", select= c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome))
exp1_H_5e6<- subset(exp1_exp2_5e6, id.outcome=="1", select= c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome))

##### Turn exp1 from outcome to exposure to merge the datasets ##### 
names(exp1_H) <- gsub("outcome", "exposure", names(exp1_H))
Exposures_H<- merge(exp1_H, exp2_H, all= TRUE)
Exposures_H["Phenotype"]<- NA
Exposures_H$Phenotype[Exposures_H$id.exposure == 1] <- "exp1"
Exposures_H$Phenotype[Exposures_H$id.exposure == 2] <- "exp2"
str(Exposures_H)

names(exp1_H_5e6) <- gsub("outcome", "exposure", names(exp1_H_5e6))
Exposures_H_5e6<- merge(exp1_H_5e6, exp2_H_5e6, all= TRUE)
Exposures_H_5e6["Phenotype"]<- NA
Exposures_H_5e6$Phenotype[Exposures_H_5e6$id.exposure == 1] <- "exp1"
Exposures_H_5e6$Phenotype[Exposures_H_5e6$id.exposure == 2] <- "exp2"
str(Exposures_H_5e6)

################################################################################
##### Extract outcome data for MVMR #####
################################################################################

outcomes<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
samples<-c('current', 'ever', 'former', 'never')
samp_sizes <- c(49721, 213341, 163620, 258056) # current, ever, former, never
setwd(out_wd)

mvmr_dfs <- list()
mvmr_dfs_5e6 <- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dfs[[outcomes[i]]][[samples[j]]] <- read_outcome_data(
      paste0(outcomes[i], "_", samples[j], "_imputed.txt/",outcomes[i], "_", samples[j], "_imputed.txt"),
      snps = Exposures_H$SNP,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      pval_col = "P_BOLT_LMM_INF",
      samplesize_col = samp_sizes[j],
      min_pval = 1e-200,
      log_pval = FALSE
    )
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dfs_5e6[[outcomes[i]]][[samples[j]]] <- read_outcome_data(
      paste0(outcomes[i], "_", samples[j], "_imputed.txt/",outcomes[i], "_", samples[j], "_imputed.txt"),
      snps = Exposures_H_5e6$SNP,
      sep = "\t",
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "A1FREQ",
      effect_allele_col = "ALLELE1",
      other_allele_col = "ALLELE0",
      pval_col = "P_BOLT_LMM_INF",
      samplesize_col = samp_sizes[j],
      min_pval = 1e-200,
      log_pval = FALSE
    )
  }
}

################################################################################
##### Organise outcome #####
################################################################################

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dfs[[outcomes[i]]][[samples[j]]]["Phenotype"]<- NA
    mvmr_dfs[[outcomes[i]]][[samples[j]]]["Phenotype"]<- outcomes[i]
    mvmr_dfs_5e6[[outcomes[i]]][[samples[j]]]["Phenotype"]<- NA
    mvmr_dfs_5e6[[outcomes[i]]][[samples[j]]]["Phenotype"]<- outcomes[i]
  }
}

################################################################################
##### Harmonise with outcome #####
################################################################################

mvmr_dat<-list()
mvmr_dat_5e6<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dat[[outcomes[i]]][[samples[j]]]<- harmonise_data(Exposures_H, mvmr_dfs[[outcomes[i]]][[samples[j]]])
    mvmr_dat[[outcomes[i]]][[samples[j]]]<- mvmr_dat[[outcomes[i]]][[samples[j]]][mvmr_dat[[outcomes[i]]][[samples[j]]]$mr_keep== TRUE, ]
    str(mvmr_dat[[outcomes[i]]][[samples[j]]])
    mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]<- harmonise_data(Exposures_H_5e6, mvmr_dfs_5e6[[outcomes[i]]][[samples[j]]])
    mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]<- mvmr_dat_5e6[[outcomes[i]]][[samples[j]]][mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$mr_keep== TRUE, ]
    str(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]])    
  }
}   

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed_mvmr<-list()
proxy_needed_mvmr_5e6<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    proxy_needed_mvmr[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(exp1_H$SNP, mvmr_dat[[outcomes[i]]][[samples[j]]]$SNP))
    proxy_needed_mvmr_5e6[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(exp1_H_5e6$SNP, mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$SNP))
    
  }
}

#("proxy_search_loop_CancerMVMR2_exp1_Rep.R")

################################################################################
##### Save dataframes #####
################################################################################
setwd(wd)

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    write.csv(mvmr_dat[[outcomes[i]]][[samples[j]]], paste0("MVMR_dat_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
    write.csv(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]], paste0("MVMR_dat_5e6_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
  }
}

##### Can run from here if data frames are unchanged.##### 
#####Can skip if data is already loaded. #####
setwd(wd)

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dat[[outcomes[i]]][[samples[j]]]<-  read.csv(paste0("MVMR_dat_", outcomes[i], "_", samples[j], ".csv"))
    mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]<-  read.csv(paste0("MVMR_dat_5e6_", outcomes[i], "_", samples[j], ".csv"))
    
  }
}

################################################################################
##### Run MVMR #####
################################################################################

bX1<-list()
bX2<-list()
bXse1<-list()
bXse2<-list()
bY<-list()
bYse<-list()
SNP<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    bX1[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat[[outcomes[i]]][[samples[j]]]$beta.exposure[mvmr_dat[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    bX2[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat[[outcomes[i]]][[samples[j]]]$beta.exposure[mvmr_dat[[outcomes[i]]][[samples[j]]]$id.exposure== 2]))
    bXse1[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat[[outcomes[i]]][[samples[j]]]$se.exposure[mvmr_dat[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    bXse2[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat[[outcomes[i]]][[samples[j]]]$se.exposure[mvmr_dat[[outcomes[i]]][[samples[j]]]$id.exposure== 2]))
    bY[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat[[outcomes[i]]][[samples[j]]]$beta.outcome[mvmr_dat[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    bYse[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat[[outcomes[i]]][[samples[j]]]$se.outcome[mvmr_dat[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    SNP[[outcomes[i]]][[samples[j]]]<- c(mvmr_dat[[outcomes[i]]][[samples[j]]]$SNP[mvmr_dat[[outcomes[i]]][[samples[j]]]$id.exposure== 1])
    
  }
}

dfs<-list()
df_mvmr<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    dfs[[outcomes[i]]][[samples[j]]]<- data.frame(bX1[[outcomes[i]]][[samples[j]]], bXse1[[outcomes[i]]][[samples[j]]], bX2[[outcomes[i]]][[samples[j]]], bXse2[[outcomes[i]]][[samples[j]]], bY[[outcomes[i]]][[samples[j]]], bYse[[outcomes[i]]][[samples[j]]], SNP[[outcomes[i]]][[samples[j]]])
    df_mvmr[[outcomes[i]]][[samples[j]]]<- format_mvmr(dfs[[outcomes[i]]][[samples[j]]][, c(1, 3)], dfs[[outcomes[i]]][[samples[j]]][, 5], dfs[[outcomes[i]]][[samples[j]]][, c(2, 4)], dfs[[outcomes[i]]][[samples[j]]][, 6], dfs[[outcomes[i]]][[samples[j]]][, 7])
  }
}

set.seed(1234)
mod.MVMR<-list()
se_theta1MI.random<-list()
mod<-list()
mod_or<-list()
mod_or_ivw<-list()
res<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod.MVMR[[outcomes[i]]][[samples[j]]]<-lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]]-1, weights=bYse[[outcomes[i]]][[samples[j]]]^-2)
    se_theta1MI.random[[outcomes[i]]][[samples[j]]] = summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]]-1, weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$coef[1,2]/
      min(summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]]-1, weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$sigma,1)
    
    mod[[outcomes[i]]][[samples[j]]]<- summary(mod.MVMR[[outcomes[i]]][[samples[j]]])
    
    mod_or[[outcomes[i]]][[samples[j]]] <- coef(summary(mod.MVMR[[outcomes[i]]][[samples[j]]]))
    colnames(mod_or[[outcomes[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_or_ivw[[outcomes[i]]][[samples[j]]]<-as.data.frame(mod_or[[outcomes[i]]][[samples[j]]])
    mod_or_ivw[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(mod_or_ivw[[outcomes[i]]][[samples[j]]])
    res[[outcomes[i]]][[samples[j]]]<- ivw_mvmr(df_mvmr[[outcomes[i]]][[samples[j]]])
  }
}

##### Orientation exp1 #####

clist<-c("bX2", "bY")

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for (var in clist) {
      eval(parse(text=paste0(var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "<-ifelse(bX1", "[['", outcomes[i], "']]", "[['", samples[j], "']]", ">0,",var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", ",", var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "*-1)")))
    }
    bX1[[outcomes[i]]][[samples[j]]]<-abs(bX1[[outcomes[i]]][[samples[j]]])
  }
}

##### MVMR Egger #####
mod.MVMRME_exp1<-list()
se_theta1ME.random<-list()
mod_ME_exp1<-list()
mod_ME_or_exp1<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod.MVMRME_exp1[[outcomes[i]]][[samples[j]]] <- summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))
    se_theta1ME.random[[outcomes[i]]][[samples[j]]] = summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$coef[2,2]/
      min(summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$sigma,1)
    mod_ME_exp1[[outcomes[i]]][[samples[j]]]<- summary(mod.MVMRME_exp1[[outcomes[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod_ME_or_exp1[[outcomes[i]]][[samples[j]]] <- data.frame(mod.MVMRME_exp1[[outcomes[i]]][[samples[j]]][["coefficients"]])
    colnames(mod_ME_or_exp1[[outcomes[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_ME_or_exp1[[outcomes[i]]][[samples[j]]]<-as.data.frame(mod_ME_or_exp1[[outcomes[i]]][[samples[j]]])
    mod_ME_or_exp1[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(mod_ME_or_exp1[[outcomes[i]]][[samples[j]]])
  }
}

##### Orientation exp2 #####

clist<-c("bX1","bY")
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for (var in clist) {
      eval(parse(text=paste0(var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "<-ifelse(bX2", "[['", outcomes[i], "']]", "[['", samples[j], "']]", ">0,",var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", ",", var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "*-1)")))
    }
    bX2[[outcomes[i]]][[samples[j]]]<-abs(bX2[[outcomes[i]]][[samples[j]]])
  }
}

##### MVMR Egger #####
mod.MVMRME_exp2<-list()
se_theta1ME.random<-list()
mod_ME_exp2<-list()
mod_ME_or_exp2<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod.MVMRME_exp2[[outcomes[i]]][[samples[j]]] <- summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))
    se_theta1ME.random[[outcomes[i]]][[samples[j]]] = summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$coef[2,2]/
      min(summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$sigma,1)
    mod_ME_exp2[[outcomes[i]]][[samples[j]]]<- summary(mod.MVMRME_exp2[[outcomes[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod_ME_or_exp2[[outcomes[i]]][[samples[j]]] <- data.frame(mod.MVMRME_exp2[[outcomes[i]]][[samples[j]]][["coefficients"]])
    colnames(mod_ME_or_exp2[[outcomes[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_ME_or_exp2[[outcomes[i]]][[samples[j]]]<-as.data.frame(mod_ME_or_exp2[[outcomes[i]]][[samples[j]]])
    mod_ME_or_exp2[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(mod_ME_or_exp2[[outcomes[i]]][[samples[j]]])
  }
}

##### IVW 5e6 #####

bX1<-list()
bX2<-list()
bXse1<-list()
bXse2<-list()
bY<-list()
bYse<-list()
SNP<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    bX1[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure[mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    bX2[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$beta.exposure[mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$id.exposure== 2]))
    bXse1[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$se.exposure[mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    bXse2[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$se.exposure[mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$id.exposure== 2]))
    bY[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$beta.outcome[mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    bYse[[outcomes[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$se.outcome[mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$id.exposure== 1]))
    SNP[[outcomes[i]]][[samples[j]]]<- c(mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$SNP[mvmr_dat_5e6[[outcomes[i]]][[samples[j]]]$id.exposure== 1])
    
  }
}

dfs_5e6<-list()
df_mvmr_5e6<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    dfs_5e6[[outcomes[i]]][[samples[j]]]<- data.frame(bX1[[outcomes[i]]][[samples[j]]], 
                                                      bXse1[[outcomes[i]]][[samples[j]]], 
                                                      bX2[[outcomes[i]]][[samples[j]]], 
                                                      bXse2[[outcomes[i]]][[samples[j]]], 
                                                      bY[[outcomes[i]]][[samples[j]]], 
                                                      bYse[[outcomes[i]]][[samples[j]]], 
                                                      SNP[[outcomes[i]]][[samples[j]]])
    
    df_mvmr_5e6[[outcomes[i]]][[samples[j]]]<- format_mvmr(dfs_5e6[[outcomes[i]]][[samples[j]]][, c(1, 3)], 
                                                           dfs_5e6[[outcomes[i]]][[samples[j]]][, 5], 
                                                           dfs_5e6[[outcomes[i]]][[samples[j]]][, c(2, 4)], 
                                                           dfs_5e6[[outcomes[i]]][[samples[j]]][, 6],
                                                           dfs_5e6[[outcomes[i]]][[samples[j]]][, 7])
  }
}

mod.MVMR.5e6<-list()
mod_5e6<-list()
mod_or_5e6<-list()
mod_or_ivw_5e6<-list()
res_5e6<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod.MVMR.5e6[[outcomes[i]]][[samples[j]]]<-lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]]-1, weights=bYse[[outcomes[i]]][[samples[j]]]^-2)
    se_theta1MI.random[[outcomes[i]]][[samples[j]]] = summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]]-1, weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$coef[1,2]/
      min(summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]]-1, weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$sigma,1)
    
    mod_5e6[[outcomes[i]]][[samples[j]]]<- summary(mod.MVMR.5e6[[outcomes[i]]][[samples[j]]])
    
    mod_or_5e6[[outcomes[i]]][[samples[j]]] <- coef(summary(mod.MVMR.5e6[[outcomes[i]]][[samples[j]]]))
    colnames(mod_or_5e6[[outcomes[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]]<-as.data.frame(mod_or_5e6[[outcomes[i]]][[samples[j]]])
    mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]])
    res_5e6[[outcomes[i]]][[samples[j]]]<- ivw_mvmr(df_mvmr_5e6[[outcomes[i]]][[samples[j]]])
    
  }
}

##### Orientation exp1 #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for Cplus3HC, and then     #####
##### orient the betas to be positive for CPD. In the paper, we           #####
##### report the result for each exposure only with the right orientation #####

clist<-c("bX2", "bY")

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for (var in clist) {
      eval(parse(text=paste0(var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "<-ifelse(bX1", "[['", outcomes[i], "']]", "[['", samples[j], "']]", ">0,",var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", ",", var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "*-1)")))
    }
    bX1[[outcomes[i]]][[samples[j]]]<-abs(bX1[[outcomes[i]]][[samples[j]]])
  }
}

##### MVMR Egger #####
mod.MVMRME_exp1_5e6<-list()
se_theta1ME.random<-list()
mod_ME_exp1_5e6<-list()
mod_ME_or_exp1_5e6<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod.MVMRME_exp1_5e6[[outcomes[i]]][[samples[j]]] <- summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))
    se_theta1ME.random[[outcomes[i]]][[samples[j]]] = summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$coef[2,2]/
      min(summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$sigma,1)
    mod_ME_exp1_5e6[[outcomes[i]]][[samples[j]]]<- summary(mod.MVMRME_exp1_5e6[[outcomes[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]] <- data.frame(mod.MVMRME_exp1_5e6[[outcomes[i]]][[samples[j]]][["coefficients"]])
    colnames(mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]]<-as.data.frame(mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]])
    mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]])
  }
}

##### Orientation exp2 #####

clist<-c("bX1","bY")
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for (var in clist) {
      eval(parse(text=paste0(var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "<-ifelse(bX2", "[['", outcomes[i], "']]", "[['", samples[j], "']]", ">0,",var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", ",", var, "[['", outcomes[i], "']]", "[['", samples[j], "']]", "*-1)")))
    }
    bX2[[outcomes[i]]][[samples[j]]]<-abs(bX2[[outcomes[i]]][[samples[j]]])
  }
}

##### MVMR Egger #####
mod.MVMRME_exp2_5e6<-list()
se_theta1ME.random<-list()
mod_ME_exp2_5e6<-list()
mod_ME_or_exp2_5e6<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod.MVMRME_exp2_5e6[[outcomes[i]]][[samples[j]]] <- summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))
    se_theta1ME.random[[outcomes[i]]][[samples[j]]] = summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$coef[2,2]/
      min(summary(lm(bY[[outcomes[i]]][[samples[j]]]~bX1[[outcomes[i]]][[samples[j]]]+bX2[[outcomes[i]]][[samples[j]]], weights=bYse[[outcomes[i]]][[samples[j]]]^-2))$sigma,1)
    mod_ME_exp2_5e6[[outcomes[i]]][[samples[j]]]<- summary(mod.MVMRME_exp2_5e6[[outcomes[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]] <- data.frame(mod.MVMRME_exp2_5e6[[outcomes[i]]][[samples[j]]][["coefficients"]])
    colnames(mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]]<-as.data.frame(mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]])
    mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]])
  }
}


################################################################################
##### Calculate F-statistic and covariance #####
# Note: >10 is strong
# correlation between exp1 and exp2 in Buchwald  = 0.468
################################################################################

covmat<-list()
Fstat<-list()
cov <- matrix(c(1,0.468,0.468,1), nrow=2, ncol=2)
covmat_5e6<-list()
Fstat_5e6<-list()



for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    covmat[[outcomes[i]]][[samples[j]]] <- phenocov_mvmr(cov, df_mvmr[[outcomes[j]]][[samples[j]]][,c(6,7)] ) 
    Fstat[[outcomes[i]]][[samples[j]]] <- strength_mvmr(df_mvmr[[outcomes[j]]][[samples[j]]], gencov = covmat[[outcomes[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    covmat_5e6[[outcomes[i]]][[samples[j]]] <- phenocov_mvmr(cov, df_mvmr_5e6[[outcomes[j]]][[samples[j]]][,c(6,7)] ) 
    Fstat_5e6[[outcomes[i]]][[samples[j]]] <- strength_mvmr(df_mvmr_5e6[[outcomes[j]]][[samples[j]]], gencov = covmat_5e6[[outcomes[i]]][[samples[j]]])
  }
}

################################################################################
##### Test for horizontal pleiotropy #####
# Note: Q should be greater than the number of SNPs included
################################################################################

ptr<-list()
ptr_5e6<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ptr[[outcomes[i]]][[samples[j]]] <- pleiotropy_mvmr(df_mvmr[[outcomes[i]]][[samples[j]]], gencov = covmat[[outcomes[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ptr_5e6[[outcomes[i]]][[samples[j]]] <- pleiotropy_mvmr(df_mvmr_5e6[[outcomes[i]]][[samples[j]]], gencov = covmat_5e6[[outcomes[i]]][[samples[j]]])
  }
}

################################################################################
##### Q-statistic minimisation estimate #####
# qhet_mvmr is used to adjust for covariance but at present, CIs take substantial time to calculate and crash R
# Compare effects with and without adjustment
################################################################################

res_qhet<-list()
res_qhet_5e6<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    res_qhet[[outcomes[i]]][[samples[j]]] <- qhet_mvmr(df_mvmr[[outcomes[i]]][[samples[j]]], cov, CI = F, iterations = 100)
  }
}
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    res_qhet_5e6[[outcomes[i]]][[samples[j]]] <- qhet_mvmr(df_mvmr_5e6[[outcomes[i]]][[samples[j]]], cov, CI = F, iterations = 100)
  }
}
################################################################################
##### Forest Plots #####
# Only plotting 5e6 because egger for 5e8 because the I2GX was below 0.6
# The results could be severely impacted by NOME violation and adjusting may 
# be worse
################################################################################

Exposure<-list()
Method<-list()
Beta<-list()
SE<-list()
LCI<-list()
UCI<-list()
OR<-list()
ORLCI<-list()
ORUCI<-list()
p<-list()
Q<-list()
I2GX<-list()
EggerI<-list()
EggerIp<-list()
F_stat<-list()
AllRes<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    Exposure[[outcomes[i]]][[samples[j]]]<- c("COT plus 3HC","COT plus 3HC","COT plus 3HC","COT plus 3HC","CPD","CPD","CPD","CPD")
    
    Method[[outcomes[i]]][[samples[j]]]<- c("MR-IVW","MVMR-IVW","MR-Egger","MVMR-Egger","MR-IVW","MVMR-IVW","MR-Egger","MVMR-Egger")
    
    Beta[[outcomes[i]]][[samples[j]]] <- as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"b"],
                                                      mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]] [1,"b"],
                                                      simexegger_exp1_5e6[[outcomes[i]]][[samples[j]]] [1,"beta_weighted"],
                                                      mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]] [2,"b"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]] [3,"b"],
                                                      mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]] [2,"b"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]] [1,"b"],
                                                      mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"b"]))
    
    LCI[[outcomes[i]]][[samples[j]]] <- as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"lo_ci"],
                                                     mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]] [1,"lo_ci"], 
                                                     simexegger_exp1_5e6[[outcomes[i]]][[samples[j]]] [1,"lowerCI_weighted"],
                                                     mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]] [2,"lo_ci"],
                                                     result_exp2[[outcomes[i]]][[samples[j]]] [3,"lo_ci"],
                                                     mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]] [2,"lo_ci"],
                                                     result_exp2[[outcomes[i]]][[samples[j]]] [1,"lo_ci"],
                                                     mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"lo_ci"]))
    
    UCI[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"up_ci"],
                                                    mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]] [1,"up_ci"], 
                                                    simexegger_exp1_5e6[[outcomes[i]]][[samples[j]]] [1,"upperCI_weighted"],
                                                    mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]] [2,"up_ci"],
                                                    result_exp2[[outcomes[i]]][[samples[j]]] [3,"up_ci"],
                                                    mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]] [2,"up_ci"],
                                                    result_exp2[[outcomes[i]]][[samples[j]]] [1,"up_ci"],
                                                    mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"up_ci"]))  
    
    
    OR[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"or"],
                                                   mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][1,"or"],
                                                   simexegger_exp1_5e6[[outcomes[i]]][[samples[j]]] [1,"OR_weighted"],
                                                   mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]] [2,"or"],
                                                   result_exp2[[outcomes[i]]][[samples[j]]][3,"or"],
                                                   mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][2,"or"],
                                                   result_exp2[[outcomes[i]]][[samples[j]]][1,"or"],
                                                   mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"or"]))
    
    ORLCI[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"or_lci95"],
                                                      mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][1,"or_lci95"],
                                                      simexegger_exp1_5e6[[outcomes[i]]][[samples[j]]][1,"lCIOR_weighted"],
                                                      mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]][2,"or_lci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][3,"or_lci95"],
                                                      mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][2,"or_lci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][1,"or_lci95"],
                                                      mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"or_lci95"])) 
    
    ORUCI[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"or_uci95"],
                                                      mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][1,"or_uci95"],
                                                      simexegger_exp1_5e6[[outcomes[i]]][[samples[j]]][1,"uCIOR_weighted"],
                                                      mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]][2,"or_uci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][3,"or_uci95"],
                                                      mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][2,"or_uci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][1,"or_uci95"],
                                                      mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"or_uci95"]))
    
    SE[[outcomes[i]]][[samples[j]]]<-c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"se"],
                                       mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][1,"se"], 
                                       ".",
                                       mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]][2,"se"],
                                       result_exp2[[outcomes[i]]][[samples[j]]][3,"se"],
                                       mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][2,"se"],
                                       result_exp2[[outcomes[i]]][[samples[j]]][1,"se"],
                                       mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"se"])
    
    p[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"pval"],
                                                  mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][1,"p"],
                                                  simexegger_exp1_5e6[[outcomes[i]]][[samples[j]]][1,"p_weighted"],
                                                  mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]][2,"p"],
                                                  result_exp2[[outcomes[i]]][[samples[j]]][3,"pval"],
                                                  mod_or_ivw_5e6[[outcomes[i]]][[samples[j]]][2,"p"],
                                                  result_exp2[[outcomes[i]]][[samples[j]]][1,"pval"],
                                                  mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][3,"p"]))
    
    I2GX[[outcomes[i]]][[samples[j]]]<-c(".",".",ISQ_exp1_5e6[1,1],".",".",".",ISQ_exp2[1,1],".")
    
    Q[[outcomes[i]]][[samples[j]]]<-c(ptr_exp1_5e6[[outcomes[i]]][[samples[j]]][1,"Q"],
                                      ptr_5e6[[outcomes[i]]][[samples[j]]][["Qstat"]],
                                      ptr_exp1_5e6[[outcomes[i]]][[samples[j]]][2,"Q"],
                                      ptr_5e6[[outcomes[i]]][[samples[j]]][["Qstat"]],
                                      ptr_exp2[[outcomes[i]]][[samples[j]]][1,"Q"],
                                      ptr_5e6[[outcomes[i]]][[samples[j]]][["Qstat"]],
                                      ptr_exp2[[outcomes[i]]][[samples[j]]][2,"Q"],
                                      ptr_5e6[[outcomes[i]]][[samples[j]]][["Qstat"]])
    
    EggerI[[outcomes[i]]][[samples[j]]]<-c(".",".",egger_exp1_5e6[[outcomes[i]]][[samples[j]]][["b_i"]],
                                           mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]][1,"b"],
                                           ".",".",egger_exp2[[outcomes[i]]][[samples[j]]][["b_i"]],
                                           mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][1,"b"])
    
    EggerIp[[outcomes[i]]][[samples[j]]]<-c(".",".",egger_exp1_5e6[[outcomes[i]]][[samples[j]]][["pval_i"]],
                                            mod_ME_or_exp1_5e6[[outcomes[i]]][[samples[j]]][1,"p"],
                                            ".",".",egger_exp2[[outcomes[i]]][[samples[j]]][["pval_i"]],
                                            mod_ME_or_exp2_5e6[[outcomes[i]]][[samples[j]]][1,"p"])
    
    F_stat[[outcomes[i]]][[samples[j]]]<-c(ptr_exp1_5e6[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat_5e6[[outcomes[i]]][[samples[j]]][1,1],
                                           ptr_exp1_5e6[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat_5e6[[outcomes[i]]][[samples[j]]][1,1],
                                           ptr_exp2[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat_5e6[[outcomes[i]]][[samples[j]]][1,2],
                                           ptr_exp2[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat_5e6[[outcomes[i]]][[samples[j]]][1,2])
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    AllRes[[outcomes[i]]][[samples[j]]]<-data.frame(
      Exposure = Exposure[[outcomes[i]]][[samples[j]]],
      Method = Method[[outcomes[i]]][[samples[j]]],
      Beta = Beta[[outcomes[i]]][[samples[j]]],
      SE = SE[[outcomes[i]]][[samples[j]]],
      LCI = LCI[[outcomes[i]]][[samples[j]]],
      UCI = UCI[[outcomes[i]]][[samples[j]]],
      OR = OR[[outcomes[i]]][[samples[j]]],
      ORLCI = ORLCI[[outcomes[i]]][[samples[j]]],
      ORUCI = ORUCI[[outcomes[i]]][[samples[j]]],
      p = p[[outcomes[i]]][[samples[j]]],
      I2GX = I2GX[[outcomes[i]]][[samples[j]]],
      Q = Q[[outcomes[i]]][[samples[j]]],
      EggerI = EggerI[[outcomes[i]]][[samples[j]]],
      EggerIp = EggerIp[[outcomes[i]]][[samples[j]]],
      F_stat = F_stat[[outcomes[i]]][[samples[j]]])
  }
}

###### Save results #####
setwd(wd)
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    write.csv(AllRes[[outcomes[i]]][[samples[j]]], paste0("Results_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
  }
}

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    AllRes[[outcomes[i]]][[samples[j]]] <- read.csv(paste0("Results_", outcomes[i], "_", samples[j], ".csv"), header = TRUE)
  }
}

nSNPs<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    nSNPs[[outcomes[i]]][[samples[j]]]<-as.numeric(c(result_exp1_5e6[[outcomes[i]]][[samples[j]]][3,"nsnp"],
                                                     nrow(mvmr_dfs_5e6[[outcomes[i]]][[samples[j]]]),
                                                     result_exp1_5e6[[outcomes[i]]][[samples[j]]][1,"nsnp"],
                                                     length(mod.MVMRME_exp1_5e6[[outcomes[i]]][[samples[j]]][["residuals"]]),
                                                     result_exp2[[outcomes[i]]][[samples[j]]][3,"nsnp"],
                                                     nrow(mvmr_dfs_5e6[[outcomes[i]]][[samples[j]]]),
                                                     result_exp2[[outcomes[i]]][[samples[j]]][1,"nsnp"],
                                                     length(mod.MVMRME_exp1_5e6[[outcomes[i]]][[samples[j]]][["residuals"]])))
    
  }
}

##### Betas #####
# Create dataframe including the effect estimates, lower and upper confidence intervals

##### For each outcome and sample, produce a  forest plot split for OR and beta #####
setwd(wd)

outcomes_OR<- c("COPD", "CHD")
outcomes_b<- c("BMI", "FEV", "FVC", "HR")

cochrane_from_rmeta_OR<-list()
OR_LCI_UCI<-list()
tabletext_OR<-list()

for (i in 1:length(outcomes_OR)) {
  for (j in 1:length(samples)) {
    cochrane_from_rmeta_OR[[outcomes_OR[i]]][[samples[j]]]<- 
      structure(list(
        mean  = c(NA,AllRes[[outcomes_OR[i]]][[samples[j]]]$OR), 
        lower = c(NA,AllRes[[outcomes_OR[i]]][[samples[j]]]$ORLCI),
        upper = c(NA,AllRes[[outcomes_OR[i]]][[samples[j]]]$ORUCI)), 
        .Names = c("mean", "lower", "upper"), 
        row.names = c(NA, -9), 
        class = "data.frame")
    
  }
}

cochrane_from_rmeta_b<-list()
b_LCI_UCI<-list()
tabletext_b<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    cochrane_from_rmeta_b[[outcomes_b[i]]][[samples[j]]]<- 
      structure(list(
        mean  = c(NA,AllRes[[outcomes_b[i]]][[samples[j]]]$Beta), 
        lower = c(NA,AllRes[[outcomes_b[i]]][[samples[j]]]$LCI),
        upper = c(NA,AllRes[[outcomes_b[i]]][[samples[j]]]$UCI)), 
        .Names = c("mean", "lower", "upper"), 
        row.names = c(NA, -9), 
        class = "data.frame")
  }
}

for (i in 1:length(outcomes_OR)) {
  for (j in 1:length(samples)) {
    OR_LCI_UCI[[outcomes_OR[i]]][[samples[j]]]<-str_c(round(AllRes[[outcomes_OR[i]]][[samples[j]]]$OR, digits=2),
                                                      ", (",
                                                      round(AllRes[[outcomes_OR[i]]][[samples[j]]]$ORLCI, digits=2),
                                                      ", ",
                                                      round(AllRes[[outcomes_OR[i]]][[samples[j]]]$ORUCI, digits = 2),")")
    
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    b_LCI_UCI[[outcomes_OR[i]]][[samples[j]]]<-str_c(round(AllRes[[outcomes_b[i]]][[samples[j]]]$Beta, digits=2),
                                                     ", (",
                                                     round(AllRes[[outcomes_b[i]]][[samples[j]]]$LCI, digits=2),
                                                     ", ",
                                                     round(AllRes[[outcomes_b[i]]][[samples[j]]]$UCI, digits = 2),")")
    
  }
}

for (i in 1:length(outcomes_OR)) {
  for (j in 1:length(samples)) {
    tabletext_OR[[outcomes_OR[i]]][[samples[j]]]<-cbind(
      c("Exposure", "3HC plus Cotinine", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
      c("N SNPs", nSNPs[[outcomes_OR[i]]][[samples[j]]]),
      c("Method", "MR IVW", "MVMR IVW", "MR-Egger (SIMEX)", "MVMR-Egger", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"),
      c("OR (95% CI)", OR_LCI_UCI[[outcomes_OR[i]]][[samples[j]]]),
      c("P value", round(AllRes[[outcomes_OR[i]]][[samples[j]]]$p, digits=3)))
    tabletext_OR[[outcomes_OR[i]]][[samples[j]]][,5][tabletext_OR[[outcomes_OR[i]]][[samples[j]]][,5]==0]<- "<0.001"
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    tabletext_b[[outcomes_b[i]]][[samples[j]]]<-cbind(
      c("Exposure", "3HC plus Cotinine", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
      c("N SNPs", nSNPs[[outcomes_b[i]]][[samples[j]]]),
      c("Method", "MR IVW", "MVMR IVW", "MR-Egger (SIMEX)", "MVMR-Egger", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"),
      c("beta (95% CI)", b_LCI_UCI[[outcomes_OR[i]]][[samples[j]]]),
      c("P value", round(AllRes[[outcomes_b[i]]][[samples[j]]]$p, digits=3)))
    tabletext_b[[outcomes_b[i]]][[samples[j]]][,5][tabletext_b[[outcomes_b[i]]][[samples[j]]][,5]==0]<- "<0.001"
  }
}

xtick_l<-function(lower_limit 
){
  if (lower_limit < 0.1) {
    xtick_lower<-0
  } else if (lower_limit < 0.2){
    xtick_lower<-0.1
  } else if (lower_limit < 0.3){
    xtick_lower<-0.2
  } else if (lower_limit < 0.4){
    xtick_lower<-0.3
  } else if (lower_limit < 0.5){
    xtick_lower<-0.4
  } else if (lower_limit < 0.6){
    xtick_lower<-0.5
  } else if (lower_limit < 0.7){
    xtick_lower<-0.6
  } else if (lower_limit < 0.8){
    xtick_lower<-0.7
  } else if (lower_limit < 0.9){
    xtick_lower<-0.8
  } else if (lower_limit < 1){
    xtick_lower<-0.9
  } else {
    xtick_lower<-floor(lower_limit)
  }
}

xtick_u<-function(upper_limit
) {
  if (upper_limit>2){
    xtick_upper <- ceiling(upper_limit)
  } else if (upper_limit < 0.1) {
    xtick_upper<-0.1
  } else if (upper_limit < 0.2) {
    xtick_upper<-0.2
  } else if (upper_limit < 0.3){
    xtick_upper<-0.3
  } else if (upper_limit < 0.4){
    xtick_upper<-0.4
  } else if (upper_limit < 0.5){
    xtick_upper<-0.5
  } else if (upper_limit < 0.6){
    xtick_upper<-0.6
  } else if (upper_limit < 0.7){
    xtick_upper<-0.7
  } else if (upper_limit < 0.8){
    xtick_upper<-0.8
  } else if (upper_limit < 0.9){
    xtick_upper<-0.9
  } else if (upper_limit < 1){
    xtick_upper<-1
  } else if (upper_limit <1.2){
    xtick_upper<-1.2
  } else if (upper_limit <1.5){
    xtick_upper<-1.5
  } else if (upper_limit <=2){
    xtick_upper<-2
  }
}

##### Create graphs and write to file #####
for (i in 1:length(outcomes_OR)) {
  for (j in 1:length(samples)) {
    tabletext<-matrix(tabletext_OR[[outcomes_OR[i]]][[samples[j]]], nrow=9, ncol=5)
    cochrane_from_rmeta<-data.frame(cochrane_from_rmeta_OR[[outcomes_OR[i]]][[samples[j]]])
    xticks<-c(xtick_l(min(cochrane_from_rmeta$lower, na.rm = TRUE)), 1, xtick_u(max(cochrane_from_rmeta$upper, na.rm = TRUE)))
    pdf.options(reset = TRUE, onefile = FALSE)
    pdfname<-paste0("MVMR_",outcomes_OR[i],"_",samples[j],".pdf")
    pdf(file=pdfname, width=15,height=15)
    plot(forestplot(tabletext,cochrane_from_rmeta,
                    new_page = TRUE,
                    is.summary=FALSE,
                    lineheight = unit(1, "cm"),
                    graphwidth=unit(10, "cm"),
                    boxsize=0.15,
                    clip=c(0.001,Inf), 
                    xticks=xticks,
                    xlog= TRUE,
                    zero=1,
                    ci.vertices=TRUE,
                    colgap=unit(4, "mm"),
                    cex=0.5,
                    txt_gp = fpTxtGp(ticks=gpar(cex=1)),
                    col=fpColors(box="royalblue",line="darkblue", summary="royalblue")
    ))
    while (!is.null(dev.list())) dev.off()
  }
}

xtick_u<-function(upper_limit) {
  ceiling(upper_limit)
}
xtick_l<-function(lower_limit) {
  floor(lower_limit)
}
##### Create graphs and write to file #####
for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    tabletext<-matrix(tabletext_b[[outcomes_b[i]]][[samples[j]]], nrow=9, ncol=5)
    cochrane_from_rmeta<-data.frame(cochrane_from_rmeta_b[[outcomes_b[i]]][[samples[j]]])
    xticks<-c(xtick_l(min(cochrane_from_rmeta$lower, na.rm = TRUE)), 0, xtick_u(max(cochrane_from_rmeta$upper, na.rm = TRUE)))
    pdf.options(reset = TRUE, onefile = FALSE)
    pdfname<-paste0("MVMR_",outcomes_b[i],"_",samples[j],".pdf")
    pdf(file=pdfname, width=15,height=15)
    plot(forestplot(tabletext,cochrane_from_rmeta,
                    new_page = TRUE,
                    is.summary=FALSE,
                    lineheight = unit(1, "cm"),
                    graphwidth=unit(10, "cm"),
                    boxsize=0.15,
                    xticks=xticks,
                    ci.vertices=TRUE,
                    colgap=unit(4, "mm"),
                    cex=0.5,
                    txt_gp = fpTxtGp(ticks=gpar(cex=1)),
                    col=fpColors(box="royalblue",line="darkblue", summary="royalblue")
    ))
    while (!is.null(dev.list())) dev.off()
  }
}
