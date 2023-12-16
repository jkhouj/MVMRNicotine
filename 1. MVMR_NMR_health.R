# This script was created by Jasmine Khouja 15.11.22. 

# The script conducts univariable and multivariable MR exploring the effects of
# nicotine and non-nicotine constituents of tobacco smoke (measured by nicotine
# metabolite ratio [NMR] and cigarettes per day [CPD]) on health outcomes 
# (chronic obstructive pulmonary disease [COPD], coronary heart disease [CHD], 
# forced expiratory volume [FEV], forced vital capacity [FVC], heart rate [HR],
# and body mass index [BMI])

# The SNPs used in this script were selected based on the finding from GSCAN 
# (Liu et al., 2019) and the cotinine consortium.

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

################################################################################
##### Set directory and memory #####
# Change as appropriate
################################################################################

wd<-"folder for working directory"
exp1_wd<-"folder which contains exposure 1 (NMR) data"
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

exp1_SNPlist<- read.table("NMR_SNPlist.txt", header=FALSE)
exp1_SNPlist<-rename(exp1_SNPlist, c("SNP" = "V1"))
exp2_SNPlist<- read.table("CPD_SNPlist.txt", header=FALSE)
exp2_SNPlist<-rename(exp2_SNPlist, c("SNP" = "V1"))

################################################################################
##### Removing rs117090198  ##### TO REMOVE FROM EXAMPLE
# Note: Removed because the SNP has a very high p-value prior to the conditional 
# independence analysis which is unusual and unexplained
################################################################################

exp1_SNPlist <- data.frame("SNP"= c(exp1_SNPlist[!exp1_SNPlist$SNP == "rs117090198", ]))

################################################################################

#####Extract exposure data for MR of exp1 and health outcomes #####
# manually check the column names in the exposure data set
# clump is false if the SNPs being inputted are (conditionally) independent 
# will need to clump (TRUE) if they are not
# sep is the separator "\t" is tab separated, " " is space separated, "," is csv 
# comma separated
# snp_col is the column with SNPs / RSIDs
# beta_col is the column where the effect size e.g., beta or odds ratio is
# se_col is the column where the standard error is
# eaf_col is the column where the effect allele frequency is
# effect_allele_col is the column where the effect or alternative allele is 
# (look at the data documentation)
# other_allele_col is the column where the other or reference allele is 
# pval_col is the column where the p value is
# samplesize_col is the column where the sample size is or this can be entered 
# manually
# min_pval selects the lowest p value that will be use - all values below this 
# will be entered as this value
# log_pval is true if the p-value is -log10(P) but the default is FALSE.

################################################################################
setwd(exp1_wd)
exp1_dat <- read_exposure_data(exp1_file,
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

setwd(wd)
write.csv(exp1_dat_mr,"exp1_dat_mr.csv", row.names = FALSE)

################################################################################
##### Extract exposure data for MR of exp2 and health outcomes #####
################################################################################

setwd(exp2_wd)
memory.limit(size = 80000)
exp2_dat <- read_exposure_data("CigarettesPerDay.WithoutUKB.txt",
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
write.csv(exp2_dat_mr,"exp2_dat_mr_noukb.csv", row.names = FALSE)

################################################################################
##### Extract outcome data for MR #####
# Set outcomes, samples (e.g., if different stratification of the data) 
# and sample sizes list for loops
# Extract outcome data
# Organise outcome names
# Adapt after paste0 to reflect the outcome folder and file name
# In this example, each outcome was within it's own folder 
# i.e., wd/outcome_sample_imptuted.txt/outcome_sample_imputed.txt
################################################################################

outcomes<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
samples<-c('current', 'ever', 'former', 'never')
samp_sizes <- c(49721, 213341, 163620, 258056) # current, ever, former, never
outcome_names<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')

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

mr_dfs_exp2 <- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mr_dfs_exp2[[outcomes[i]]][[samples[j]]] <- read_outcome_data(
      paste0(outcomes[i], "_", samples[j], "_imputed.txt/", outcomes[i], "_", samples[j], "_imputed.txt"),
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
##### Convert odds ratios to log odds if necessary #####
# not necessary, log odds already
################################################################################

################################################################################
##### Harmonising #####
################################################################################

harmonise_mr_dfs_exp1<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]] <- harmonise_data(
      exp1_dat_mr, 
      mr_dfs_exp1[[outcomes[i]]][[samples[j]]])
  }
}

harmonise_mr_dfs_exp2<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]] <- harmonise_data(
      exp2_dat_mr, 
      mr_dfs_exp2[[outcomes[i]]][[samples[j]]])
  }
}

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed_exp1<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    proxy_needed_exp1[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(
      exp1_SNPlist$SNP, 
      harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$SNP))
  }
}

proxy_needed_exp2<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    proxy_needed_exp2[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(
      exp2_SNPlist$SNP, 
      harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$SNP))
  }
}

################################################################################
##### Search for proxies if needed and reload data as required #####
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
    write.csv(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]], paste0("MR_exp2_noukb_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
  }
}

outcomes<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
samples<-c('current', 'ever', 'former', 'never')
samp_sizes <- c(49721, 213341, 163620, 258056) # current, ever, former, never
outcome_names<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')

################################################################################
# Create new variable with standardised continuous outcomes by dividing the beta 
# and se by the std dev 
################################################################################

outcomes_b<- c("BMI", "FEV", "FVC", "HR")
cont_out_sd_df<- read.csv("cont_outcome_sd.csv")
rownames(cont_out_sd_df)<-samples
harmonise_mr_dfs_exp1_std<-list()
harmonise_mr_dfs_exp2_std<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp1_std[[outcomes_b[i]]][[samples[j]]] <- harmonise_mr_dfs_exp1[[outcomes_b[i]]][[samples[j]]]
    harmonise_mr_dfs_exp1_std[[outcomes_b[i]]][[samples[j]]][["beta.outcome"]] <- 
     harmonise_mr_dfs_exp1[[outcomes_b[i]]][[samples[j]]][["beta.outcome"]]/cont_out_sd_df[samples[j],outcomes_b[i]]
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp1_std[[outcomes_b[i]]][[samples[j]]][["se.outcome"]] <- 
      harmonise_mr_dfs_exp1[[outcomes_b[i]]][[samples[j]]][["se.outcome"]]/cont_out_sd_df[samples[j],outcomes_b[i]]
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp2_std[[outcomes_b[i]]][[samples[j]]] <- harmonise_mr_dfs_exp2[[outcomes_b[i]]][[samples[j]]]
    harmonise_mr_dfs_exp2_std[[outcomes_b[i]]][[samples[j]]][["beta.outcome"]] <- 
      harmonise_mr_dfs_exp2[[outcomes_b[i]]][[samples[j]]][["beta.outcome"]]/cont_out_sd_df[samples[j],outcomes_b[i]]
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp2_std[[outcomes_b[i]]][[samples[j]]][["se.outcome"]] <- 
      harmonise_mr_dfs_exp2[[outcomes_b[i]]][[samples[j]]][["se.outcome"]]/cont_out_sd_df[samples[j],outcomes_b[i]]
  }
}

################################################################################
################################## MR ##########################################
################################################################################

################################################################################
##### Generate results inc. F and Q stats for heterogeneity #####
# Results are saved to result_exp[n] ivw_exp[n] and egger_exp_n
# IVW Q statistic saved to the first row in ptr_exp[n]
# Egger Q statistic saved to the second row in ptr_exp[n]
# F statistic saved to the third row in ptr_exp[n]
################################################################################

result_exp1<-list()
result_exp1_std<-list()
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

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    result_exp1_std[[outcomes_b[i]]][[samples[j]]] <- mr(
      harmonise_mr_dfs_exp1_std[[outcomes_b[i]]][[samples[j]]],
      parameters = default_parameters(),
      method_list = subset(mr_method_list(), use_by_default)$obj
)
result_exp1_std[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(result_exp1_std[[outcomes[i]]][[samples[j]]])
  }
}

result_exp2<-list()
result_exp2_std<-list()
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

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    harmonise_mr_dfs_exp2_std[[outcomes[i]]][[samples[j]]]$se.exposure<- as.numeric(harmonise_mr_dfs_exp2_std[[outcomes[i]]][[samples[j]]]$se.exposure)
    result_exp2_std[[outcomes_b[i]]][[samples[j]]] <- mr(
     harmonise_mr_dfs_exp2_std[[outcomes_b[i]]][[samples[j]]],
     parameters = default_parameters(),
     method_list = subset(mr_method_list(), use_by_default)$obj
)
result_exp2_std[[outcomes[i]]][[samples[j]]]<-generate_odds_ratios(result_exp2_std[[outcomes[i]]][[samples[j]]])
  }
}

################################################################################
##### Calculate I2GX #####
# Isq(y, s) where y = vector of effects and s = vector of standard errors
# Apply simulation extrapolation SIMEX corrections to MR-Egger analysis where I2GX estimates < 0.9 
# This would indicate the effect estimate is biased by 10% due to measurement error
# (Bowden, Del Greco, et al., 2016)
# (Lederer & Kichenhoff, 2006)
# Run SIMEX if below 0.9 for exp1 or exp2 - contact Jasmine for script
################################################################################

ISQ_exp1<-data.frame(matrix(NA, nrow = length(outcomes), ncol = length(samples)))

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ISQ_exp1[i,j]<-Isq(harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$beta.exposure, harmonise_mr_dfs_exp1[[outcomes[i]]][[samples[j]]]$se.exposure)
  }
}

ISQ_exp2<-data.frame(matrix(NA, nrow = length(outcomes), ncol = length(samples)))

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ISQ_exp2[i,j]<-Isq(harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$beta.exposure, harmonise_mr_dfs_exp2[[outcomes[i]]][[samples[j]]]$se.exposure)
  }
}

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

exp1_instr = exp1_dat_mr[,c("SNP","effect_allele.exposure","se.exposure","pval.exposure","beta.exposure","exposure","mr_keep.exposure","id.exposure","samplesize.exposure")]
exp2_instr = exp2_dat_mr[,c("SNP","effect_allele.exposure","se.exposure","pval.exposure","beta.exposure","exposure","mr_keep.exposure","id.exposure","samplesize.exposure")]

INSTR<- do.call("rbind", list(exp1_instr, exp2_instr))

################################################################################
##### Check for overlapping SNPs between the exposure instruments #####
# Do not need to drop these but good to be aware because it will impact the SNPs N
################################################################################

n_occur<- data.frame(table(INSTR$SNP))
n_occur[n_occur$Freq >1, ]
INSTR[INSTR$SNP %in% n_occur$Var1[n_occur$Freq >1], ]

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
##### Merge them to extract from outcome #####
################################################################################

exp1_dat_mvmr_1 = subset(exp1_dat_mvmr, select = c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "exposure", "mr_keep.exposure", "id.exposure", "eaf.exposure", "samplesize.exposure"))
exp2_dat_mvmr_1 = subset(exp2_dat_mvmr, select = c("SNP", "effect_allele.exposure", "other_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure", "exposure", "mr_keep.exposure", "id.exposure", "eaf.exposure", "samplesize.exposure"))

##### Check structure is the same #####
str(exp1_dat_mvmr_1)
str(exp2_dat_mvmr_1)

################################################################################
##### Merge ####
################################################################################

exposures<- do.call("rbind", list(exp1_dat_mvmr_1, exp2_dat_mvmr_1))

##### Save dataframes #####
setwd(wd)
write.csv(exposures,"exp1_exp2_noukb.csv", row.names = FALSE)

################################################################################
##### Find proxies missing from either the exp1 or exp2 dataset to add to missing outcome SNPs #####
# First highlights where a SNP is only present in one dataset, then identifies which
# Proxy files are split into a and b where some SNPs are missing from the exp2 dataset and some are missing from the exp1 dataset
# 1a = exp1 snps missing from the exp2 data 
# 1b = exp2 snps missing from the exp1 data 
# Either find proxies manually and replace in the 
################################################################################

n_occur<- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq <2, ]
exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq <2], ]

proxy_needed1 <- data.frame(exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq <2], ])
proxy_needed1a<- subset(proxy_needed1, id.exposure==id.exposure[2])
proxy_needed1b<- subset(proxy_needed1, id.exposure==id.exposure[1])
proxy_needed1 <- data.frame(proxy_needed1[1])
proxy_needed1a <- data.frame(proxy_needed1a[1])
proxy_needed1b <- data.frame(proxy_needed1b[1])

################################################################################
##### Search for proxies if needed and reload data as required #####
################################################################################

################################################################################
##### Clumping ####
# If one exposure (e.g., exp1) has many fewer associated SNPs:
# Change all p-values for exp1 to 1e-200 for clumping so that none are dropped
# Save original p-values first 
# Clump the data
# Add ID's back
# Revert all p-values for exp1 from 1e-200 back to original
# Split again to harmonise based on exposure id
################################################################################

exposures$oldpvalues <-exposures$pval.exposure

exposures<- exposures %>% 
  mutate(pval.exposure = if_else(exposures$SNP %in% exp1_instr$SNP, 1e-201, pval.exposure))

exposures$id.exposure[exposures$id.exposure == "2"] <- "1"
exposures<- clump_data(exposures, clump_kb=500, clump_r2=0.1) 
str(exposures) 

exposures$id.exposure[exposures$samplesize.exposure<6000] <- "1"
exposures$id.exposure[exposures$samplesize.exposure>6000] <- "2"

exposures$pval.exposure<-exposures$oldpvalues
exposures<-select(exposures,-c(oldpvalues))

# Added rs56113850 back in (removed due to LD with other NMR SNPs but is conditionally independent) #####
rs56113850 <- do.call("rbind", list(exp1_dat_mvmr_1, exp2_dat_mvmr_1))
rs56113850 <- rs56113850[rs56113850$SNP=="rs56113850", ]
exposures <- do.call("rbind", list(exposures, rs56113850))

exp1 = split(exposures, exposures$id.exposure)[['1']]
exp2 = split(exposures, exposures$id.exposure)[['2']]

################################################################################
##### Harmonise exp1 on exp2 #####
################################################################################

names(exp1) = gsub( "exposure", "outcome", names(exp1))
exp1_exp2 = harmonise_data(exp2, exp1)

################################################################################
##### Keep only snps that are present across both exposures ####
#Note: they would have frequency 1 if only available in one dataset
################################################################################

n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq == 2,]
exposures<- exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq == 2],]
str(exposures)

################################################################################
##### Format exposures #####
################################################################################

##### Keep only snps MrKeep= TRUE #####
exp1_exp2 = exp1_exp2[exp1_exp2$mr_keep== TRUE, ]
str(exp1_exp2) 

##### Split the tables - exp2 ##### 
exp2_H<- subset(exp1_exp2, id.exposure=="2", select= c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure))

##### Split the tables - exp1 ##### 
exp1_H<- subset(exp1_exp2, id.outcome=="1", select= c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome))

##### Turn exp1 from outcome to exposure to merge the datasets ##### 
names(exp1_H) <- gsub("outcome", "exposure", names(exp1_H))
Exposures_H<- merge(exp1_H, exp2_H, all= TRUE)
Exposures_H["Phenotype"]<- NA
Exposures_H$Phenotype[Exposures_H$id.exposure == 1] <- "exp1"
Exposures_H$Phenotype[Exposures_H$id.exposure == 2] <- "exp2"
str(Exposures_H)

################################################################################
##### Extract outcome data for MVMR #####
################################################################################

outcomes<-c('BMI', 'FEV', 'FVC', 'HR', 'CHD', 'COPD')
samples<-c('current', 'ever', 'former', 'never')
samp_sizes <- c(49721, 213341, 163620, 258056) # current, ever, former, never
setwd("")

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

################################################################################
##### Organise outcome #####
################################################################################

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dfs[[outcomes[i]]][[samples[j]]]["Phenotype"]<- NA
    mvmr_dfs[[outcomes[i]]][[samples[j]]]["Phenotype"]<- outcomes[i]
  }
}

################################################################################
##### Harmonise with outcome #####
################################################################################

mvmr_dat<-list()
for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dat[[outcomes[i]]][[samples[j]]]<- harmonise_data(Exposures_H, mvmr_dfs[[outcomes[i]]][[samples[j]]])
    mvmr_dat[[outcomes[i]]][[samples[j]]]<- mvmr_dat[[outcomes[i]]][[samples[j]]][mvmr_dat[[outcomes[i]]][[samples[j]]]$mr_keep== TRUE, ]
    str(mvmr_dat[[outcomes[i]]][[samples[j]]])
  }
}   

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed_mvmr<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    proxy_needed_mvmr[[outcomes[i]]][[samples[j]]]<- data.frame(setdiff(exp1_H$SNP, mvmr_dat[[outcomes[i]]][[samples[j]]]$SNP))
    
  }
}

################################################################################
##### Search for proxies if needed and reload data as required #####
################################################################################

################################################################################
##### Save dataframes #####
################################################################################

setwd(wd)

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    write.csv(mvmr_dat[[outcomes[i]]][[samples[j]]], paste0("MVMR_dat_", outcomes[i], "_", samples[j], ".csv"), row.names = FALSE)
  }
}

##### Can run from here if data frames are unchanged.##### 
#####Can skip if data is already loaded. #####
setwd(wd)
mvmr_dat <- list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    mvmr_dat[[outcomes[i]]][[samples[j]]]<-  read.csv(paste0("MVMR_dat_", outcomes[i], "_", samples[j], ".csv"))
  }
}

################################################################################
# Create new variable with standardised continuous outcomes by dividing the beta and se by
# the std dev 
################################################################################

outcomes_b<- c("BMI", "FEV", "FVC", "HR")
cont_out_sd_df<- read.csv("cont_outcome_sd.csv")
rownames(cont_out_sd_df)<-samples
mvmr_dat_std <- list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    mvmr_dat_std[[outcomes_b[i]]][[samples[j]]] <- mvmr_dat[[outcomes_b[i]]][[samples[j]]]
    mvmr_dat_std[[outcomes_b[i]]][[samples[j]]][["beta.outcome"]] <- 
      mvmr_dat[[outcomes_b[i]]][[samples[j]]][["beta.outcome"]]/cont_out_sd_df[samples[j],outcomes_b[i]]
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    mvmr_dat_std[[outcomes_b[i]]][[samples[j]]][["se.outcome"]] <- 
      mvmr_dat[[outcomes_b[i]]][[samples[j]]][["se.outcome"]]/cont_out_sd_df[samples[j],outcomes_b[i]]
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
    dfs[[outcomes[i]]][[samples[j]]]<- data.frame(bX1[[outcomes[i]]][[samples[j]]], 
                                                  bXse1[[outcomes[i]]][[samples[j]]], 
                                                  bX2[[outcomes[i]]][[samples[j]]], 
                                                  bXse2[[outcomes[i]]][[samples[j]]], 
                                                  bY[[outcomes[i]]][[samples[j]]], 
                                                  bYse[[outcomes[i]]][[samples[j]]], 
                                                  SNP[[outcomes[i]]][[samples[j]]])
    df_mvmr[[outcomes[i]]][[samples[j]]]<- format_mvmr(dfs[[outcomes[i]]][[samples[j]]][, c(1, 3)], 
                                                       dfs[[outcomes[i]]][[samples[j]]][, 5], 
                                                       dfs[[outcomes[i]]][[samples[j]]][, c(2, 4)], 
                                                       dfs[[outcomes[i]]][[samples[j]]][, 6], 
                                                       dfs[[outcomes[i]]][[samples[j]]][, 7])
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
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for NMR, and then          #####
##### orient the betas to be positive for CPD. In the paper, we           #####
##### report the result for each exposure only with the right orientation #####

clist<-c("bX2", "bY")

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    for (var in clist) {
      eval(parse(text=paste0(var, 
                             "[['", outcomes[i], "']]", 
                             "[['", samples[j], "']]", 
                             "<-ifelse(bX1", 
                             "[['", outcomes[i], "']]",
                             "[['", samples[j], "']]", 
                             ">0,",var, 
                             "[['", outcomes[i], "']]", 
                             "[['", samples[j], "']]", 
                             ",", var, 
                             "[['", outcomes[i], "']]", 
                             "[['", samples[j], "']]", 
                             "*-1)")))
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

##### Repeat for Standardised ##### 

bX1<-list()
bX2<-list()
bXse1<-list()
bXse2<-list()
bY<-list()
bYse<-list()
SNP<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    bX1[[outcomes_b[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$beta.exposure[mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$id.exposure== 1]))
    bX2[[outcomes_b[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$beta.exposure[mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$id.exposure== 2]))
    bXse1[[outcomes_b[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$se.exposure[mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$id.exposure== 1]))
    bXse2[[outcomes_b[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$se.exposure[mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$id.exposure== 2]))
    bY[[outcomes_b[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$beta.outcome[mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$id.exposure== 1]))
    bYse[[outcomes_b[i]]][[samples[j]]]<- as.numeric(c(mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$se.outcome[mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$id.exposure== 1]))
    SNP[[outcomes_b[i]]][[samples[j]]]<- c(mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$SNP[mvmr_dat_std[[outcomes_b[i]]][[samples[j]]]$id.exposure== 1])
    
  }
}

dfs_std<-list()
df_mvmr_std<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    dfs_std[[outcomes_b[i]]][[samples[j]]]<- data.frame(bX1[[outcomes_b[i]]][[samples[j]]], 
                                                        bXse1[[outcomes_b[i]]][[samples[j]]], 
                                                        bX2[[outcomes_b[i]]][[samples[j]]], 
                                                        bXse2[[outcomes_b[i]]][[samples[j]]], 
                                                        bY[[outcomes_b[i]]][[samples[j]]], 
                                                        bYse[[outcomes_b[i]]][[samples[j]]], 
                                                        SNP[[outcomes_b[i]]][[samples[j]]])
    df_mvmr_std[[outcomes_b[i]]][[samples[j]]]<- format_mvmr(dfs_std[[outcomes_b[i]]][[samples[j]]][, c(1, 3)], 
                                                             dfs_std[[outcomes_b[i]]][[samples[j]]][, 5], 
                                                             dfs_std[[outcomes_b[i]]][[samples[j]]][, c(2, 4)], 
                                                             dfs_std[[outcomes_b[i]]][[samples[j]]][, 6], 
                                                             dfs_std[[outcomes_b[i]]][[samples[j]]][, 7])
  }
}

set.seed(1234)
mod.MVMR_std<-list()
se_theta1MI.random_std<-list()
mod_std<-list()
mod_or_std<-list()
mod_or_ivw_std<-list()
res_std<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    mod.MVMR_std[[outcomes_b[i]]][[samples[j]]]<-lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]]-1, weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2)
    se_theta1MI.random_std[[outcomes_b[i]]][[samples[j]]] = summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]]-1, weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))$coef[1,2]/
      min(summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]]-1, weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))$sigma,1)
    
    mod_std[[outcomes_b[i]]][[samples[j]]]<- summary(mod.MVMR_std[[outcomes_b[i]]][[samples[j]]])
    
    mod_or_std[[outcomes_b[i]]][[samples[j]]] <- coef(summary(mod.MVMR_std[[outcomes_b[i]]][[samples[j]]]))
    colnames(mod_or_std[[outcomes_b[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]]<-as.data.frame(mod_or_std[[outcomes_b[i]]][[samples[j]]])
    mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]]<-generate_odds_ratios(mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]])
    res_std[[outcomes_b[i]]][[samples[j]]]<- ivw_mvmr(df_mvmr_std[[outcomes_b[i]]][[samples[j]]])
  }
}

##### Orientation exp1 #####

clist<-c("bX2", "bY")

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    for (var in clist) {
      eval(parse(text=paste0(var, 
                             "[['", outcomes_b[i], "']]", 
                             "[['", samples[j], "']]", 
                             "<-ifelse(bX1", 
                             "[['", outcomes_b[i], "']]",
                             "[['", samples[j], "']]", 
                             ">0,",var, 
                             "[['", outcomes_b[i], "']]", 
                             "[['", samples[j], "']]", 
                             ",", var, 
                             "[['", outcomes_b[i], "']]", 
                             "[['", samples[j], "']]", 
                             "*-1)")))
    }
    bX1[[outcomes_b[i]]][[samples[j]]]<-abs(bX1[[outcomes_b[i]]][[samples[j]]])
  }
}

##### MVMR Egger #####
mod.MVMRME_exp1_std<-list()
se_theta1ME.random_std<-list()
mod_ME_exp1_std<-list()
mod_ME_or_exp1_std<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    mod.MVMRME_exp1_std[[outcomes_b[i]]][[samples[j]]] <- summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]], weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))
    se_theta1ME.random_std[[outcomes_b[i]]][[samples[j]]] = summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]], weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))$coef[2,2]/
      min(summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]], weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))$sigma,1)
    mod_ME_exp1_std[[outcomes_b[i]]][[samples[j]]]<- summary(mod.MVMRME_exp1_std[[outcomes_b[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]] <- data.frame(mod.MVMRME_exp1_std[[outcomes_b[i]]][[samples[j]]][["coefficients"]])
    colnames(mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]]<-as.data.frame(mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]])
    mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]]<-generate_odds_ratios(mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]])
  }
}

##### Orientation exp2 #####

clist<-c("bX1","bY")
for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    for (var in clist) {
      eval(parse(text=paste0(var, "[['", outcomes_b[i], "']]", 
                             "[['", samples[j], "']]", 
                             "<-ifelse(bX2", 
                             "[['", outcomes_b[i], "']]", 
                             "[['", samples[j], "']]", 
                             ">0,",var, 
                             "[['", outcomes_b[i], "']]", 
                             "[['", samples[j], "']]", 
                             ",", var, 
                             "[['", outcomes_b[i], "']]", 
                             "[['", samples[j], "']]", 
                             "*-1)")))
    }
    bX2[[outcomes_b[i]]][[samples[j]]]<-abs(bX2[[outcomes_b[i]]][[samples[j]]])
  }
}

##### MVMR Egger #####
mod.MVMRME_exp2_std<-list()
se_theta1ME.random_std<-list()
mod_ME_exp2_std<-list()
mod_ME_or_exp2_std<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    mod.MVMRME_exp2_std[[outcomes_b[i]]][[samples[j]]] <- summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]], weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))
    se_theta1ME.random_std[[outcomes_b[i]]][[samples[j]]] = summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]], weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))$coef[2,2]/
      min(summary(lm(bY[[outcomes_b[i]]][[samples[j]]]~bX1[[outcomes_b[i]]][[samples[j]]]+bX2[[outcomes_b[i]]][[samples[j]]], weights=bYse[[outcomes_b[i]]][[samples[j]]]^-2))$sigma,1)
    mod_ME_exp2_std[[outcomes_b[i]]][[samples[j]]]<- summary(mod.MVMRME_exp2_std[[outcomes_b[i]]][[samples[j]]])
  }
}

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]] <- data.frame(mod.MVMRME_exp2_std[[outcomes_b[i]]][[samples[j]]][["coefficients"]])
    colnames(mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]]) <- c("b", "se", "t", "p")
    mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]]<-as.data.frame(mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]])
    mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]]<-generate_odds_ratios(mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]])
  }
}

################################################################################
##### Calculate F-statistic and covariance #####
# Note: >10 is strong
# correlation between exp1 and exp2 in the example = 0.458 (found in a paper - Buchwald)
# cov should be updated with the relevant genetic correlation between exp1 and exp2
# cov <- matrix(c(1,insert genetic correlation here,and here,1), nrow=2, ncol=2)
################################################################################

covmat<-list()
Fstat<-list()
cov <- matrix(c(1,0.458,0.458,1), nrow=2, ncol=2)

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    covmat[[outcomes[i]]][[samples[j]]] <- phenocov_mvmr(cov, df_mvmr[[outcomes[j]]][[samples[j]]][,c(6,7)] ) 
    Fstat[[outcomes[i]]][[samples[j]]] <- strength_mvmr(df_mvmr[[outcomes[j]]][[samples[j]]], gencov = covmat[[outcomes[i]]][[samples[j]]])
  }
}

################################################################################
##### Test for horizontal pleiotropy #####
# Note: Q should be greater than the number of SNPs included
################################################################################

ptr<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    ptr[[outcomes[i]]][[samples[j]]] <- pleiotropy_mvmr(df_mvmr[[outcomes[i]]][[samples[j]]], gencov = covmat[[outcomes[i]]][[samples[j]]])
  }
}

################################################################################
##### Q-statistic minimisation estimate #####
# qhet_mvmr is used to adjust for covariance but at present, CIs take substantial time to calculate and crash R
# Compare effects with and without adjustment
################################################################################

res_qhet<-list()

for (i in 1:length(outcomes)) {
  for (j in 1:length(samples)) {
    res_qhet[[outcomes[i]]][[samples[j]]] <- qhet_mvmr(df_mvmr[[outcomes[i]]][[samples[j]]], cov, CI = F, iterations = 100)
  }
}

################################################################################
##### Forest Plots #####
# Only plotting SIME because egger for 5e8 because the I2GX was below 0.6
# The results could be severely impacted by NOME violation and adjusting may be worse
# Update exp1 and exp2 to reflect the name of your exposures in this line e.g.:
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
    Exposure[[outcomes[i]]][[samples[j]]]<- c("exp1","exp1","exp1","exp1","exp2","exp2","exp2","exp2")
    
    Method[[outcomes[i]]][[samples[j]]]<- c("MR-IVW","MVMR-IVW","MR-Egger","MVMR-Egger","MR-IVW","MVMR-IVW","MR-Egger","MVMR-Egger")
    
    Beta[[outcomes[i]]][[samples[j]]] <- as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"b"],
                                                      mod_or_ivw[[outcomes[i]]][[samples[j]]] [1,"b"],
                                                      result_exp1[[outcomes[i]]][[samples[j]]] [1,"b"],
                                                      mod_ME_or_exp1[[outcomes[i]]][[samples[j]]] [2,"b"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]] [3,"b"],
                                                      mod_or_ivw[[outcomes[i]]][[samples[j]]] [2,"b"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]] [1,"b"],
                                                      mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"b"]))
    
    LCI[[outcomes[i]]][[samples[j]]] <- as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"lo_ci"],
                                                     mod_or_ivw[[outcomes[i]]][[samples[j]]] [1,"lo_ci"], 
                                                     result_exp1[[outcomes[i]]][[samples[j]]] [1,"lo_ci"],
                                                     mod_ME_or_exp1[[outcomes[i]]][[samples[j]]] [2,"lo_ci"],
                                                     result_exp2[[outcomes[i]]][[samples[j]]] [3,"lo_ci"],
                                                     mod_or_ivw[[outcomes[i]]][[samples[j]]] [2,"lo_ci"],
                                                     result_exp2[[outcomes[i]]][[samples[j]]] [1,"lo_ci"],
                                                     mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"lo_ci"]))
    
    UCI[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"up_ci"],
                                                    mod_or_ivw[[outcomes[i]]][[samples[j]]] [1,"up_ci"], 
                                                    result_exp1[[outcomes[i]]][[samples[j]]] [1,"up_ci"],
                                                    mod_ME_or_exp1[[outcomes[i]]][[samples[j]]] [2,"up_ci"],
                                                    result_exp2[[outcomes[i]]][[samples[j]]] [3,"up_ci"],
                                                    mod_or_ivw[[outcomes[i]]][[samples[j]]] [2,"up_ci"],
                                                    result_exp2[[outcomes[i]]][[samples[j]]] [1,"up_ci"],
                                                    mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"up_ci"]))  
    
    
    OR[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"or"],
                                                   mod_or_ivw[[outcomes[i]]][[samples[j]]][1,"or"],
                                                   result_exp1[[outcomes[i]]][[samples[j]]] [1,"or"],
                                                   mod_ME_or_exp1[[outcomes[i]]][[samples[j]]] [2,"or"],
                                                   result_exp2[[outcomes[i]]][[samples[j]]][3,"or"],
                                                   mod_or_ivw[[outcomes[i]]][[samples[j]]][2,"or"],
                                                   result_exp2[[outcomes[i]]][[samples[j]]][1,"or"],
                                                   mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"or"]))
    
    ORLCI[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"or_lci95"],
                                                      mod_or_ivw[[outcomes[i]]][[samples[j]]][1,"or_lci95"],
                                                      result_exp1[[outcomes[i]]][[samples[j]]][1,"or_lci95"],
                                                      mod_ME_or_exp1[[outcomes[i]]][[samples[j]]][2,"or_lci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][3,"or_lci95"],
                                                      mod_or_ivw[[outcomes[i]]][[samples[j]]][2,"or_lci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][1,"or_lci95"],
                                                      mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"or_lci95"])) 
    
    ORUCI[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"or_uci95"],
                                                      mod_or_ivw[[outcomes[i]]][[samples[j]]][1,"or_uci95"],
                                                      result_exp1[[outcomes[i]]][[samples[j]]][1,"or_uci95"],
                                                      mod_ME_or_exp1[[outcomes[i]]][[samples[j]]][2,"or_uci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][3,"or_uci95"],
                                                      mod_or_ivw[[outcomes[i]]][[samples[j]]][2,"or_uci95"],
                                                      result_exp2[[outcomes[i]]][[samples[j]]][1,"or_uci95"],
                                                      mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"or_uci95"]))
    
    SE[[outcomes[i]]][[samples[j]]]<-c(result_exp1[[outcomes[i]]][[samples[j]]][3,"se"],
                                       mod_or_ivw[[outcomes[i]]][[samples[j]]][1,"se"], 
                                       result_exp1[[outcomes[i]]][[samples[j]]][1,"se"],
                                       mod_ME_or_exp1[[outcomes[i]]][[samples[j]]][2,"se"],
                                       result_exp2[[outcomes[i]]][[samples[j]]][3,"se"],
                                       mod_or_ivw[[outcomes[i]]][[samples[j]]][2,"se"],
                                       result_exp2[[outcomes[i]]][[samples[j]]][1,"se"],
                                       mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"se"])
    
    p[[outcomes[i]]][[samples[j]]] <-as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"pval"],
                                                  mod_or_ivw[[outcomes[i]]][[samples[j]]][1,"p"],
                                                  result_exp1[[outcomes[i]]][[samples[j]]][1,"pval"],
                                                  mod_ME_or_exp1[[outcomes[i]]][[samples[j]]][2,"p"],
                                                  result_exp2[[outcomes[i]]][[samples[j]]][3,"pval"],
                                                  mod_or_ivw[[outcomes[i]]][[samples[j]]][2,"p"],
                                                  result_exp2[[outcomes[i]]][[samples[j]]][1,"pval"],
                                                  mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][3,"p"]))
    
    I2GX[[outcomes[i]]][[samples[j]]]<-c(".",".",ISQ_exp1[1,1],".",".",".",ISQ_exp2[1,1],".")
    
    Q[[outcomes[i]]][[samples[j]]]<-c(ptr_exp1[[outcomes[i]]][[samples[j]]][1,"Q"],
                                      ptr[[outcomes[i]]][[samples[j]]][["Qstat"]],
                                      ptr_exp1[[outcomes[i]]][[samples[j]]][2,"Q"],
                                      ptr[[outcomes[i]]][[samples[j]]][["Qstat"]],
                                      ptr_exp2[[outcomes[i]]][[samples[j]]][1,"Q"],
                                      ptr[[outcomes[i]]][[samples[j]]][["Qstat"]],
                                      ptr_exp2[[outcomes[i]]][[samples[j]]][2,"Q"],
                                      ptr[[outcomes[i]]][[samples[j]]][["Qstat"]])
    
    EggerI[[outcomes[i]]][[samples[j]]]<-c(".",".",egger_exp1[[outcomes[i]]][[samples[j]]][["b_i"]],
                                           mod_ME_or_exp1[[outcomes[i]]][[samples[j]]][1,"b"],
                                           ".",".",egger_exp2[[outcomes[i]]][[samples[j]]][["b_i"]],
                                           mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][1,"b"])
    
    EggerIp[[outcomes[i]]][[samples[j]]]<-c(".",".",egger_exp1[[outcomes[i]]][[samples[j]]][["pval_i"]],
                                            mod_ME_or_exp1[[outcomes[i]]][[samples[j]]][1,"p"],
                                            ".",".",egger_exp2[[outcomes[i]]][[samples[j]]][["pval_i"]],
                                            mod_ME_or_exp2[[outcomes[i]]][[samples[j]]][1,"p"])
    
    F_stat[[outcomes[i]]][[samples[j]]]<-c(ptr_exp1[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat[[outcomes[i]]][[samples[j]]][1,1],
                                           ptr_exp1[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat[[outcomes[i]]][[samples[j]]][1,1],
                                           ptr_exp2[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat[[outcomes[i]]][[samples[j]]][1,2],
                                           ptr_exp2[[outcomes[i]]][[samples[j]]][3,1],
                                           Fstat[[outcomes[i]]][[samples[j]]][1,2])
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


for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
AllRes[[outcomes_b[i]]][[samples[j]]]["Beta_std"]<- as.numeric(c(result_exp1_std[[outcomes_b[i]]][[samples[j]]][3,"b"],
                                                                  mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]] [1,"b"],
                                                                  result_exp1_std[[outcomes_b[i]]][[samples[j]]] [1,"b"],
                                                                  mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]] [2,"b"],
                                                                  result_exp2_std[[outcomes_b[i]]][[samples[j]]] [3,"b"],
                                                                  mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]] [2,"b"],
                                                                  result_exp2_std[[outcomes_b[i]]][[samples[j]]] [1,"b"],
                                                                  mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]][3,"b"]))
AllRes[[outcomes_b[i]]][[samples[j]]]["LCI_std"]<- as.numeric(c(result_exp1_std[[outcomes_b[i]]][[samples[j]]][3,"lo_ci"],
                                                                mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]] [1,"lo_ci"], 
                                                                result_exp1_std[[outcomes_b[i]]][[samples[j]]] [1,"lo_ci"],
                                                                mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]] [2,"lo_ci"],
                                                                result_exp2_std[[outcomes_b[i]]][[samples[j]]] [3,"lo_ci"],
                                                                mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]] [2,"lo_ci"],
                                                                result_exp2_std[[outcomes_b[i]]][[samples[j]]] [1,"lo_ci"],
                                                                mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]][3,"lo_ci"]))
AllRes[[outcomes_b[i]]][[samples[j]]]["UCI_std"]<-as.numeric(c(result_exp1_std[[outcomes_b[i]]][[samples[j]]][3,"up_ci"],
                                                               mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]] [1,"up_ci"], 
                                                               result_exp1_std[[outcomes_b[i]]][[samples[j]]] [1,"up_ci"],
                                                               mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]] [2,"up_ci"],
                                                               result_exp2_std[[outcomes_b[i]]][[samples[j]]] [3,"up_ci"],
                                                               mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]] [2,"up_ci"],
                                                               result_exp2_std[[outcomes_b[i]]][[samples[j]]] [1,"up_ci"],
                                                               mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]][3,"up_ci"]))  
AllRes[[outcomes_b[i]]][[samples[j]]]["SE_std"]<- c(result_exp1_std[[outcomes_b[i]]][[samples[j]]][3,"se"],
                                                    mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]][1,"se"], 
                                                    result_exp1_std[[outcomes_b[i]]][[samples[j]]][1,"se"],
                                                    mod_ME_or_exp1_std[[outcomes_b[i]]][[samples[j]]][2,"se"],
                                                    result_exp2_std[[outcomes_b[i]]][[samples[j]]][3,"se"],
                                                    mod_or_ivw_std[[outcomes_b[i]]][[samples[j]]][2,"se"],
                                                    result_exp2_std[[outcomes_b[i]]][[samples[j]]][1,"se"],
                                                    mod_ME_or_exp2_std[[outcomes_b[i]]][[samples[j]]][3,"se"])
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
    nSNPs[[outcomes[i]]][[samples[j]]]<-as.numeric(c(result_exp1[[outcomes[i]]][[samples[j]]][3,"nsnp"],
                                                     nrow(mvmr_dfs[[outcomes[i]]][[samples[j]]]),
                                                     result_exp1[[outcomes[i]]][[samples[j]]][1,"nsnp"],
                                                     length(mod.MVMRME_exp1[[outcomes[i]]][[samples[j]]][["residuals"]]),
                                                     result_exp2[[outcomes[i]]][[samples[j]]][3,"nsnp"],
                                                     nrow(mvmr_dfs[[outcomes[i]]][[samples[j]]]),
                                                     result_exp2[[outcomes[i]]][[samples[j]]][1,"nsnp"],
                                                     length(mod.MVMRME_exp1[[outcomes[i]]][[samples[j]]][["residuals"]])))
    
  }
}

##### Betas #####
# Create dataframe including the effect estimates, lower and upper confidence intervals
# For each outcome and sample, produce a  forest plot split for OR and beta
# insert names of results which are odds ratios versus betas in the appropriate place to create the appropriate forest plot

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


cochrane_from_rmeta_b_std<-list()
b_LCI_UCI_std<-list()
tabletext_b_std<-list()

for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    cochrane_from_rmeta_b_std[[outcomes_b[i]]][[samples[j]]]<- 
      structure(list(
        mean  = c(NA,AllRes[[outcomes_b[i]]][[samples[j]]]$Beta_std), 
        lower = c(NA,AllRes[[outcomes_b[i]]][[samples[j]]]$LCI_std),
        upper = c(NA,AllRes[[outcomes_b[i]]][[samples[j]]]$UCI_std)), 
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
    b_LCI_UCI[[outcomes_b[i]]][[samples[j]]]<-str_c(round(AllRes[[outcomes_b[i]]][[samples[j]]]$Beta, digits=2),
                                                     ", (",
                                                     round(AllRes[[outcomes_b[i]]][[samples[j]]]$LCI, digits=2),
                                                     ", ",
                                                     round(AllRes[[outcomes_b[i]]][[samples[j]]]$UCI, digits = 2),")")
    
  }
}


for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    b_LCI_UCI_std[[outcomes_b[i]]][[samples[j]]]<-str_c(round(AllRes[[outcomes_b[i]]][[samples[j]]]$Beta_std, digits=2),
                                                     ", (",
                                                     round(AllRes[[outcomes_b[i]]][[samples[j]]]$LCI_std, digits=2),
                                                     ", ",
                                                     round(AllRes[[outcomes_b[i]]][[samples[j]]]$UCI_std, digits = 2),")")
    
  }
}

##### Update exposure labels appropriately i.e., change NMR to your exp1 and Smoking Heaviness to your exp2

for (i in 1:length(outcomes_OR)) {
  for (j in 1:length(samples)) {
    tabletext_OR[[outcomes_OR[i]]][[samples[j]]]<-cbind(
      c("Exposure", "NMR", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
      c("N SNPs", nSNPs[[outcomes_OR[i]]][[samples[j]]]),
      c("Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"),
      c("OR (95% CI)", OR_LCI_UCI[[outcomes_OR[i]]][[samples[j]]]),
      c("P value", round(AllRes[[outcomes_OR[i]]][[samples[j]]]$p, digits=3)))
    tabletext_OR[[outcomes_OR[i]]][[samples[j]]][,5][tabletext_OR[[outcomes_OR[i]]][[samples[j]]][,5]==0]<- "<0.001"
  }
}


for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    tabletext_b[[outcomes_b[i]]][[samples[j]]]<-cbind(
      c("Exposure", "NMR", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
      c("N SNPs", nSNPs[[outcomes_b[i]]][[samples[j]]]),
      c("Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"),
      c("beta (95% CI)", b_LCI_UCI[[outcomes_b[i]]][[samples[j]]]),
      c("P value", round(AllRes[[outcomes_b[i]]][[samples[j]]]$p, digits=3)))
    tabletext_b[[outcomes_b[i]]][[samples[j]]][,5][tabletext_b[[outcomes_b[i]]][[samples[j]]][,5]==0]<- "<0.001"
  }
}



for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    tabletext_b_std[[outcomes_b[i]]][[samples[j]]]<-cbind(
      c("Exposure", "NMR", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
      c("N SNPs",  NA, NA, NA, NA, NA, NA, NA, NA),
      c("Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"),
      c("beta (95% CI)", b_LCI_UCI_std[[outcomes_b[i]]][[samples[j]]]),
      c("P value", round(AllRes[[outcomes_b[i]]][[samples[j]]]$p, digits=3)))
    tabletext_b_std[[outcomes_b[i]]][[samples[j]]][,5][tabletext_b_std[[outcomes_b[i]]][[samples[j]]][,5]==0]<- "<0.001"
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

##### Create graphs and write to file #####
for (i in 1:length(outcomes_b)) {
  for (j in 1:length(samples)) {
    tabletext_std<-matrix(tabletext_b_std[[outcomes_b[i]]][[samples[j]]], nrow=9, ncol=5)
    cochrane_from_rmeta<-data.frame(cochrane_from_rmeta_b_std[[outcomes_b[i]]][[samples[j]]])
    xticks<-c(xtick_l(min(cochrane_from_rmeta$lower, na.rm = TRUE)), 0, xtick_u(max(cochrane_from_rmeta$upper, na.rm = TRUE)))
    pdf.options(reset = TRUE, onefile = FALSE)
    pdfname<-paste0("MVMR_",outcomes_b[i],"_std_",samples[j],".pdf")
    pdf(file=pdfname, width=15,height=15)
    plot(forestplot(tabletext_std,cochrane_from_rmeta,
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
