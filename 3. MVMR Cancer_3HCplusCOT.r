# This script was created by Jasmine Khouja 25.05.22.

# The script conducts univariable and multivariable MR exploring the effects of
# nicotine and non-nicotine constituents of tobacco smoke (measured by 3HC plus
# cotinine [Cplus3HC] and cigarettes per day [CPD]) on cancer.

# The SNPs used in this script were selected based on the finding from GSCAN
# (Liu et al., 2019) and the Buchwald et al.

# The script calls a proxy searching script - contact Jasmine for access.

################################################################################
##### Load packages #####
################################################################################
library(usethis)
library(TwoSampleMR)
library(ieugwasr)
library(googleAuthR)
library(tidyverse)
library(stringr)
library(dplyr)
library(forestplot)
library(plyr)
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
################################################################################
setwd("")
memory.limit(size = 80000)

################################################################################
##### Set SNP lists for Cplus3HC, CPD and both #####
################################################################################

Cplus3HC_SNPlist <- read.table("COTPLUS3HC_SNPlist.txt", header = FALSE)
Cplus3HC_SNPlist <- rename(Cplus3HC_SNPlist, c("SNP" = "V1"))

CPD_SNPlist <- read.table("CPD_SNPlist.txt", header = FALSE)
CPD_SNPlist <- rename(CPD_SNPlist, c("SNP" = "V1"))

# SNPlist<- read.table("summarydata.txt", header=FALSE)
# SNPlist<-rename(SNPlist, c( "SNP"= "V1"))

################################################################################
##### Extract exposure data for MR of Cplus3HC and lung cancer #####
################################################################################

Cplus3HC_dat <- read_exposure_data("META_cotplus3hc_Extended.ma",
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

##### Creating a SNP list for Cplus3HC whereby p<5e-6 because only 3 #####
##### conditionally independent SNPS - 16 SNPs remain after clumping #####

Cplus3HC_SNPlist_5e6 <- Cplus3HC_dat[Cplus3HC_dat$pval.exposure <= 0.000005, ]
Cplus3HC_SNPlist_5e6 <- clump_data(
  Cplus3HC_SNPlist_5e6,
  clump_kb = 500,
  clump_r2 = 0.1,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)

Cplus3HC_dat_mr <- format_data(
  Cplus3HC_dat,
  type = "exposure",
  snps = Cplus3HC_SNPlist$SNP,
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

Cplus3HC_dat_mr_5e6 <- format_data(
  Cplus3HC_dat,
  type = "exposure",
  snps = Cplus3HC_SNPlist_5e6$SNP,
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

write.csv(
  Cplus3HC_dat_mr,
  "Cplus3HC_dat_mr.csv",
  row.names = FALSE
)
write.csv(
  Cplus3HC_dat_mr_5e6,
  "Cplus3HC_dat_mr_5e6.csv",
  row.names = FALSE
)
write.csv(
  Cplus3HC_dat,
  "Cplus3HC_dat.csv",
  row.names = FALSE
)

################################################################################
##### Extract exposure data for MR of CPD and lung cancer #####
################################################################################

setwd("")
memory.limit(size = 80000)
# Extract exposure data for MR of CPD and lung cancer
cpd_dat <- read_exposure_data("CigarettesPerDay.txt",
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

cpd_dat_mr <- format_data(
  cpd_dat,
  type = "exposure",
  snps = CPD_SNPlist$SNP,
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

write.csv(
  cpd_dat_mr,
  "cpd_dat_mr.csv",
  row.names = FALSE
)
write.csv(
  cpd_dat,
  "cpd_dat.csv",
  row.names = FALSE
)

################################################################################
##### Extract outcome data for MR #####
################################################################################

setwd("")
memory.limit(size = 80000)
outcome_dat_never_Cplus3HC <- read_outcome_data(
  "Never.txt",
  snps = Cplus3HC_dat_mr$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat_ever_Cplus3HC <- read_outcome_data(
  "Ever.txt",
  snps = Cplus3HC_dat_mr$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 40187,
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd("")
memory.limit(size = 80000)
outcome_dat_never_Cplus3HC_5e6 <- read_outcome_data(
  "Never.txt",
  snps = Cplus3HC_dat_mr_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = "N_study",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat_ever_Cplus3HC_5e6 <- read_outcome_data(
  "Ever.txt",
  snps = Cplus3HC_dat_mr_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = "N_study",
  min_pval = 1e-200,
  log_pval = FALSE
)


outcome_dat_never_cpd <- read_outcome_data(
  "Never.txt",
  snps = cpd_dat_mr$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat_ever_cpd <- read_outcome_data(
  "Ever.txt",
  snps = cpd_dat_mr$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 40187,
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Convert odds ratios to log odds #####
# Note: the exposures are already on the log scale.
################################################################################

outcome_dat_never_Cplus3HC$beta.outcome <- as.numeric(as.character(
  outcome_dat_never_Cplus3HC$beta.outcome
))
outcome_dat_ever_Cplus3HC$beta.outcome <- as.numeric(as.character(
  outcome_dat_ever_Cplus3HC$beta.outcome
))
outcome_dat_never_Cplus3HC_5e6$beta.outcome <- as.numeric(as.character(
  outcome_dat_never_Cplus3HC_5e6$beta.outcome
))
outcome_dat_ever_Cplus3HC_5e6$beta.outcome <- as.numeric(as.character(
  outcome_dat_ever_Cplus3HC_5e6$beta.outcome
))
outcome_dat_never_cpd$beta.outcome <- as.numeric(as.character(
  outcome_dat_never_cpd$beta.outcome
))
outcome_dat_ever_cpd$beta.outcome <- as.numeric(as.character(
  outcome_dat_ever_cpd$beta.outcome
))
outcome_dat_never_Cplus3HC["beta.outcome"] <- log(
  outcome_dat_never_Cplus3HC["beta.outcome"]
)
outcome_dat_ever_Cplus3HC["beta.outcome"] <- log(
  outcome_dat_ever_Cplus3HC["beta.outcome"]
)
outcome_dat_never_Cplus3HC_5e6["beta.outcome"] <- log(
  outcome_dat_never_Cplus3HC_5e6["beta.outcome"]
)
outcome_dat_ever_Cplus3HC_5e6["beta.outcome"] <- log(
  outcome_dat_ever_Cplus3HC_5e6["beta.outcome"]
)
outcome_dat_never_cpd["beta.outcome"] <- log(
  outcome_dat_never_cpd["beta.outcome"]
)
outcome_dat_ever_cpd["beta.outcome"] <- log(
  outcome_dat_ever_cpd["beta.outcome"]
)

# Outcome se is already se for logodds rather than odds:
# Checked using
# tail <- 2
# se <- abs(outcome_dat_never_Cplus3HC$beta.outcome[1]/ 
# qnorm(outcome_dat_never_Cplus3HC$pval.outcome[1]/tail))
# outcome_dat_never_Cplus3HC$se[1]

################################################################################
##### Harmonising #####
################################################################################

dat_Cplus3HC_never <- harmonise_data(
  Cplus3HC_dat_mr,
  outcome_dat_never_Cplus3HC,
  action = 2
)
dat_Cplus3HC_ever <- harmonise_data(
  Cplus3HC_dat_mr,
  outcome_dat_ever_Cplus3HC,
  action = 2
)
dat_Cplus3HC_never_5e6 <- harmonise_data(
  Cplus3HC_dat_mr_5e6,
  outcome_dat_never_Cplus3HC_5e6,
  action = 2
)
dat_Cplus3HC_ever_5e6 <- harmonise_data(
  Cplus3HC_dat_mr_5e6,
  outcome_dat_ever_Cplus3HC_5e6,
  action = 2
)
dat_cpd_never <- harmonise_data(
  cpd_dat_mr,
  outcome_dat_never_cpd,
  action = 2
)
dat_cpd_ever <- harmonise_data(
  cpd_dat_mr,
  outcome_dat_ever_cpd,
  action = 2
)

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed1 <- data.frame(setdiff(
  Cplus3HC_SNPlist$SNP,
  dat_Cplus3HC_never$SNP
))
proxy_needed2 <- data.frame(setdiff(
  Cplus3HC_SNPlist$SNP,
  dat_Cplus3HC_ever$SNP
))
proxy_needed3 <- data.frame(setdiff(
  Cplus3HC_SNPlist_5e6$SNP,
  dat_Cplus3HC_never_5e6$SNP
))
proxy_needed4 <- data.frame(setdiff(
  Cplus3HC_SNPlist_5e6$SNP,
  dat_Cplus3HC_ever_5e6$SNP
))
proxy_needed5 <- data.frame(setdiff(
  CPD_SNPlist$SNP,
  dat_cpd_never$SNP
))
proxy_needed6 <- data.frame(setdiff(
  CPD_SNPlist$SNP,
  dat_cpd_ever$SNP
))

################################################################################
##### Run to find proxies #####
################################################################################

source("proxy_search_loop_CancerMVMR_Cplus3HC.R")

################################################################################
##### Re-read outcome data with proxies #####
# Note: SNPs that were missing in the outcome have now been replaced with proxy 
# SNPs in the exposure
# Need to reload the outcome data looking for the proxy SNPs
################################################################################

setwd("")
memory.limit(size = 80000)
outcome_dat_never_Cplus3HC <- read_outcome_data(
  "Never.txt",
  snps = Cplus3HC_dat_mr_never$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = "N_study",
  min_pval = 1e-200,
  log_pval = FALSE
)
setwd("")
memory.limit(size = 80000)
outcome_dat_ever_Cplus3HC <- read_outcome_data(
  "Ever.txt",
  snps = Cplus3HC_dat_mr_ever$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = "N_study",
  min_pval = 1e-200,
  log_pval = FALSE
)
setwd("")
memory.limit(size = 80000)
outcome_dat_never_Cplus3HC_5e6 <- read_outcome_data(
  "Never.txt",
  snps = Cplus3HC_dat_mr_never_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = "N_study",
  min_pval = 1e-200,
  log_pval = FALSE
)
setwd("")
memory.limit(size = 80000)
outcome_dat_ever_Cplus3HC_5e6 <- read_outcome_data(
  "Ever.txt",
  snps = Cplus3HC_dat_mr_ever_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = "N_study",
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat_never_cpd <- read_outcome_data(
  "Never.txt",
  snps = cpd_dat_mr_never$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

outcome_dat_ever_cpd <- read_outcome_data(
  "Ever.txt",
  snps = cpd_dat_mr_ever$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 40187,
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Convert odds ratios to log odds again #####
################################################################################

outcome_dat_never_Cplus3HC$beta.outcome <- as.numeric(
  as.character(outcome_dat_never_Cplus3HC$beta.outcome)
)
outcome_dat_ever_Cplus3HC$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever_Cplus3HC$beta.outcome)
)
outcome_dat_never_Cplus3HC_5e6$beta.outcome <- as.numeric(
  as.character(outcome_dat_never_Cplus3HC_5e6$beta.outcome)
)
outcome_dat_ever_Cplus3HC_5e6$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever_Cplus3HC_5e6$beta.outcome)
)
outcome_dat_never_cpd$beta.outcome <- as.numeric(
  as.character(outcome_dat_never_cpd$beta.outcome)
)
outcome_dat_ever_cpd$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever_cpd$beta.outcome)
)

outcome_dat_never_Cplus3HC["beta.outcome"] <- log(
  outcome_dat_never_Cplus3HC["beta.outcome"]
)
outcome_dat_ever_Cplus3HC["beta.outcome"] <- log(
  outcome_dat_ever_Cplus3HC["beta.outcome"]
)
outcome_dat_never_Cplus3HC_5e6["beta.outcome"] <- log(
  outcome_dat_never_Cplus3HC_5e6["beta.outcome"]
)
outcome_dat_ever_Cplus3HC_5e6["beta.outcome"] <- log(
  outcome_dat_ever_Cplus3HC_5e6["beta.outcome"]
)
outcome_dat_never_cpd["beta.outcome"] <- log(
  outcome_dat_never_cpd["beta.outcome"]
)
outcome_dat_ever_cpd["beta.outcome"] <- log(
  outcome_dat_ever_cpd["beta.outcome"]
)

################################################################################
##### Harmonising data after proxies found #####
################################################################################

dat_Cplus3HC_never <- harmonise_data(
  Cplus3HC_dat_mr_never,
  outcome_dat_never_Cplus3HC,
  action = 2
)
dat_Cplus3HC_ever <- harmonise_data(
  Cplus3HC_dat_mr_ever,
  outcome_dat_ever_Cplus3HC,
  action = 2
)
dat_Cplus3HC_never_5e6 <- harmonise_data(
  Cplus3HC_dat_mr_never_5e6,
  outcome_dat_never_Cplus3HC_5e6,
  action = 2
)
dat_Cplus3HC_ever_5e6 <- harmonise_data(
  Cplus3HC_dat_mr_ever_5e6,
  outcome_dat_ever_Cplus3HC_5e6,
  action = 2
)
dat_cpd_never <- harmonise_data(
  cpd_dat_mr_never,
  outcome_dat_never_cpd,
  action = 2
)
dat_cpd_ever <- harmonise_data(
  cpd_dat_mr_ever,
  outcome_dat_ever_cpd,
  action = 2
)

################################################################################
##### Harmonising negative effect alleles #####
# Note: Must ensure the beta is positive for the MR results to be correct
# This means you need to flip the beta and alleles so that the effect is positive
################################################################################

for (i in 1:length(dat_Cplus3HC_never$beta.exposure)) {
  if (dat_Cplus3HC_never$beta.exposure[i] < 0) {
    dat_Cplus3HC_never$beta.outcome[i] <- -1 * dat_Cplus3HC_never$beta.outcome[i]
  }
  if (dat_Cplus3HC_never$beta.exposure[i] < 0) {
    dat_Cplus3HC_never$beta.exposure[i] <- -1 * dat_Cplus3HC_never$beta.exposure[i]
  }
}
for (i in 1:length(dat_Cplus3HC_ever$beta.exposure)) {
  if (dat_Cplus3HC_ever$beta.exposure[i] < 0) {
    dat_Cplus3HC_ever$beta.outcome[i] <- -1 * dat_Cplus3HC_ever$beta.outcome[i]
  }
  if (dat_Cplus3HC_ever$beta.exposure[i] < 0) {
    dat_Cplus3HC_ever$beta.exposure[i] <- -1 * dat_Cplus3HC_ever$beta.exposure[i]
  }
}
for (i in 1:length(dat_Cplus3HC_never_5e6$beta.exposure)) {
  if (dat_Cplus3HC_never_5e6$beta.exposure[i] < 0) {
    dat_Cplus3HC_never_5e6$beta.outcome[i] <- -1 * dat_Cplus3HC_never_5e6$beta.outcome[i]
  }
  if (dat_Cplus3HC_never_5e6$beta.exposure[i] < 0) {
    dat_Cplus3HC_never_5e6$beta.exposure[i] <- -1 * dat_Cplus3HC_never_5e6$beta.exposure[i]
  }
}
for (i in 1:length(dat_Cplus3HC_ever_5e6$beta.exposure)) {
  if (dat_Cplus3HC_ever_5e6$beta.exposure[i] < 0) {
    dat_Cplus3HC_ever_5e6$beta.outcome[i] <- -1 * dat_Cplus3HC_ever_5e6$beta.outcome[i]
  }
  if (dat_Cplus3HC_ever_5e6$beta.exposure[i] < 0) {
    dat_Cplus3HC_ever_5e6$beta.exposure[i] <- -1 * dat_Cplus3HC_ever_5e6$beta.exposure[i]
  }
}
for (i in 1:length(dat_cpd_never$beta.exposure)) {
  if (dat_cpd_never$beta.exposure[i] < 0) {
    dat_cpd_never$beta.outcome[i] <- -1 * dat_cpd_never$beta.outcome[i]
  }
  if (dat_cpd_never$beta.exposure[i] < 0) {
    dat_cpd_never$beta.exposure[i] <- -1 * dat_cpd_never$beta.exposure[i]
  }
}
for (i in 1:length(dat_cpd_ever$beta.exposure)) {
  if (dat_cpd_ever$beta.exposure[i] < 0) {
    dat_cpd_ever$beta.outcome[i] <- -1 * dat_cpd_ever$beta.outcome[i]
  }
  if (dat_cpd_ever$beta.exposure[i] < 0) {
    dat_cpd_ever$beta.exposure[i] <- -1 * dat_cpd_ever$beta.exposure[i]
  }
}

################################################################################
##### Save dataframes #####
# Note: This is not necessary but can be useful for sharing the datasets
# Also can avoid rerunning time-consuming SNP extraction steps above
################################################################################

write.csv(dat_Cplus3HC_never, "dat_Cplus3HC_never.csv", row.names = FALSE)
write.csv(dat_Cplus3HC_ever, "dat_Cplus3HC_ever.csv", row.names = FALSE)
write.csv(dat_Cplus3HC_never_5e6, "dat_Cplus3HC_never_5e6.csv", row.names = FALSE)
write.csv(dat_Cplus3HC_ever_5e6, "dat_Cplus3HC_ever_5e6.csv", row.names = FALSE)
write.csv(dat_cpd_never, "dat_cpd_never.csv", row.names = FALSE)
write.csv(dat_cpd_ever, "dat_cpd_ever.csv", row.names = FALSE)

################################################################################
################################## MR ##########################################
################################################################################

################################################################################
##### Generate results inc. F and Q stats for heterogeneity for: #####
# Cplus3HC never
# Cplus3HC ever
# Cplus3HC never(5e6)
# Cplus3HC ever (5e6)
# CPD never
# CPD ever
################################################################################

# Cplus3HC never
result_Cplus3HC_never <- mr(
  dat_Cplus3HC_never,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_Cplus3HC_never <- generate_odds_ratios(result_Cplus3HC_never)

dat_Cplus3HC_never <- subset(dat_Cplus3HC_never, mr_keep)
mr_n_Cplus3HC <- TwoSampleMR::mr_ivw(
  dat_Cplus3HC_never$beta.exposure,
  dat_Cplus3HC_never$beta.outcome,
  dat_Cplus3HC_never$se.exposure,
  dat_Cplus3HC_never$se.outcome,
  parameters = default_parameters()
)
ptr_n_Cplus3HC <- data.frame(mr_n_Cplus3HC["Q"])
egger_n_Cplus3HC <- mr_egger_regression(
  dat_Cplus3HC_never$beta.exposure,
  dat_Cplus3HC_never$beta.outcome,
  dat_Cplus3HC_never$se.exposure,
  dat_Cplus3HC_never$se.outcome,
  parameters
)
ptr_n_Cplus3HC[2, 1] <- egger_n_Cplus3HC["Q"]
F <- dat_Cplus3HC_never$beta.exposure^2 / dat_Cplus3HC_never$se.exposure^2
mF <- mean(F)
ptr_n_Cplus3HC[3, 1] <- mF

# Cplus3HC ever
result_Cplus3HC_ever <- mr(
  dat_Cplus3HC_ever,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_Cplus3HC_ever <- generate_odds_ratios(result_Cplus3HC_ever)

dat_Cplus3HC_ever <- subset(dat_Cplus3HC_ever, mr_keep)
mr_e_Cplus3HC <- TwoSampleMR::mr_ivw(
  dat_Cplus3HC_ever$beta.exposure,
  dat_Cplus3HC_ever$beta.outcome,
  dat_Cplus3HC_ever$se.exposure,
  dat_Cplus3HC_ever$se.outcome,
  parameters = default_parameters()
)
ptr_e_Cplus3HC <- data.frame(mr_e_Cplus3HC["Q"])
egger_e_Cplus3HC <- mr_egger_regression(
  dat_Cplus3HC_ever$beta.exposure,
  dat_Cplus3HC_ever$beta.outcome,
  dat_Cplus3HC_ever$se.exposure,
  dat_Cplus3HC_ever$se.outcome,
  parameters
)
ptr_e_Cplus3HC[2, 1] <- egger_e_Cplus3HC["Q"]
F <- dat_Cplus3HC_ever$beta.exposure^2 / dat_Cplus3HC_ever$se.exposure^2
mF <- mean(F)
ptr_e_Cplus3HC[3, 1] <- mF

# Cplus3HC never(5e6)
result_Cplus3HC_never_5e6 <- mr(
  dat_Cplus3HC_never_5e6,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_Cplus3HC_never_5e6 <- generate_odds_ratios(result_Cplus3HC_never_5e6)

dat_Cplus3HC_never_5e6 <- subset(dat_Cplus3HC_never_5e6, mr_keep)
mr_n_Cplus3HC_5e6 <- TwoSampleMR::mr_ivw(
  dat_Cplus3HC_never_5e6$beta.exposure,
  dat_Cplus3HC_never_5e6$beta.outcome,
  dat_Cplus3HC_never_5e6$se.exposure,
  dat_Cplus3HC_never_5e6$se.outcome,
  parameters = default_parameters()
)
ptr_n_Cplus3HC_5e6 <- data.frame(mr_n_Cplus3HC_5e6["Q"])
egger_n_Cplus3HC_5e6 <- mr_egger_regression(
  dat_Cplus3HC_never_5e6$beta.exposure,
  dat_Cplus3HC_never_5e6$beta.outcome,
  dat_Cplus3HC_never_5e6$se.exposure,
  dat_Cplus3HC_never_5e6$se.outcome,
  parameters
)
ptr_n_Cplus3HC_5e6[2, 1] <- egger_n_Cplus3HC_5e6["Q"]
F <- dat_Cplus3HC_never_5e6$beta.exposure^2 / dat_Cplus3HC_never_5e6$se.exposure^2
mF <- mean(F)
ptr_n_Cplus3HC_5e6[3, 1] <- mF

# Cplus3HC ever (5e6)
result_Cplus3HC_ever_5e6 <- mr(
  dat_Cplus3HC_ever_5e6,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_Cplus3HC_ever_5e6 <- generate_odds_ratios(result_Cplus3HC_ever_5e6)

dat_Cplus3HC_ever_5e6 <- subset(dat_Cplus3HC_ever_5e6, mr_keep)
mr_e_Cplus3HC_5e6 <- TwoSampleMR::mr_ivw(
  dat_Cplus3HC_ever_5e6$beta.exposure, 
  dat_Cplus3HC_ever_5e6$beta.outcome, 
  dat_Cplus3HC_ever_5e6$se.exposure, 
  dat_Cplus3HC_ever_5e6$se.outcome,
  parameters = default_parameters()
  )
ptr_e_Cplus3HC_5e6 <- data.frame(mr_e_Cplus3HC_5e6["Q"])
egger_e_Cplus3HC_5e6 <- mr_egger_regression(
  dat_Cplus3HC_ever_5e6$beta.exposure, 
  dat_Cplus3HC_ever_5e6$beta.outcome, 
  dat_Cplus3HC_ever_5e6$se.exposure, 
  dat_Cplus3HC_ever_5e6$se.outcome,
  parameters
  )
ptr_e_Cplus3HC_5e6[2, 1] <- egger_e_Cplus3HC_5e6["Q"]
F <- dat_Cplus3HC_ever_5e6$beta.exposure^2 / dat_Cplus3HC_ever_5e6$se.exposure^2
mF <- mean(F)
ptr_e_Cplus3HC_5e6[3, 1] <- mF

# CPD never
result_cpd_never <- mr(
  dat_cpd_never,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)

result_cpd_never <- generate_odds_ratios(result_cpd_never)

dat_cpd_never <- subset(dat_cpd_never, mr_keep)
mr_n_cpd <- TwoSampleMR::mr_ivw(
  dat_cpd_never$beta.exposure,
  dat_cpd_never$beta.outcome,
  dat_cpd_never$se.exposure,
  dat_cpd_never$se.outcome,
  parameters = default_parameters()
)
ptr_n_cpd <- data.frame(mr_n_cpd["Q"])
egger_n_cpd <- mr_egger_regression(
  dat_cpd_never$beta.exposure,
  dat_cpd_never$beta.outcome,
  dat_cpd_never$se.exposure,
  dat_cpd_never$se.outcome,
  parameters
)
ptr_n_cpd[2, 1] <- egger_n_cpd["Q"]
F <- dat_cpd_never$beta.exposure^2 / dat_cpd_never$se.exposure^2
mF <- mean(F)
ptr_n_cpd[3, 1] <- mF

# CPD ever
result_cpd_ever <- mr(
  dat_cpd_ever,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_cpd_ever <- generate_odds_ratios(result_cpd_ever)

dat_cpd_ever <- subset(dat_cpd_ever, mr_keep)
mr_e_cpd <- TwoSampleMR::mr_ivw(
  dat_cpd_ever$beta.exposure,
  dat_cpd_ever$beta.outcome,
  dat_cpd_ever$se.exposure,
  dat_cpd_ever$se.outcome,
  parameters = default_parameters()
)
ptr_e_cpd <- data.frame(mr_e_cpd["Q"])
egger_e_cpd <- mr_egger_regression(
  dat_cpd_ever$beta.exposure,
  dat_cpd_ever$beta.outcome,
  dat_cpd_ever$se.exposure,
  dat_cpd_ever$se.outcome,
  parameters
)
ptr_e_cpd[2, 1] <- egger_e_cpd["Q"]
F <- dat_cpd_ever$beta.exposure^2 / dat_cpd_ever$se.exposure^2
mF <- mean(F)
ptr_e_cpd[3, 1] <- mF

################################################################################
##### Calculate I2GX #####
# Isq(y, s) where y = vector of effects and s = vector of standard errors
# Apply simulation extrapolation SIMEX corrections to MR-Egger analysis where 
# I2GX estimates < 0.9
# This would indicate the effect estimate is biased by 10% due to measurement 
# error
# (Bowden, Del Greco, et al., 2016)
# (Lederer & K?chenhoff, 2006)
################################################################################

ISQ <- data.frame(c(1:6))
ISQ[1, 1] <- Isq(
  dat_Cplus3HC_never$beta.exposure,
  dat_Cplus3HC_never$se.exposure
)

ISQ[2, 1] <- Isq(
  dat_Cplus3HC_ever$beta.exposure,
  dat_Cplus3HC_ever$se.exposure
)

ISQ[3, 1] <- Isq(
  dat_Cplus3HC_never_5e6$beta.exposure,
  dat_Cplus3HC_never_5e6$se.exposure
)

ISQ[4, 1] <- Isq(
  dat_Cplus3HC_ever_5e6$beta.exposure,
  dat_Cplus3HC_ever_5e6$se.exposure
)

ISQ[5, 1] <- Isq(
  dat_cpd_never$beta.exposure,
  dat_cpd_never$se.exposure
)

ISQ[6, 1] <- Isq(
  dat_cpd_ever$beta.exposure,
  dat_cpd_ever$se.exposure
)

source("SIMEX_Cplus3HC.R")

################################################################################
##### Forest plot results #####
##### Note: doesn't include simex corrections
################################################################################

forest_Cplus3HC_never <- mr_forest(
  mr_input(
    bx = dat_Cplus3HC_never$beta.exposure,
    bxse = dat_Cplus3HC_never$se.exposure,
    by = dat_Cplus3HC_never$beta.outcome,
    byse = dat_Cplus3HC_never$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)
forest_Cplus3HC_ever <- mr_forest(
  mr_input(
    bx = dat_Cplus3HC_ever$beta.exposure,
    bxse = dat_Cplus3HC_ever$se.exposure,
    by = dat_Cplus3HC_ever$beta.outcome,
    byse = dat_Cplus3HC_ever$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)
forest <- mr_forest(mr_input(
  bx = dat_Cplus3HC_ever$beta.exposure,
  bxse = dat_Cplus3HC_ever$se.exposure,
  by = dat_Cplus3HC_ever$beta.outcome,
  byse = dat_Cplus3HC_ever$se.outcome
))
forest2 <- forest + coord_cartesian(xlim = c(-2, 2))
forest2

forest_Cplus3HC_never_5e6 <- mr_forest(
  mr_input(
    bx = dat_Cplus3HC_never_5e6$beta.exposure,
    bxse = dat_Cplus3HC_never_5e6$se.exposure,
    by = dat_Cplus3HC_never_5e6$beta.outcome,
    byse = dat_Cplus3HC_never_5e6$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)
forest_Cplus3HC_ever_5e6 <- mr_forest(
  mr_input(
    bx = dat_Cplus3HC_ever_5e6$beta.exposure,
    bxse = dat_Cplus3HC_ever_5e6$se.exposure,
    by = dat_Cplus3HC_ever_5e6$beta.outcome,
    byse = dat_Cplus3HC_ever_5e6$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)
forest_5e6 <- mr_forest(mr_input(
  bx = dat_Cplus3HC_ever_5e6$beta.exposure,
  bxse = dat_Cplus3HC_ever_5e6$se.exposure,
  by = dat_Cplus3HC_ever_5e6$beta.outcome,
  byse = dat_Cplus3HC_ever_5e6$se.outcome
))
forest2_5e6 <- forest_5e6 + coord_cartesian(xlim = c(-2, 2))
forest2_5e6

forest_cpd_never <- mr_forest(
  mr_input(
    bx = dat_cpd_never$beta.exposure,
    bxse = dat_cpd_never$se.exposure,
    by = dat_cpd_never$beta.outcome,
    byse = dat_cpd_never$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)
forest_cpd_ever <- mr_forest(
  mr_input(
    bx = dat_cpd_ever$beta.exposure,
    bxse = dat_cpd_ever$se.exposure,
    by = dat_cpd_ever$beta.outcome,
    byse = dat_cpd_ever$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)

forest_cpd <- mr_forest(mr_input(
  bx = dat_cpd_ever$beta.exposure,
  bxse = dat_cpd_ever$se.exposure,
  by = dat_cpd_ever$beta.outcome,
  byse = dat_cpd_ever$se.outcome
))
forest2_cpd <- forest_cpd + coord_cartesian(xlim = c(-2, 2))
forest2_cpd

mr_funnel(mr_input(
  bx = dat_Cplus3HC_ever$beta.exposure,
  bxse = dat_Cplus3HC_ever$se.exposure,
  by = dat_Cplus3HC_ever$beta.outcome,
  byse = dat_Cplus3HC_ever$se.outcome
))
mr_funnel(mr_input(
  bx = dat_Cplus3HC_ever_5e6$beta.exposure,
  bxse = dat_Cplus3HC_ever_5e6$se.exposure,
  by = dat_Cplus3HC_ever_5e6$beta.outcome,
  byse = dat_Cplus3HC_ever_5e6$se.outcome
))
mr_funnel(mr_input(
  bx = dat_cpd_ever$beta.exposure,
  bxse = dat_cpd_ever$se.exposure,
  by = dat_cpd_ever$beta.outcome,
  byse = dat_cpd_ever$se.outcome
))

################################################################################
################################   MVMR  #######################################
################################################################################

################################################################################
###### Merge exposure SNP lists #####
# Replace rs114612145 with rs77107237 to reflect RSID change
# Load data if not already loaded
################################################################################

Cplus3HC_dat_mr <- read.csv(
  "Cplus3HC_dat_mr.csv", 
  header = TRUE)
Cplus3HC_dat_mr_5e6 <- read.csv(
  "Cplus3HC_dat_mr_5e6.csv", 
  header = TRUE)
cpd_dat_mr <- read.csv(
  "cpd_dat_mr.csv",
  header = TRUE)
Cplus3HC_dat <- read.csv(
  "Cplus3HC_dat.csv", 
  header = TRUE)
Cplus3HC_instr <- Cplus3HC_dat_mr[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "samplesize.exposure"
)]
Cplus3HC_5e6_instr <- Cplus3HC_dat_mr_5e6[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "samplesize.exposure"
)]
CPD_instr <- cpd_dat_mr[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "samplesize.exposure"
)]

INSTR <- do.call("rbind", list(Cplus3HC_instr, CPD_instr))
INSTR_5e6 <- do.call("rbind", list(Cplus3HC_5e6_instr, CPD_instr))

################################################################################
##### Check for overlapping SNPs between the exposure instruments #####
################################################################################

n_occur <- data.frame(table(INSTR$SNP))
n_occur[n_occur$Freq > 1, ]
INSTR[INSTR$SNP %in% n_occur$Var1[n_occur$Freq > 1], ]

n_occur_5e6 <- data.frame(table(INSTR_5e6$SNP))
n_occur_5e6[n_occur_5e6$Freq > 1, ]
INSTR_5e6[INSTR_5e6$SNP %in% n_occur_5e6$Var1[n_occur_5e6$Freq > 1], ]

################################################################################
##### Extract instruments from the Cplus3HC data set #####
################################################################################

Cplus3HC_dat_mvmr <- format_data(
  Cplus3HC_dat,
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
Cplus3HC_dat_mvmr$id.exposure <- "1"
str(Cplus3HC_dat_mvmr)


################################################################################
##### Extract 5e6 instruments from the Cplus3HC data set #####
################################################################################

Cplus3HC_dat_mvmr_5e6 <- format_data(
  Cplus3HC_dat,
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
Cplus3HC_dat_mvmr_5e6$id.exposure <- "1"
str(Cplus3HC_dat_mvmr_5e6)

################################################################################
##### Extract instruments from the CPD data set #####
################################################################################

cpd_dat_mvmr <- format_data(
  cpd_dat,
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
cpd_dat_mvmr$id.exposure <- "2"
str(cpd_dat_mvmr)

################################################################################
##### Extract 5e6 instruments from the CPD data set #####
################################################################################

cpd_dat_mvmr_5e6 <- format_data(
  cpd_dat,
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
cpd_dat_mvmr_5e6$id.exposure <- "2"
str(cpd_dat_mvmr_5e6)

################################################################################
##### Merge them to extract from outcome #####
################################################################################

Cplus3HC_dat_mvmr_1 <- subset(Cplus3HC_dat_mvmr, select = c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure",
  "samplesize.exposure"
))
cpd_dat_mvmr_1 <- subset(cpd_dat_mvmr, select = c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure",
  "samplesize.exposure"
))

Cplus3HC_dat_mvmr_5e6_1 <- subset(Cplus3HC_dat_mvmr_5e6, select = c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure",
  "samplesize.exposure"
))
cpd_dat_mvmr_5e6_1 <- subset(cpd_dat_mvmr_5e6, select = c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure",
  "samplesize.exposure"
))


##### Check structure is the same #####
str(Cplus3HC_dat_mvmr_1)
str(cpd_dat_mvmr_1)
str(Cplus3HC_dat_mvmr_5e6_1)
str(cpd_dat_mvmr_5e6_1)

################################################################################
##### Merge ####
################################################################################

exposures <- do.call("rbind", list(Cplus3HC_dat_mvmr_1, cpd_dat_mvmr_1))
exposures_5e6 <- do.call("rbind", list(
  Cplus3HC_dat_mvmr_5e6_1, cpd_dat_mvmr_5e6_1
))

##### Save dataframes #####
write.csv(
  exposures,
  "Cplus3HC_CPD.csv",
  row.names = FALSE
)
write.csv(
  exposures_5e6,
  "Cplus3HC_CPD_5e6.csv",
  row.names = FALSE
)

################################################################################
##### Find proxies missing from either the Cplus3HC or CPD dataset to add to 
# missing outcome SNPs #####
# Note: rs4886550 is not available in the ref panel so do not need to search for 
# proxies needed.
################################################################################

n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq < 2, ]
exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq < 2], ]

n_occur_5e6 <- data.frame(table(exposures_5e6$SNP))
n_occur_5e6[n_occur_5e6$Freq < 2, ]
exposures_5e6[exposures_5e6$SNP %in% n_occur_5e6$Var1[n_occur_5e6$Freq < 2], ]

##### Note: No proxies needed but rs4886550 is not in reference panel. #####
##### Deleted and skiped proxy search #####

# source("proxy_search_loop_CancerMVMR2a_Cplus3HC.R")

exposures <- exposures[!(exposures$SNP == "rs4886550"), ]
exposures_5e6 <- exposures_5e6[!(exposures_5e6$SNP == "rs4886550"), ]

################################################################################
##### Re-read exposure data with proxies #####
# Note: Skipped because no proxies needed.
################################################################################
# cpd_proxies <- format_data(
#  cpd_dat,
#  type = "exposure",
#  snps = SNPs.proxies8b$RS_Number,
#  header = TRUE,
#  phenotype_col = "Phenotype",
#  snp_col = "SNP",
#  beta_col = "beta.exposure",
#  se_col = "se.exposure",
#  eaf_col = "eaf.exposure",
#  effect_allele_col = "effect_allele.exposure",
#  other_allele_col = "other_allele.exposure",
#  pval_col = "pval.exposure",
#  samplesize_col = "samplesize.exposure",
#  min_pval = 1e-200,
#  log_pval = FALSE
# )

# Cplus3HC_proxies <- format_data(
#  Cplus3HC_dat,
#  type = "exposure",
#  snps = SNPs.proxies8a$RS_Number,
#  header = TRUE,
#  phenotype_col = "Phenotype",
#  snp_col = "SNP",
#  beta_col = "beta.exposure",
#  se_col = "se.exposure",
#  eaf_col = "eaf.exposure",
#  effect_allele_col = "effect_allele.exposure",
#  other_allele_col = "other_allele.exposure",
#  pval_col = "pval.exposure",
#  samplesize_col = "samplesize.exposure",
#  min_pval = 1e-200,
#  log_pval = FALSE
# )

# Cplus3HC_dat_mvmr_1 <- merge(Cplus3HC_dat_mvmr_1, Cplus3HC_proxies, all=TRUE)
# Cplus3HC_dat_mvmr_5e6_1 <- merge(Cplus3HC_dat_mvmr_5e6_1, Cplus3HC_proxies, all=TRUE)
# cpd_dat_mvmr_5e6_1 <- merge(cpd_dat_mvmr_5e6_1, cpd_proxies, all=TRUE)

# Cplus3HC_dat_mvmr_1$id.exposure<-"1"
# Cplus3HC_dat_mvmr_5e6_1$id.exposure<-"2"
# cpd_dat_mvmr_5e6_1$id.exposure<-"2"


################################################################################
##### Merge the datasets with proxies in #####
# Add pval_origin.exposure to dataframe to make columns equal
################################################################################

# cpd_dat_mvmr_1$pval_origin.exposure<-NA
# exposures<- do.call("rbind", list(Cplus3HC_dat_mvmr_1, cpd_dat_mvmr_1))
# exposures_5e6<- do.call("rbind", list(Cplus3HC_dat_mvmr_5e6_1, cpd_dat_mvmr_5e6_1))

################################################################################
##### Clumping ####
################################################################################

##### Change all p-values for Cplus3HC to 1e-200 for clumping so that none are dropped #####
##### Save old p-values first #####
exposures$oldpvalues <- exposures$pval.exposure
exposures_5e6$oldpvalues <- exposures_5e6$pval.exposure
exposures <- exposures %>%
  mutate(pval.exposure = if_else(
    exposures$SNP %in% Cplus3HC_instr$SNP,
    1e-201,
    pval.exposure
  ))

exposures_5e6 <- exposures_5e6 %>%
  mutate(pval.exposure = if_else(
    exposures_5e6$SNP %in% Cplus3HC_5e6_instr$SNP,
    1e-201, pval.exposure
  ))

##### Clump the data #####
exposures$id.exposure[exposures$id.exposure == "2"] <- "1"
exposures <- clump_data(exposures, clump_kb = 500, clump_r2 = 0.1)
str(exposures)

exposures_5e6$id.exposure[exposures_5e6$id.exposure == "2"] <- "1"
exposures_5e6 <- clump_data(exposures_5e6, clump_kb = 500, clump_r2 = 0.1)
str(exposures_5e6)

##### Add ID's back #####
exposures$id.exposure[exposures$samplesize.exposure < 6000] <- "1"
exposures$id.exposure[exposures$samplesize.exposure > 6000] <- "2"

exposures_5e6$id.exposure[exposures_5e6$samplesize.exposure < 6000] <- "1"
exposures_5e6$id.exposure[exposures_5e6$samplesize.exposure > 6000] <- "2"

##### Revert all p-values for Cplus3HC from 1e-200 #####
exposures$pval.exposure <- exposures$oldpvalues
exposures <- select(exposures, -c(oldpvalues))

exposures_5e6$pval.exposure <- exposures_5e6$oldpvalues
exposures_5e6 <- select(exposures_5e6, -c(oldpvalues))

##### Split again to harmonise based on exposure id #####
Cplus3HC <- split(exposures, exposures$id.exposure)[["1"]]
CPD <- split(exposures, exposures$id.exposure)[["2"]]

Cplus3HC_5e6 <- split(exposures_5e6, exposures_5e6$id.exposure)[["1"]]
CPD_5e6 <- split(exposures_5e6, exposures_5e6$id.exposure)[["2"]]

################################################################################
##### Harmonise Cplus3HC on CPD #####
################################################################################

names(Cplus3HC) <- gsub("exposure", "outcome", names(Cplus3HC))
Cplus3HC_CPD <- harmonise_data(CPD, Cplus3HC)

names(Cplus3HC_5e6) <- gsub("exposure", "outcome", names(Cplus3HC_5e6))
Cplus3HC_CPD_5e6 <- harmonise_data(CPD_5e6, Cplus3HC_5e6)

################################################################################
##### Keep only snps that are present across both exposures ####
# Note: they would have frequency 1 if only available in one dataset
################################################################################

n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq == 2, ]
exposures <- exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq == 2], ]
str(exposures)

n_occur_5e6 <- data.frame(table(exposures_5e6$SNP))
n_occur_5e6[n_occur_5e6$Freq == 2, ]
exposures_5e6 <- exposures_5e6[exposures_5e6$SNP %in% n_occur_5e6$Var1[n_occur_5e6$Freq == 2], ]
str(exposures_5e6)

################################################################################
##### Format exposures #####
################################################################################

##### Keep only snps MrKeep= TRUE #####
Cplus3HC_CPD <- Cplus3HC_CPD[Cplus3HC_CPD$mr_keep == TRUE, ]
str(Cplus3HC_CPD)

Cplus3HC_CPD_5e6 <- Cplus3HC_CPD_5e6[Cplus3HC_CPD_5e6$mr_keep == TRUE, ]
str(Cplus3HC_CPD_5e6)

##### Split the tables - CPD #####
CPD_H <- subset(Cplus3HC_CPD, id.exposure == "2", select = c(
  SNP,
  exposure,
  id.exposure,
  effect_allele.exposure,
  other_allele.exposure,
  beta.exposure,
  se.exposure,
  pval.exposure,
  eaf.exposure
))
CPD_H_5e6 <- subset(Cplus3HC_CPD_5e6, id.exposure == "2", select = c(
  SNP,
  exposure,
  id.exposure,
  effect_allele.exposure,
  other_allele.exposure,
  beta.exposure,
  se.exposure,
  pval.exposure,
  eaf.exposure
))

##### Split the tables - Cplus3HC #####
Cplus3HC_H <- subset(Cplus3HC_CPD, id.outcome == "1", select = c(
  SNP,
  outcome,
  id.outcome,
  effect_allele.outcome,
  other_allele.outcome,
  beta.outcome,
  se.outcome,
  pval.outcome,
  eaf.outcome
))
Cplus3HC_H_5e6 <- subset(Cplus3HC_CPD_5e6, id.outcome == "1", select = c(
  SNP,
  outcome,
  id.outcome,
  effect_allele.outcome,
  other_allele.outcome,
  beta.outcome,
  se.outcome,
  pval.outcome,
  eaf.outcome
))

##### Turn Cplus3HC from outcome to exposure to merge the datasets #####
names(Cplus3HC_H) <- gsub("outcome", "exposure", names(Cplus3HC_H))
Exposures_H <- merge(Cplus3HC_H, CPD_H, all = TRUE)
Exposures_H["Phenotype"] <- NA
Exposures_H$Phenotype[Exposures_H$id.exposure == 1] <- "Cplus3HC"
Exposures_H$Phenotype[Exposures_H$id.exposure == 2] <- "CPD"
str(Exposures_H)

names(Cplus3HC_H_5e6) <- gsub("outcome", "exposure", names(Cplus3HC_H_5e6))
Exposures_H_5e6 <- merge(Cplus3HC_H_5e6, CPD_H_5e6, all = TRUE)
Exposures_H_5e6["Phenotype"] <- NA
Exposures_H_5e6$Phenotype[Exposures_H_5e6$id.exposure == 1] <- "Cplus3HC"
Exposures_H_5e6$Phenotype[Exposures_H_5e6$id.exposure == 2] <- "CPD"
str(Exposures_H_5e6)

################################################################################
##### Extract outcome data for MVMR #####
################################################################################

Exposures_H[Exposures_H == "rs77107237"] <- "rs114612145"
Exposures_H_5e6[Exposures_H_5e6 == "rs77107237"] <- "rs114612145"

setwd("")
memory.limit(size = 80000)
outcome_dat_never <- read_outcome_data(
  "Never.txt",
  snps = Exposures_H$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd("")
memory.limit(size = 80000)
outcome_dat_never_5e6 <- read_outcome_data(
  "Never.txt",
  snps = Exposures_H_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd("")
memory.limit(size = 80000)
outcome_dat_ever <- read_outcome_data(
  "Ever.txt",
  snps = Exposures_H$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 40187,
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd("")
memory.limit(size = 80000)
outcome_dat_ever_5e6 <- read_outcome_data(
  "Ever.txt",
  snps = Exposures_H_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 40187,
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Convert odds ratios to log odds #####
################################################################################

outcome_dat_never$beta.outcome <- as.numeric(
  as.character(outcome_dat_never$beta.outcome)
)
outcome_dat_ever$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever$beta.outcome)
)

outcome_dat_never["beta.outcome"] <- log(
  outcome_dat_never["beta.outcome"]
)
outcome_dat_ever["beta.outcome"] <- log(
  outcome_dat_ever["beta.outcome"]
)

outcome_dat_never_5e6$beta.outcome <- as.numeric(
  as.character(outcome_dat_never_5e6$beta.outcome)
)
outcome_dat_ever_5e6$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever_5e6$beta.outcome)
)

outcome_dat_never_5e6["beta.outcome"] <- log(
  outcome_dat_never_5e6["beta.outcome"]
)
outcome_dat_ever_5e6["beta.outcome"] <- log(
  outcome_dat_ever_5e6["beta.outcome"]
)

################################################################################
##### Organise outcome #####
################################################################################

outcome_dat_never["Phenotype"] <- NA
outcome_dat_never$Phenotype <- "LungCancer"

outcome_dat_ever["Phenotype"] <- NA
outcome_dat_ever$Phenotype <- "LungCancer"

outcome_dat_never_5e6["Phenotype"] <- NA
outcome_dat_never_5e6$Phenotype <- "LungCancer"

outcome_dat_ever_5e6["Phenotype"] <- NA
outcome_dat_ever$Phenotype_5e6 <- "LungCancer"

################################################################################
##### Harmonise with outcome #####
################################################################################

mvdat_never <- mv_harmonise_data(
  Exposures_H,
  outcome_dat_never,
  harmonise_strictness = 2
)
mvdat_ever <- mv_harmonise_data(
  Exposures_H,
  outcome_dat_ever,
  harmonise_strictness = 2
)

mvdat_never1 <- harmonise_data(
  Exposures_H,
  outcome_dat_never
)
mvdat_never1 <- mvdat_never1[mvdat_never1$mr_keep == TRUE, ]
str(mvdat_never1)

mvdat_ever1 <- harmonise_data(Exposures_H, outcome_dat_ever)
mvdat_ever1 <- mvdat_ever1[mvdat_ever1$mr_keep == TRUE, ]
str(mvdat_ever1)

mvdat_never_5e6 <- mv_harmonise_data(
  Exposures_H_5e6,
  outcome_dat_never_5e6,
  harmonise_strictness = 2
)
mvdat_ever_5e6 <- mv_harmonise_data(
  Exposures_H_5e6,
  outcome_dat_ever_5e6,
  harmonise_strictness = 2
)

mvdat_never1_5e6 <- harmonise_data(
  Exposures_H_5e6,
  outcome_dat_never_5e6
)
mvdat_never1_5e6 <- mvdat_never1_5e6[mvdat_never1_5e6$mr_keep == TRUE, ]
str(mvdat_never1_5e6)

mvdat_ever1_5e6 <- harmonise_data(Exposures_H_5e6, outcome_dat_ever_5e6)
mvdat_ever1_5e6 <- mvdat_ever1_5e6[mvdat_ever1_5e6$mr_keep == TRUE, ]
str(mvdat_ever1_5e6)

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed9 <- data.frame(setdiff(Cplus3HC_H$SNP, outcome_dat_never$SNP))
proxy_needed10 <- data.frame(setdiff(Cplus3HC_H$SNP, outcome_dat_ever$SNP))
proxy_needed11 <- data.frame(setdiff(
  Cplus3HC_H_5e6$SNP,
  outcome_dat_never_5e6$SNP
))
proxy_needed12 <- data.frame(setdiff(
  Cplus3HC_H_5e6$SNP,
  outcome_dat_ever_5e6$SNP
))

source("proxy_search_loop_CancerMVMR2_Cplus3HC.R")

################################################################################
##### Re-merge the data #####
################################################################################

Exposures_H_never <- merge(Cplus3HC_never, CPD_never, all = TRUE)
Exposures_H_never["Phenotype"] <- NA
Exposures_H_never$Phenotype[Exposures_H_never$id.exposure == 1] <- "Cplus3HC"
Exposures_H_never$Phenotype[Exposures_H_never$id.exposure == 2] <- "CPD"
str(Exposures_H_never)

Exposures_H_never_5e6 <- merge(Cplus3HC_never_5e6, CPD_never_5e6, all = TRUE)
Exposures_H_never_5e6["Phenotype"] <- NA
Exposures_H_never_5e6$Phenotype[Exposures_H_never_5e6$id.exposure == 1] <- "Cplus3HC"
Exposures_H_never_5e6$Phenotype[Exposures_H_never_5e6$id.exposure == 2] <- "CPD"
str(Exposures_H_never_5e6)

Exposures_H_ever <- merge(Cplus3HC_ever, CPD_ever, all = TRUE)
Exposures_H_ever["Phenotype"] <- NA
Exposures_H_ever$Phenotype[Exposures_H_ever$id.exposure == 1] <- "Cplus3HC"
Exposures_H_ever$Phenotype[Exposures_H_ever$id.exposure == 2] <- "CPD"
str(Exposures_H_ever)

Exposures_H_ever_5e6 <- merge(Cplus3HC_ever_5e6, CPD_ever_5e6, all = TRUE)
Exposures_H_ever_5e6["Phenotype"] <- NA
Exposures_H_ever_5e6$Phenotype[Exposures_H_ever_5e6$id.exposure == 1] <- "Cplus3HC"
Exposures_H_ever_5e6$Phenotype[Exposures_H_ever_5e6$id.exposure == 2] <- "CPD"
str(Exposures_H_ever_5e6)

################################################################################
##### Re-read the outcome data #####
################################################################################

Exposures_H_never[Exposures_H_never == "rs77107237"] <- "rs114612145"
Exposures_H_never_5e6[Exposures_H_never_5e6 == "rs77107237"] <- "rs114612145"
Exposures_H_ever[Exposures_H_ever == "rs77107237"] <- "rs114612145"
Exposures_H_ever_5e6[Exposures_H_ever_5e6 == "rs77107237"] <- "rs114612145"

setwd("")
memory.limit(size = 80000)
outcome_dat_never <- read_outcome_data(
  "Never.txt",
  snps = Exposures_H_never$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd("")
memory.limit(size = 80000)
outcome_dat_never_5e6 <- read_outcome_data(
  "Never.txt",
  snps = Exposures_H_never_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd("")
memory.limit(size = 80000)
outcome_dat_ever <- read_outcome_data(
  "Ever.txt",
  snps = Exposures_H_ever$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 40187,
  min_pval = 1e-200,
  log_pval = FALSE
)

setwd("")
memory.limit(size = 80000)
outcome_dat_ever_5e6 <- read_outcome_data(
  "Ever.txt",
  snps = Exposures_H_ever_5e6$SNP,
  sep = ",",
  snp_col = "rs_number",
  beta_col = "OR_random",
  se_col = "StdError_random",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "reference_allele",
  pval_col = "Pvalue_random",
  samplesize_col = 40187,
  min_pval = 1e-200,
  log_pval = FALSE
)

################################################################################
##### Convert odds ratios to log odds #####
################################################################################

outcome_dat_never$beta.outcome <- as.numeric(
  as.character(outcome_dat_never$beta.outcome)
)
outcome_dat_ever$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever$beta.outcome)
)

outcome_dat_never["beta.outcome"] <- log(outcome_dat_never["beta.outcome"])
outcome_dat_ever["beta.outcome"] <- log(outcome_dat_ever["beta.outcome"])

outcome_dat_never_5e6$beta.outcome <- as.numeric(
  as.character(outcome_dat_never_5e6$beta.outcome)
)
outcome_dat_ever_5e6$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever_5e6$beta.outcome)
)

outcome_dat_never_5e6["beta.outcome"] <- log(
  outcome_dat_never_5e6["beta.outcome"]
)
outcome_dat_ever_5e6["beta.outcome"] <- log(
  outcome_dat_ever_5e6["beta.outcome"]
)

################################################################################
##### Organise outcome #####
################################################################################

outcome_dat_never["Phenotype"] <- NA
outcome_dat_never$Phenotype <- "LungCancer"

outcome_dat_ever["Phenotype"] <- NA
outcome_dat_ever$Phenotype <- "LungCancer"

outcome_dat_never_5e6["Phenotype"] <- NA
outcome_dat_never_5e6$Phenotype <- "LungCancer"

outcome_dat_ever_5e6["Phenotype"] <- NA
outcome_dat_ever$Phenotype_5e6 <- "LungCancer"

################################################################################
##### Harmonise with outcome #####
################################################################################

mvdat_never <- mv_harmonise_data(
  Exposures_H_never,
  outcome_dat_never,
  harmonise_strictness = 2
)
mvdat_ever <- mv_harmonise_data(
  Exposures_H_ever,
  outcome_dat_ever,
  harmonise_strictness = 2
)

mvdat_never1 <- harmonise_data(Exposures_H_never, outcome_dat_never)
mvdat_never1 <- mvdat_never1[mvdat_never1$mr_keep == TRUE, ]
str(mvdat_never1)

mvdat_ever1 <- harmonise_data(Exposures_H_ever, outcome_dat_ever)
mvdat_ever1 <- mvdat_ever1[mvdat_ever1$mr_keep == TRUE, ]
str(mvdat_ever1)

mvdat_never_5e6 <- mv_harmonise_data(
  Exposures_H_never_5e6,
  outcome_dat_never_5e6,
  harmonise_strictness = 2
)
mvdat_ever_5e6 <- mv_harmonise_data(
  Exposures_H_ever_5e6,
  outcome_dat_ever_5e6,
  harmonise_strictness = 2
)

mvdat_never1_5e6 <- harmonise_data(Exposures_H_never_5e6, outcome_dat_never_5e6)
mvdat_never1_5e6 <- mvdat_never1_5e6[mvdat_never1_5e6$mr_keep == TRUE, ]
str(mvdat_never1_5e6)

mvdat_ever1_5e6 <- harmonise_data(Exposures_H_ever_5e6, outcome_dat_ever_5e6)
mvdat_ever1_5e6 <- mvdat_ever1_5e6[mvdat_ever1_5e6$mr_keep == TRUE, ]
str(mvdat_ever1_5e6)

################################################################################
##### Save dataframes #####
################################################################################

write.csv(
  mvdat_never1,
  "mvdat_never_Cplus3HC.csv",
  row.names = FALSE
)
write.csv(
  mvdat_ever1,
  "mvdat_ever_Cplus3HC.csv",
  row.names = FALSE
)

write.csv(
  mvdat_never1_5e6,
  "mvdat_never_Cplus3HC5e6.csv",
  row.names = FALSE
)
write.csv(
  mvdat_ever1_5e6,
  "mvdat_ever_Cplus3HC5e6.csv",
  row.names = FALSE
)

################################################################################
##### Run MVMR #####
################################################################################

##### IVW never 5e8 #####
bX1 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 1])
bX2 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 2])
bY <- c(mvdat_never1$beta.outcome[mvdat_never1$id.exposure == 1])
bYse <- c(mvdat_never1$se.outcome[mvdat_never1$id.exposure == 1])

set.seed(1234)
mod.MVMR_never <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$sigma, 1)

mod_n <- summary(mod.MVMR_never)

mod_n_or <- coef(summary(mod.MVMR_never))
colnames(mod_n_or) <- c("b", "se", "t", "p")
mod_n_or <- as.data.frame(mod_n_or)
mod_n_or <- generate_odds_ratios(mod_n_or)

##### Orientation Cplus3HC #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for Cplus3HC, and then     #####
##### orient the betas to be positive for CPD. In the paper, we           #####
##### report the result for each exposure only with the right orientation #####
clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_never <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_n <- summary(mod.MVMRME_never)

mod_ME_n_or <- data.frame(mod.MVMRME_never[["coefficients"]])
colnames(mod_ME_n_or) <- c("b", "se", "t", "p")
mod_ME_n_or <- as.data.frame(mod_ME_n_or)
mod_ME_n_or <- generate_odds_ratios(mod_ME_n_or)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_never_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_n_2 <- summary(mod.MVMRME_never_2)

mod_ME_n_2_or <- data.frame(mod.MVMRME_never_2[["coefficients"]])
colnames(mod_ME_n_2_or) <- c("b", "se", "t", "p")
mod_ME_n_2_or <- as.data.frame(mod_ME_n_2_or)
mod_ME_n_2_or <- generate_odds_ratios(mod_ME_n_2_or)

##### IVW never 5e6 #####
bX1 <- c(mvdat_never1_5e6$beta.exposure[mvdat_never1_5e6$id.exposure == 1])
bX2 <- c(mvdat_never1_5e6$beta.exposure[mvdat_never1_5e6$id.exposure == 2])
bY <- c(mvdat_never1_5e6$beta.outcome[mvdat_never1_5e6$id.exposure == 1])
bYse <- c(mvdat_never1_5e6$se.outcome[mvdat_never1_5e6$id.exposure == 1])
set.seed(1234)
mod.MVMR_never_5e6 <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$sigma, 1)

mod_n_5e6 <- summary(mod.MVMR_never_5e6)

mod_n_or_5e6 <- coef(summary(mod.MVMR_never_5e6))
colnames(mod_n_or_5e6) <- c("b", "se", "t", "p")
mod_n_or_5e6 <- as.data.frame(mod_n_or_5e6)
mod_n_or_5e6 <- generate_odds_ratios(mod_n_or_5e6)

##### Orientation Cplus3HC #####
clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_never_5e6 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_n_5e6 <- summary(mod.MVMRME_never_5e6)

mod_ME_n_or_5e6 <- data.frame(mod.MVMRME_never_5e6[["coefficients"]])
colnames(mod_ME_n_or_5e6) <- c("b", "se", "t", "p")
mod_ME_n_or_5e6 <- as.data.frame(mod_ME_n_or)
mod_ME_n_or_5e6 <- generate_odds_ratios(mod_ME_n_or_5e6)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_never_2_5e6 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_n_2_5e6 <- summary(mod.MVMRME_never_2_5e6)

mod_ME_n_2_or_5e6 <- data.frame(mod.MVMRME_never_2_5e6[["coefficients"]])
colnames(mod_ME_n_2_or_5e6) <- c("b", "se", "t", "p")
mod_ME_n_2_or_5e6 <- as.data.frame(mod_ME_n_2_or_5e6)
mod_ME_n_2_or_5e6 <- generate_odds_ratios(mod_ME_n_2_or_5e6)

##### IVW ever 5e8 #####
bX1 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 1])
bX2 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 2])
bY <- c(mvdat_ever1$beta.outcome[mvdat_ever1$id.exposure == 1])
bYse <- c(mvdat_ever1$se.outcome[mvdat_ever1$id.exposure == 1])

set.seed(1234)
mod.MVMR_ever <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$sigma, 1)

mod_e <- summary(mod.MVMR_ever)

mod_e_or <- coef(summary(mod.MVMR_ever))
colnames(mod_e_or) <- c("b", "se", "t", "p")
mod_e_or <- as.data.frame(mod_e_or)
mod_e_or <- generate_odds_ratios(mod_e_or)

##### Orientation Cplus3HC #####
clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_ever <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_e <- summary(mod.MVMRME_ever)

mod_ME_e_or <- data.frame(mod.MVMRME_ever[["coefficients"]])
colnames(mod_ME_e_or) <- c("b", "se", "t", "p")
mod_ME_e_or <- as.data.frame(mod_ME_e_or)
mod_ME_e_or <- generate_odds_ratios(mod_ME_e_or)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_ever_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_e_2 <- summary(mod.MVMRME_ever_2)

mod_ME_e_2_or <- data.frame(mod.MVMRME_ever_2[["coefficients"]])
colnames(mod_ME_e_2_or) <- c("b", "se", "t", "p")
mod_ME_e_2_or <- as.data.frame(mod_ME_e_2_or)
mod_ME_e_2_or <- generate_odds_ratios(mod_ME_e_2_or)

##### IVW ever 5e6 #####
bX1 <- c(mvdat_ever1_5e6$beta.exposure[mvdat_ever1_5e6$id.exposure == 1])
bX2 <- c(mvdat_ever1_5e6$beta.exposure[mvdat_ever1_5e6$id.exposure == 2])
bY <- c(mvdat_ever1_5e6$beta.outcome[mvdat_ever1_5e6$id.exposure == 1])
bYse <- c(mvdat_ever1_5e6$se.outcome[mvdat_ever1_5e6$id.exposure == 1])

set.seed(1234)
mod.MVMR_ever_5e6 <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$sigma, 1)

mod_e_5e6 <- summary(mod.MVMR_ever_5e6)

mod_e_or_5e6 <- coef(summary(mod.MVMR_ever_5e6))
colnames(mod_e_or_5e6) <- c("b", "se", "t", "p")
mod_e_or_5e6 <- as.data.frame(mod_e_or_5e6)
mod_e_or_5e6 <- generate_odds_ratios(mod_e_or_5e6)

##### Orientation Cplus3HC #####
clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_ever_5e6 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_e_5e6 <- summary(mod.MVMRME_ever_5e6)

mod_ME_e_or_5e6 <- data.frame(mod.MVMRME_ever_5e6[["coefficients"]])
colnames(mod_ME_e_or_5e6) <- c("b", "se", "t", "p")
mod_ME_e_or_5e6 <- as.data.frame(mod_ME_e_or_5e6)
mod_ME_e_or_5e6 <- generate_odds_ratios(mod_ME_e_or_5e6)

##### Orientation CPD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_ever_2_5e6 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_e_2_5e6 <- summary(mod.MVMRME_ever_2_5e6)


mod_ME_e_2_or_5e6 <- data.frame(mod.MVMRME_ever_2_5e6[["coefficients"]])
colnames(mod_ME_e_2_or_5e6) <- c("b", "se", "t", "p")
mod_ME_e_2_or_5e6 <- as.data.frame(mod_ME_e_2_or_5e6)
mod_ME_e_2_or_5e6 <- generate_odds_ratios(mod_ME_e_2_or_5e6)

##### Format to analyse with MVMR package #####

bX1 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 1])
bX2 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 2])
bY <- c(mvdat_never1$beta.outcome[mvdat_never1$id.exposure == 1])
bYse <- c(mvdat_never1$se.outcome[mvdat_never1$id.exposure == 1])
bXse1 <- c(mvdat_never1$se.exposure[mvdat_never1$id.exposure == 1])
bXse2 <- c(mvdat_never1$se.exposure[mvdat_never1$id.exposure == 2])
df_never <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

bX1 <- c(mvdat_never1_5e6$beta.exposure[mvdat_never1_5e6$id.exposure == 1])
bX2 <- c(mvdat_never1_5e6$beta.exposure[mvdat_never1_5e6$id.exposure == 2])
bY <- c(mvdat_never1_5e6$beta.outcome[mvdat_never1_5e6$id.exposure == 1])
bYse <- c(mvdat_never1_5e6$se.outcome[mvdat_never1_5e6$id.exposure == 1])
bXse1 <- c(mvdat_never1_5e6$se.exposure[mvdat_never1_5e6$id.exposure == 1])
bXse2 <- c(mvdat_never1_5e6$se.exposure[mvdat_never1_5e6$id.exposure == 2])
df_never_5e6 <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

bX1 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 1])
bX2 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 2])
bY <- c(mvdat_ever1$beta.outcome[mvdat_ever1$id.exposure == 1])
bYse <- c(mvdat_ever1$se.outcome[mvdat_ever1$id.exposure == 1])
bXse1 <- c(mvdat_ever1$se.exposure[mvdat_ever1$id.exposure == 1])
bXse2 <- c(mvdat_ever1$se.exposure[mvdat_ever1$id.exposure == 2])
df_ever <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

bX1 <- c(mvdat_ever1_5e6$beta.exposure[mvdat_ever1_5e6$id.exposure == 1])
bX2 <- c(mvdat_ever1_5e6$beta.exposure[mvdat_ever1_5e6$id.exposure == 2])
bY <- c(mvdat_ever1_5e6$beta.outcome[mvdat_ever1_5e6$id.exposure == 1])
bYse <- c(mvdat_ever1_5e6$se.outcome[mvdat_ever1_5e6$id.exposure == 1])
bXse1 <- c(mvdat_ever1_5e6$se.exposure[mvdat_ever1_5e6$id.exposure == 1])
bXse2 <- c(mvdat_ever1_5e6$se.exposure[mvdat_ever1_5e6$id.exposure == 2])
df_ever_5e6 <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)
df_mvmr_never <- format_mvmr(
  df_never[, c(1, 3)],
  df_never[, 5],
  df_never[, c(2, 4)],
  df_never[, 6]
)
df_mvmr_ever <- format_mvmr(
  df_ever[, c(1, 3)],
  df_ever[, 5],
  df_ever[, c(2, 4)],
  df_ever[, 6]
)
df_mvmr_never_5e6 <- format_mvmr(
  df_never_5e6[, c(1, 3)],
  df_never_5e6[, 5],
  df_never_5e6[, c(2, 4)],
  df_never_5e6[, 6]
)
df_mvmr_ever_5e6 <- format_mvmr(
  df_ever_5e6[, c(1, 3)],
  df_ever_5e6[, 5],
  df_ever_5e6[, c(2, 4)],
  df_ever_5e6[, 6]
)

##### Cross check my result with MVMR package #####
res_n <- ivw_mvmr(df_mvmr_never)
res_e <- ivw_mvmr(df_mvmr_ever)
res_n_5e6 <- ivw_mvmr(df_mvmr_never_5e6)
res_e_5e6 <- ivw_mvmr(df_mvmr_ever_5e6)

################################################################################
##### Calculate F-statistic and covariance #####
# Note: >10 is strong
# correlation between Cplus3HC and CPD in Buchwald  = 0.468
# qhet_mvmr is used to adjust for covariance but at present, CIs take substantial time to calculate
# Compare effects with and without adjustment
################################################################################

cov <- matrix(c(1, 0.468, 0.468, 1), nrow = 2, ncol = 2)
Xcovmat_n <- phenocov_mvmr(cov, df_mvmr_never[, c(6, 7)])
Fstat_n <- strength_mvmr(df_mvmr_never, gencov = Xcovmat_n)
Xcovmat_e <- phenocov_mvmr(cov, df_mvmr_ever[, c(6, 7)])
Fstat_e <- strength_mvmr(df_mvmr_ever, gencov = Xcovmat_e)
Xcovmat_n_5e6 <- phenocov_mvmr(cov, df_mvmr_never_5e6[, c(6, 7)])
Fstat_n_5e6 <- strength_mvmr(df_mvmr_never_5e6, gencov = Xcovmat_n_5e6)
Xcovmat_e_5e6 <- phenocov_mvmr(cov, df_mvmr_ever_5e6[, c(6, 7)])
Fstat_e_5e6 <- strength_mvmr(df_mvmr_ever_5e6, gencov = Xcovmat_e_5e6)

cov_adj_mvmr_n <- qhet_mvmr(
  df_mvmr_never,
  cov,
  CI = FALSE,
  iterations = 1000
)
cov_adj_mvmr_e <- qhet_mvmr(
  df_mvmr_ever,
  cov,
  CI = FALSE,
  iterations = 1000
)
cov_adj_mvmr_n_5e6 <- qhet_mvmr(
  df_mvmr_never_5e6,
  cov,
  CI = FALSE,
  iterations = 1000
)
cov_adj_mvmr_e_5e6 <- qhet_mvmr(
  df_mvmr_ever_5e6,
  cov,
  CI = FALSE,
  iterations = 1000
)

################################################################################
##### Test for horizontal pleiotropy #####
# Note: Q should be greater than the number of SNPs included
################################################################################

ptr_n <- pleiotropy_mvmr(df_mvmr_never, gencov = Xcovmat_n)
ptr_e <- pleiotropy_mvmr(df_mvmr_ever, gencov = Xcovmat_e)
ptr_n_5e6 <- pleiotropy_mvmr(df_mvmr_never_5e6, gencov = Xcovmat_n_5e6)
ptr_e_5e6 <- pleiotropy_mvmr(df_mvmr_ever_5e6, gencov = Xcovmat_e_5e6)

################################################################################
##### Forest Plots #####
# Only plotting 5e6 because egger for 5e8 because the I2GX was below 0.6
# The results could be severely impacted by NOME violation and adjusting may be worse
################################################################################

Exposure <- c(
  "Cplus3HC", "Cplus3HC", "Cplus3HC", "Cplus3HC",
  "CPD", "CPD", "CPD", "CPD",
  "Cplus3HC", "Cplus3HC", "Cplus3HC", "Cplus3HC",
  "CPD", "CPD", "CPD", "CPD"
)
Smoking <- c(
  "Never", "Never", "Never", "Never", "Never", "Never", "Never", "Never",
  "Ever", "Ever", "Ever", "Ever", "Ever", "Ever", "Ever", "Ever"
)
res_n <- data.frame(res_n)
res_e <- data.frame(res_e)
simexegger_Cplus3HCn5e6 <- data.frame(simexegger_Cplus3HCn5e6)
simexegger_Cplus3HCe5e6 <- data.frame(simexegger_Cplus3HCe5e6)

Method <- c(
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger"
)
OR <- as.numeric(c(
  result_Cplus3HC_never_5e6[3, "or"], mod_n_or_5e6[1, 7],
  simexegger_Cplus3HCn5e6[1, "OR_weighted"], mod_ME_n_or_5e6[2, 7],
  result_cpd_never[3, "or"], mod_n_or_5e6[2, 7],
  result_cpd_never[1, "or"], mod_ME_n_2_or_5e6[3, 7],
  result_Cplus3HC_ever_5e6[3, "or"], mod_e_or_5e6[1, 7],
  simexegger_Cplus3HCe5e6[1, "OR_weighted"], mod_ME_e_or_5e6[2, 7],
  result_cpd_ever[3, "or"], mod_e_or_5e6[2, 7],
  result_cpd_ever[1, "or"], mod_ME_e_2_or_5e6[3, 7]
))

LCI <- as.numeric(c(
  result_Cplus3HC_never_5e6[3, "or_lci95"], mod_n_or_5e6[1, 8],
  simexegger_Cplus3HCn5e6[1, "lCIOR_weighted"], mod_ME_n_or_5e6[2, 8],
  result_cpd_never[3, "or_lci95"], mod_n_or_5e6[2, 8],
  result_cpd_never[1, "or_lci95"], mod_ME_n_2_or_5e6[3, 8],
  result_Cplus3HC_ever_5e6[3, "or_lci95"], mod_e_or_5e6[1, 8],
  simexegger_Cplus3HCe5e6[1, "lCIOR_weighted"], mod_ME_e_or_5e6[2, 8],
  result_cpd_ever[3, "or_lci95"], mod_e_or_5e6[2, 8],
  result_cpd_ever[1, "or_lci95"], mod_ME_e_2_or_5e6[3, 8]
))

UCI <- as.numeric(c(
  result_Cplus3HC_never_5e6[3, "or_uci95"], mod_n_or_5e6[1, 9],
  simexegger_Cplus3HCn5e6[1, "uCIOR_weighted"], mod_ME_n_or_5e6[2, 9],
  result_cpd_never[3, "or_uci95"], mod_n_or_5e6[2, 9],
  result_cpd_never[1, "or_uci95"], mod_ME_n_2_or_5e6[3, 9],
  result_Cplus3HC_ever_5e6[3, "or_uci95"], mod_e_or_5e6[1, 9],
  simexegger_Cplus3HCe5e6[1, "uCIOR_weighted"], mod_ME_e_or_5e6[2, 9],
  result_cpd_ever[3, "or_uci95"], mod_e_or_5e6[2, 9],
  result_cpd_ever[1, "or_uci95"], mod_ME_e_2_or_5e6[3, 9]
))

p <- as.numeric(c(
  result_Cplus3HC_never_5e6[3, "pval"], mod_n_or_5e6[1, 4],
  simexegger_Cplus3HCn5e6[1, "p_weighted"], mod_ME_n_or_5e6[2, 4],
  result_cpd_never[3, "pval"], mod_n_or_5e6[2, 4],
  result_cpd_never[1, "pval"], mod_ME_n_2_or_5e6[3, 4],
  result_Cplus3HC_ever_5e6[3, "pval"], mod_e_or_5e6[1, 4],
  simexegger_Cplus3HCe5e6[1, "p_weighted"], mod_ME_e_or_5e6[2, 4],
  result_cpd_ever[3, "pval"], mod_e_or_5e6[2, 4],
  result_cpd_ever[1, "pval"], mod_ME_e_2_or_5e6[3, 4]
))

I2Gx <- c(
  ".", ".", ISQ[3, 1], ".",
  ".", ".", ISQ[5, 1], ".",
  ".", ".", ISQ[4, 1], ".",
  ".", ".", ISQ[6, 1], "."
)
Q <- c(
  ptr_n_Cplus3HC_5e6[1, 1], ptr_n_5e6[["Qstat"]],
  ptr_n_Cplus3HC_5e6[2, 1], ptr_n_5e6[["Qstat"]],
  ptr_n_cpd[1, 1], ptr_n_5e6[["Qstat"]],
  ptr_n_cpd[2, 1], ptr_n_5e6[["Qstat"]],
  ptr_e_Cplus3HC_5e6[1, 1], ptr_e_5e6[["Qstat"]],
  ptr_e_Cplus3HC_5e6[2, 1], ptr_e_5e6[["Qstat"]],
  ptr_e_cpd[1, 1], ptr_e_5e6[["Qstat"]],
  ptr_e_cpd[2, 1], ptr_e_5e6[["Qstat"]]
)

EggerI <- c(
  ".", ".", egger_n_Cplus3HC_5e6[["b_i"]], mod_ME_n_or_5e6[1, 7],
  ".", ".", egger_n_cpd[["b_i"]], mod_ME_n_2_or_5e6[1, 7],
  ".", ".", egger_e_Cplus3HC_5e6[["b_i"]], mod_ME_e_or_5e6[1, 7],
  ".", ".", egger_e_cpd[["b_i"]], mod_ME_e_2_or_5e6[1, 7]
)

EggerIp <- c(
  ".", ".", egger_n_Cplus3HC_5e6[["pval_i"]], mod_ME_n_or_5e6[1, 4],
  ".", ".", egger_n_cpd[["pval_i"]], mod_ME_n_2_or_5e6[1, 4],
  ".", ".", egger_e_Cplus3HC_5e6[["pval_i"]], mod_ME_e_or_5e6[1, 4],
  ".", ".", egger_e_cpd[["pval_i"]], mod_ME_e_2_or_5e6[1, 4]
)

F_stat <- c(
  ptr_n_Cplus3HC_5e6[3, 1], Fstat_n_5e6[1, 1],
  ptr_n_Cplus3HC_5e6[3, 1], Fstat_n_5e6[1, 1],
  ptr_n_cpd[3, 1], Fstat_n_5e6[1, 2],
  ptr_n_cpd[3, 1], Fstat_n_5e6[1, 2],
  ptr_e_Cplus3HC_5e6[3, 1], Fstat_e_5e6[1, 1],
  ptr_e_Cplus3HC_5e6[3, 1], Fstat_e_5e6[1, 1],
  ptr_e_cpd[3, 1], Fstat_e_5e6[1, 2],
  ptr_e_cpd[3, 1], Fstat_e_5e6[1, 2]
)

all_results <- data.frame(
  Exposure,
  Smoking,
  Method,
  OR,
  LCI,
  UCI,
  p,
  I2Gx,
  Q,
  EggerI,
  EggerIp,
  F_stat
)

write.csv(
  all_results,
  "Results_3HCplusCOT_CPD_lung.csv",
  row.names = FALSE
)

MVMRproxies_Cplus3HC5e6_n <- subset(
  SNPs.proxies11,
  target %in% Cplus3HC_SNPlist_5e6$SNP
)
MVMRproxies_cpd_n <- subset(SNPs.proxies11, target %in% CPD_SNPlist$SNP)
MVMRproxies_Cplus3HC5e6_e <- subset(
  SNPs.proxies12,
  target %in% Cplus3HC_SNPlist_5e6$SNP
)
MVMRproxies_cpd_e <- subset(SNPs.proxies12, target %in% CPD_SNPlist$SNP)

nSNPs_Cplus3HCn_5e6 <- length(intersect(
  Cplus3HC_SNPlist_5e6$SNP, mvdat_never1_5e6$SNP
)) + length(intersect(
  MVMRproxies_Cplus3HC5e6_n$RS_Number, mvdat_never1_5e6$SNP
))
nSNPs_CPDn <- length(intersect(
  CPD_SNPlist$SNP, mvdat_never1_5e6$SNP
)) + length(intersect(
  MVMRproxies_cpd_n$RS_Number, mvdat_never1_5e6$SNP
))
nSNPs_Cplus3HCe_5e6 <- length(intersect(
  Cplus3HC_SNPlist_5e6$SNP, mvdat_ever1_5e6$SNP
)) + length(intersect(
  MVMRproxies_Cplus3HC5e6_e$RS_Number, mvdat_ever1_5e6$SNP
))
nSNPs_CPDe <- length(intersect(
  CPD_SNPlist$SNP, mvdat_ever1_5e6$SNP
)) + length(
  intersect(MVMRproxies_cpd_e$RS_Number, mvdat_ever1_5e6$SNP)
)

nSNPs_Cplus3HCn_mr <- mr_n_Cplus3HC_5e6$nsnp
nSNPs_CPDn_mr <- mr_n_cpd$nsnp
nSNPs_Cplus3HCe_mr <- mr_e_Cplus3HC_5e6$nsnp
nSNPs_CPDe_mr <- mr_e_cpd$nsnp

##### Betas #####
# Create dataframe including the effect estimates, lower and upper confidence intervals

##### Ever smokers #####
setwd("MVMR_Cancer")

cochrane_from_rmeta_ever <-
  structure(
    list(
      mean  = c(NA, all_results$OR[9:16]),
      lower = c(NA, all_results$LCI[9:16]),
      upper = c(NA, all_results$UCI[9:16])
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -9L),
    class = "data.frame"
  )

OR_LCI_UCI <- str_c(
  round(all_results$OR, digits = 2), ", (",
  round(all_results$LCI, digits = 2), ", ",
  round(all_results$UCI, digits = 2), ")"
)

tabletext_ever <- cbind(
  c(
    "Exposure",
    "3HC plus Cotinine", NA, NA, NA,
    "Smoking Heaviness", NA, NA, NA
  ),
  c(
    "N SNPs",
    nSNPs_Cplus3HCe_mr,
    nSNPs_Cplus3HCe_5e6,
    nSNPs_Cplus3HCe_mr,
    nSNPs_Cplus3HCe_5e6,
    nSNPs_CPDe_mr,
    nSNPs_CPDe,
    nSNPs_CPDe_mr,
    nSNPs_CPDe
  ),
  c(
    "Method",
    "MR IVW", "MVMR IVW", "MR-Egger (SIMEX)", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[9:16]),
  c("P value", round(all_results$p[9:16], digits = 3))
)
tabletext_ever[, 5][tabletext_ever[, 5] == 0] <- "<0.001"

##### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_ever_Cplus3HC.pdf", width = 15, height = 15)
forestplot(tabletext_ever,
  cochrane_from_rmeta_ever,
  new_page = TRUE,
  is.summary = c(TRUE, rep(FALSE, 9)),
  lineheight = unit(1, "cm"),
  graphwidth = unit(10, "cm"),
  boxsize = 0.15,
  clip = c(0.1, 20),
  xticks = c(0.5, 1, 1.5, 3, 5, 10, 20),
  xlog = TRUE,
  zero = 1,
  ci.vertices = TRUE,
  colgap = unit(4, "mm"),
  cex = 0.5,
  txt_gp = fpTxtGp(ticks = gpar(cex = 1)),
  col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue")
)
dev.off()


##### Never smokers #####
setwd("MVMR_Cancer")

cochrane_from_rmeta_never <-
  structure(
    list(
      mean  = c(NA, all_results$OR[1:8]),
      lower = c(NA, all_results$LCI[1:8]),
      upper = c(NA, all_results$UCI[1:8])
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -9L),
    class = "data.frame"
  )

OR_LCI_UCI <- str_c(
  round(all_results$OR, digits = 2), " (",
  round(all_results$LCI, digits = 2), ", ",
  round(all_results$UCI, digits = 2), ")"
)

tabletext_never <- cbind(
  c(
    "Exposure",
    "3HC plus Cotinine", NA, NA, NA,
    "Smoking Heaviness", NA, NA, NA
  ),
  c(
    "N SNPs",
    nSNPs_Cplus3HCn_mr,
    nSNPs_Cplus3HCn_5e6,
    nSNPs_Cplus3HCn_mr,
    nSNPs_Cplus3HCn_5e6,
    nSNPs_CPDn_mr,
    nSNPs_CPDn,
    nSNPs_CPDn_mr,
    nSNPs_CPDn
  ),
  c(
    "Method",
    "MR IVW", "MVMR IVW", "MR-Egger (SIMEX)", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[1:8]),
  c("P value", round(all_results$p[1:8], digits = 3))
)
tabletext_never[, 5][tabletext_never[, 5] == 0] <- "<0.001"

###### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_never_Cplus3HC.pdf", width = 15, height = 15)
forestplot(tabletext_never,
  cochrane_from_rmeta_never,
  new_page = TRUE,
  is.summary = c(TRUE, rep(FALSE, 9)),
  lineheight = unit(1, "cm"),
  graphwidth = unit(10, "cm"),
  boxsize = 0.15,
  clip = c(0.1, 20),
  xticks = c(0.1, 0.5, 1, 1.5, 3, 5, 10, 20),
  xlog = TRUE, zero = 1,
  ci.vertices = TRUE,
  colgap = unit(4, "mm"),
  cex = 0.5,
  txt_gp = fpTxtGp(ticks = gpar(cex = 1)),
  col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue")
)
dev.off()
