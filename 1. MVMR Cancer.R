# This script was created by Jasmine Khouja 27.04.22.

# The script conducts univariable and multivariable MR exploring the effects of
# nicotine and non-nicotine constituents of tobacco smoke (measured by nicotine
# metabolite ratio [NMR] and cigarettes per day [CPD]) on cancer.

# The SNPs used in this script were selected based on the finding from GSCAN
# (Liu et al., 2019) and Buchwald and colleagues.

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
##### Set SNP lists for NMR, CPD and both #####
################################################################################

NMR_SNPlist <- read.table("NMR_SNPlist.txt", header = FALSE)
NMR_SNPlist <- rename(NMR_SNPlist, c("SNP" = "V1"))
NMR_SNPlist_fm <- read.table("NMR_SNPlist_putatively causal.txt", header = TRUE)

CPD_SNPlist <- read.table("CPD_SNPlist.txt", header = FALSE)
CPD_SNPlist <- rename(CPD_SNPlist, c("SNP" = "V1"))

SNPlist <- read.table("SNPlist.txt", header = FALSE)
SNPlist <- rename(SNPlist, c("SNP" = "V1"))

################################################################################
##### Extract exposure data for MR of NMR and lung cancer#####
################################################################################

nmr_dat <- read_exposure_data("META_nmr_Extended_r.ma",
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

nmr_dat_mr <- format_data(
  nmr_dat,
  type = "exposure",
  snps = NMR_SNPlist$SNP,
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
  nmr_dat_mr,
  "nmr_dat_mr.csv",
  row.names = FALSE
)
write.csv(
  nmr_dat,
  "nmr_dat.csv",
  row.names = FALSE
)

################################################################################
##### Extract exposure data for MR of CPD and lung cancer #####
################################################################################

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
outcome_dat_never_nmr <- read_outcome_data(
  "Never.txt",
  snps = nmr_dat_mr$SNP,
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

outcome_dat_ever_nmr <- read_outcome_data(
  "Ever.txt",
  snps = nmr_dat_mr$SNP,
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

outcome_dat_never_nmr$beta.outcome <- as.numeric(as.character(
  outcome_dat_never_nmr$beta.outcome
))
outcome_dat_ever_nmr$beta.outcome <- as.numeric(as.character(
  outcome_dat_ever_nmr$beta.outcome
))
outcome_dat_never_cpd$beta.outcome <- as.numeric(as.character(
  outcome_dat_never_cpd$beta.outcome
))
outcome_dat_ever_cpd$beta.outcome <- as.numeric(as.character(
  outcome_dat_ever_cpd$beta.outcome
))

outcome_dat_never_nmr["beta.outcome"] <- log(
  outcome_dat_never_nmr["beta.outcome"]
)
outcome_dat_ever_nmr["beta.outcome"] <- log(
  outcome_dat_ever_nmr["beta.outcome"]
)
outcome_dat_never_cpd["beta.outcome"] <- log(
  outcome_dat_never_cpd["beta.outcome"]
)
outcome_dat_ever_cpd["beta.outcome"] <- log(
  outcome_dat_ever_cpd["beta.outcome"]
)

# Note: Outcome se is already se for logodds rather than odds:
# Checked using
# tail <- 2
#> se <- abs(outcome_dat_never_nmr$beta.outcome[1]/ qnorm(outcome_dat_never_nmr$pval.outcome[1]/tail))
#> outcome_dat_never_nmr$se[1]

################################################################################
##### Harmonising #####
################################################################################

dat_nmr_never <- harmonise_data(nmr_dat_mr, outcome_dat_never_nmr, action = 2)
dat_nmr_ever <- harmonise_data(nmr_dat_mr, outcome_dat_ever_nmr, action = 2)
dat_cpd_never <- harmonise_data(cpd_dat_mr, outcome_dat_never_cpd, action = 2)
dat_cpd_ever <- harmonise_data(cpd_dat_mr, outcome_dat_ever_cpd, action = 2)

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed <- data.frame(setdiff(NMR_SNPlist$SNP, dat_nmr_never$SNP))
proxy_needed1 <- data.frame(setdiff(NMR_SNPlist$SNP, dat_nmr_ever$SNP))
proxy_needed2 <- data.frame(setdiff(CPD_SNPlist$SNP, dat_cpd_never$SNP))
proxy_needed3 <- data.frame(setdiff(CPD_SNPlist$SNP, dat_cpd_ever$SNP))

################################################################################
##### Run to find proxies #####
################################################################################

source("proxy_search_loop_CancerMVMR.R")

################################################################################
##### Re-read outcome data with proxies #####
# Note: Because the SNPs that were missing in the outcome have now been replaced 
# with proxy SNPs in the exposure, need to reload the outcome data looking for 
# the proxy SNPs.
################################################################################

setwd("")
memory.limit(size = 80000)
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

memory.limit(size = 80000)
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

outcome_dat_never_cpd$beta.outcome <- as.numeric(
  as.character(outcome_dat_never_cpd$beta.outcome)
)
outcome_dat_ever_cpd$beta.outcome <- as.numeric(
  as.character(outcome_dat_ever_cpd$beta.outcome)
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

for (i in 1:length(dat_nmr_never$beta.exposure)) {
  if (dat_nmr_never$beta.exposure[i] < 0) {
    dat_nmr_never$beta.outcome[i] <- -1 * dat_nmr_never$beta.outcome[i]
  }
  if (dat_nmr_never$beta.exposure[i] < 0) {
    dat_nmr_never$beta.exposure[i] <- -1 * dat_nmr_never$beta.exposure[i]
  }
}
for (i in 1:length(dat_nmr_ever$beta.exposure)) {
  if (dat_nmr_ever$beta.exposure[i] < 0) {
    dat_nmr_ever$beta.outcome[i] <- -1 * dat_nmr_ever$beta.outcome[i]
  }
  if (dat_nmr_ever$beta.exposure[i] < 0) {
    dat_nmr_ever$beta.exposure[i] <- -1 * dat_nmr_ever$beta.exposure[i]
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

write.csv(dat_nmr_never, "dat_nmr_never.csv", row.names = FALSE)
write.csv(dat_nmr_ever, "dat_nmr_ever.csv", row.names = FALSE)
write.csv(dat_cpd_never, "dat_cpd_never.csv", row.names = FALSE)
write.csv(dat_cpd_ever, "dat_cpd_ever.csv", row.names = FALSE)
write.csv(nmr_dat_mr, "nmr_dat_mr.csv", row.names = FALSE)
write.csv(cpd_dat_mr, "cpd_dat_mr.csv", row.names = FALSE)

################################################################################
##### Removing rs117090198  #####
# Note: Removed because the SNP has a very high p-value prior to the conditional 
# independence analysis which is unusual and unexplained
################################################################################

dat_nmr_never <- dat_nmr_never[!dat_nmr_never$SNP == "rs117090198", ]
dat_nmr_ever <- dat_nmr_ever[!dat_nmr_ever$SNP == "rs117090198", ]

################################################################################
################################## MR ##########################################
################################################################################

################################################################################
##### Generate results inc. F and Q stats for heterogeneity for: #####
# NMR never
# NMR ever
# CPD never 
# CPD ever
################################################################################

# NMR never
result_nmr_never <- mr(
  dat_nmr_never,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_nmr_never <- generate_odds_ratios(result_nmr_never)
dat_nmr_never <- subset(dat_nmr_never, mr_keep)
mr_n_nmr <- TwoSampleMR::mr_ivw(
  dat_nmr_never$beta.exposure,
  dat_nmr_never$beta.outcome,
  dat_nmr_never$se.exposure,
  dat_nmr_never$se.outcome,
  parameters = default_parameters()
)
ptr_n_nmr <- data.frame(mr_n_nmr["Q"])
egger_n_nmr <- mr_egger_regression(
  dat_nmr_never$beta.exposure,
  dat_nmr_never$beta.outcome,
  dat_nmr_never$se.exposure,
  dat_nmr_never$se.outcome,
  parameters
)
ptr_n_nmr[2, 1] <- egger_n_nmr["Q"]
F <- abs(dat_nmr_never$beta.exposure)^2 / dat_nmr_never$se.exposure^2
mF <- mean(F)
ptr_n_nmr[3, 1] <- mF

# NMR ever
result_nmr_ever <- mr(
  dat_nmr_ever,
  parameters = default_parameters(),
  method_list = subset(mr_method_list(), use_by_default)$obj
)
result_nmr_ever <- generate_odds_ratios(result_nmr_ever)
dat_nmr_ever <- subset(dat_nmr_ever, mr_keep)
mr_e_nmr <- TwoSampleMR::mr_ivw(
  dat_nmr_ever$beta.exposure,
  dat_nmr_ever$beta.outcome,
  dat_nmr_ever$se.exposure,
  dat_nmr_ever$se.outcome,
  parameters = default_parameters()
)
ptr_e_nmr <- data.frame(mr_e_nmr["Q"])
egger_e_nmr <- mr_egger_regression(
  dat_nmr_ever$beta.exposure,
  dat_nmr_ever$beta.outcome,
  dat_nmr_ever$se.exposure,
  dat_nmr_ever$se.outcome,
  parameters
)
ptr_e_nmr[2, 1] <- egger_e_nmr["Q"]
F <- dat_nmr_ever$beta.exposure^2 / dat_nmr_ever$se.exposure^2
mF <- mean(F)
ptr_e_nmr[3, 1] <- mF

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
# This indicates the effect estimate is biased by 10% due to measurement error
# (Bowden, Del Greco, et al., 2016)
# (Lederer & K?chenhoff, 2006)
################################################################################

ISQ <- data.frame(c(1:4))
ISQ[1, 1] <- Isq(dat_nmr_never$beta.exposure, dat_nmr_never$se.exposure)
ISQ[2, 1] <- Isq(dat_nmr_ever$beta.exposure, dat_nmr_ever$se.exposure)
ISQ[3, 1] <- Isq(dat_cpd_never$beta.exposure, dat_cpd_never$se.exposure)
ISQ[4, 1] <- Isq(dat_cpd_ever$beta.exposure, dat_cpd_ever$se.exposure)

################################################################################
##### Forest plot results #####
################################################################################

forest_nmr_never <- mr_forest(
  mr_input(
    bx = dat_nmr_never$beta.exposure,
    bxse = dat_nmr_never$se.exposure,
    by = dat_nmr_never$beta.outcome,
    byse = dat_nmr_never$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)
forest_nmr_ever <- mr_forest(
  mr_input(
    bx = dat_nmr_ever$beta.exposure,
    bxse = dat_nmr_ever$se.exposure,
    by = dat_nmr_ever$beta.outcome,
    byse = dat_nmr_ever$se.outcome
  ),
  methods = c("ivw", "wmedian", "egger"), snp_estimates = FALSE
)
forest <- mr_forest(
  mr_input(
    bx = dat_nmr_ever$beta.exposure,
    bxse = dat_nmr_ever$se.exposure,
    by = dat_nmr_ever$beta.outcome,
    byse = dat_nmr_ever$se.outcome
  )
)
forest2 <- forest + coord_cartesian(xlim = c(-2, 2))
forest2
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
  bx = dat_nmr_ever$beta.exposure,
  bxse = dat_nmr_ever$se.exposure,
  by = dat_nmr_ever$beta.outcome,
  byse = dat_nmr_ever$se.outcome
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
# Including removing rs117090198
# Load data if not already loaded
################################################################################
NMR_instr <- nmr_dat_mr[, c(
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
NMR_instr <- NMR_instr[!NMR_instr$SNP == "rs117090198", ]
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

INSTR <- do.call("rbind", list(NMR_instr, CPD_instr))

################################################################################
##### Check for overlapping SNPs between the exposure instruments #####
################################################################################

n_occur <- data.frame(table(INSTR$SNP))
n_occur[n_occur$Freq > 1, ]
INSTR[INSTR$SNP %in% n_occur$Var1[n_occur$Freq > 1], ]

################################################################################
##### Extract instruments from the NMR data set #####
################################################################################

nmr_dat_mvmr <- format_data(
  nmr_dat,
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

##### Change name of GWAS and check n SNPs ####

nmr_dat_mvmr$id.exposure <- "1"
str(nmr_dat_mvmr)

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
##### Merge them to extract from outcome #####
################################################################################

nmr_dat_mvmr_1 <- subset(nmr_dat_mvmr, select = c(
  "SNP", "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure", "se.exposure",
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

##### check structure is the same #####
str(nmr_dat_mvmr_1)
str(cpd_dat_mvmr_1)

################################################################################
##### Merge ####
################################################################################

exposures <- do.call("rbind", list(nmr_dat_mvmr_1, cpd_dat_mvmr_1))

##### Save dataframe #####
write.csv(
  exposures,
  "NMR_CPD.csv",
  row.names = FALSE
)

################################################################################
##### Find proxies missing from either the NMR or CPD data set #####
# Note: rs4886550 is not available in the ref panel and is the only missing SNP
################################################################################

n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq < 2, ]
exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq < 2], ]

################################################################################
##### Clumping ####
################################################################################

##### Change all p-values for NMR to 1e-200 for clumping so that none are dropped #####
##### Save old p-values first #####
exposures$oldpvalues <- exposures$pval.exposure
exposures <- exposures %>%
  mutate(pval.exposure = if_else(
    exposures$SNP %in% NMR_instr$SNP,
    1e-201, pval.exposure
  ))

##### Clump the data #####
exposures$id.exposure[exposures$id.exposure == "2"] <- "1"
exposures <- clump_data(exposures, clump_kb = 500, clump_r2 = 0.1)
str(exposures)

##### Add ID's back #####
exposures$id.exposure[exposures$samplesize.exposure < 6000] <- "1"
exposures$id.exposure[exposures$samplesize.exposure > 6000] <- "2"

##### Revert all p-values for NMR from 1e-200 #####
exposures$pval.exposure <- exposures$oldpvalues
exposures <- select(exposures, -c(oldpvalues))

##### Added rs56113850 back in #####
##### Removed due to LD with other NMR SNPs but conditionally independent #####
rs56113850 <- do.call("rbind", list(nmr_dat_mvmr_1, cpd_dat_mvmr_1))
rs56113850 <- rs56113850[rs56113850$SNP == "rs56113850", ]
exposures <- do.call("rbind", list(exposures, rs56113850))

##### Split again to harmonise based on exposure id #####
NMR <- split(exposures, exposures$id.exposure)[["1"]]
CPD <- split(exposures, exposures$id.exposure)[["2"]]

################################################################################
##### Harmonise NMR on CPD #####
################################################################################

names(NMR) <- gsub("exposure", "outcome", names(NMR))
NMR_CPD <- harmonise_data(CPD, NMR)

################################################################################
##### Keep only snps that are present across both exposures ####
# Note: they would have frequency 1 if only available in one dataset
################################################################################

n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq == 2, ]
exposures <- exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq == 2], ]
str(exposures)

################################################################################
##### Format exposures #####
################################################################################

##### Keep onlySNPs where MrKeep = TRUE #####
NMR_CPD <- NMR_CPD[NMR_CPD$mr_keep == TRUE, ]
str(NMR_CPD)

##### Split the tables - CPD #####
CPD_H <- subset(
  NMR_CPD,
  id.exposure == "2",
  select = c(
    SNP, exposure,
    id.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    beta.exposure,
    se.exposure,
    pval.exposure,
    eaf.exposure
  )
)

##### Split the tables - NMR #####
NMR_H <- subset(NMR_CPD,
  id.outcome == "1",
  select = c(
    SNP,
    outcome,
    id.outcome,
    effect_allele.outcome,
    other_allele.outcome,
    beta.outcome,
    se.outcome,
    pval.outcome,
    eaf.outcome
  )
)

##### Turn NMR from outcome to exposure to merge the data sets #####
names(NMR_H) <- gsub("outcome", "exposure", names(NMR_H))
Exposures_H <- merge(NMR_H, CPD_H, all = TRUE)
Exposures_H["Phenotype"] <- NA
Exposures_H$Phenotype[Exposures_H$id.exposure == 1] <- "NMR"
Exposures_H$Phenotype[Exposures_H$id.exposure == 2] <- "CPD"
str(Exposures_H)

################################################################################
##### Extract outcome data for MVMR #####
################################################################################

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

################################################################################
##### Convert odds ratios to log odds #####
################################################################################
outcome_dat_never$beta.outcome <- as.numeric(
  as.character(
    outcome_dat_never$beta.outcome
  )
)
outcome_dat_ever$beta.outcome <- as.numeric(
  as.character(
    outcome_dat_ever$beta.outcome
  )
)
outcome_dat_never["beta.outcome"] <- log(outcome_dat_never["beta.outcome"])
outcome_dat_ever["beta.outcome"] <- log(outcome_dat_ever["beta.outcome"])

################################################################################
##### Organise outcome #####
################################################################################

outcome_dat_never["Phenotype"] <- NA
outcome_dat_never$Phenotype <- "LungCancer"

outcome_dat_ever["Phenotype"] <- NA
outcome_dat_ever$Phenotype <- "LungCancer"

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

mvdat_never1 <- harmonise_data(Exposures_H, outcome_dat_never)
mvdat_never1 <- mvdat_never1[mvdat_never1$mr_keep == TRUE, ]
str(mvdat_never1)

mvdat_ever1 <- harmonise_data(Exposures_H, outcome_dat_ever)
mvdat_ever1 <- mvdat_ever1[mvdat_ever1$mr_keep == TRUE, ]
str(mvdat_ever1)

################################################################################
##### Find proxies to add to missing outcome SNPs #####
################################################################################

proxy_needed4 <- data.frame(setdiff(NMR_H$SNP, outcome_dat_never$SNP))
proxy_needed5 <- data.frame(setdiff(NMR_H$SNP, outcome_dat_ever$SNP))

source("proxy_search_loop_CancerMVMR2.R")

################################################################################
##### Re-merge the data #####
################################################################################

Exposures_H_never <- merge(NMR_never, CPD_never, all = TRUE)
Exposures_H_never["Phenotype"] <- NA
Exposures_H_never$Phenotype[Exposures_H_never$id.exposure == 1] <- "NMR"
Exposures_H_never$Phenotype[Exposures_H_never$id.exposure == 2] <- "CPD"
str(Exposures_H_never)


Exposures_H_ever <- merge(NMR_ever, CPD_ever, all = TRUE)
Exposures_H_ever["Phenotype"] <- NA
Exposures_H_ever$Phenotype[Exposures_H_ever$id.exposure == 1] <- "NMR"
Exposures_H_ever$Phenotype[Exposures_H_ever$id.exposure == 2] <- "CPD"
str(Exposures_H_ever)

################################################################################
##### Re-read the outcome data #####
################################################################################

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

################################################################################
##### Harmonise data #####
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

################################################################################
##### Save dataframes #####
################################################################################

write.csv(
  mvdat_never1,
  "mvdat_never.csv",
  row.names = FALSE
)
write.csv(
  mvdat_ever1,
  "mvdat_ever.csv",
  row.names = FALSE
)

##### Can run from here if data frames are unchanged.#####
##### Can skip if data is already loaded. #####
setwd("")
mvdat_never1 <- read.csv(
  "mvdat_never.csv",
  header = TRUE
)
mvdat_ever1 <- read.csv(
  "mvdat_ever.csv",
  header = TRUE
)

################################################################################
##### Run MVMR #####
################################################################################

##### IVW never #####

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
  min(summary(lm(bY ~ bX1 + bX2 - 1,
    weights = bYse^-2
  ))$sigma, 1)

mod_n <- summary(mod.MVMR_never)

mod_n_or <- coef(summary(mod.MVMR_never))
colnames(mod_n_or) <- c("b", "se", "t", "p")
mod_n_or <- as.data.frame(mod_n_or)
mod_n_or <- generate_odds_ratios(mod_n_or)

##### Orientation NMR #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for NMR, and then          #####
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
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_n_2 <- summary(mod.MVMRME_never_2)

mod_ME_n_2_or <- data.frame(mod.MVMRME_never_2[["coefficients"]])
colnames(mod_ME_n_2_or) <- c("b", "se", "t", "p")
mod_ME_n_2_or <- as.data.frame(mod_ME_n_2_or)
mod_ME_n_2_or <- generate_odds_ratios(mod_ME_n_2_or)

##### IVW ever #####
bX1 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 1])
bX2 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 2])
bY <- c(mvdat_ever1$beta.outcome[mvdat_ever1$id.exposure == 1])
bYse <- c(mvdat_ever1$se.outcome[mvdat_ever1$id.exposure == 1])

set.seed(1234)
mod.MVMR_ever <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2))$sigma, 1)

mod_e <- summary(mod.MVMR_ever)

mod_e_or <- coef(summary(mod.MVMR_ever))
colnames(mod_e_or) <- c("b", "se", "t", "p")
mod_e_or <- as.data.frame(mod_e_or)
mod_e_or <- generate_odds_ratios(mod_e_or)
##### Orientation NMR #####
clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_ever <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
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
se_theta1ME.random <- summary(lm
(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_e_2 <- summary(mod.MVMRME_ever_2)


mod_ME_e_2_or <- data.frame(mod.MVMRME_ever_2[["coefficients"]])
colnames(mod_ME_e_2_or) <- c("b", "se", "t", "p")
mod_ME_e_2_or <- as.data.frame(mod_ME_e_2_or)
mod_ME_e_2_or <- generate_odds_ratios(mod_ME_e_2_or)

##### Format to analyse / cross check with MVMR package #####
bX1 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 1])
bX2 <- c(mvdat_never1$beta.exposure[mvdat_never1$id.exposure == 2])
bY <- c(mvdat_never1$beta.outcome[mvdat_never1$id.exposure == 1])
bYse <- c(mvdat_never1$se.outcome[mvdat_never1$id.exposure == 1])
bXse1 <- c(mvdat_never1$se.exposure[mvdat_never1$id.exposure == 1])
bXse2 <- c(mvdat_never1$se.exposure[mvdat_never1$id.exposure == 2])
df_never <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

bX1 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 1])
bX2 <- c(mvdat_ever1$beta.exposure[mvdat_ever1$id.exposure == 2])
bY <- c(mvdat_ever1$beta.outcome[mvdat_ever1$id.exposure == 1])
bYse <- c(mvdat_ever1$se.outcome[mvdat_ever1$id.exposure == 1])
bXse1 <- c(mvdat_ever1$se.exposure[mvdat_ever1$id.exposure == 1])
bXse2 <- c(mvdat_ever1$se.exposure[mvdat_ever1$id.exposure == 2])
df_ever <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

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
##### Cross check result with MVMR package #####
res_n <- ivw_mvmr(df_mvmr_never)
res_e <- ivw_mvmr(df_mvmr_ever)

################################################################################
##### Calculate F-statistic and covariance #####
# Note: >10 is strong
# correlation between NMR and CPD in Buchwald  = -0.019
# qhet_mvmr is used to adjust for covariance
# At present, CIs take substantial time to calculate
# Compare effects with and without adjustment
################################################################################

cov <- matrix(c(1, -0.019, -0.019, 1), nrow = 2, ncol = 2)

Xcovmat_n <- phenocov_mvmr(cov, df_mvmr_never[, c(6, 7)])
Fstat_n <- strength_mvmr(df_mvmr_never, gencov = Xcovmat_n)

Xcovmat_e <- phenocov_mvmr(cov, df_mvmr_ever[, c(6, 7)])
Fstat_e <- strength_mvmr(df_mvmr_ever, gencov = Xcovmat_e)

cov_adj_mvmr_n <- qhet_mvmr(df_mvmr_never, cov, CI = FALSE, iterations = 1000)

cov_adj_mvmr_e <- qhet_mvmr(df_mvmr_ever, cov, CI = FALSE, iterations = 1000)

################################################################################
##### Test for horizontal pleiotropy #####
# Note: Q should be greater than the number of SNPs included
################################################################################

ptr_n <- pleiotropy_mvmr(df_mvmr_never, gencov = Xcovmat_n)
ptr_e <- pleiotropy_mvmr(df_mvmr_ever, gencov = Xcovmat_e)

################################################################################
##### Forest Plots #####
################################################################################

Exposure <- c(
  "NMR", "NMR", "NMR", "NMR",
  "CPD", "CPD", "CPD", "CPD",
  "NMR", "NMR", "NMR", "NMR",
  "CPD", "CPD", "CPD", "CPD"
)
Smoking <- c(
  "Never", "Never", "Never", "Never", "Never", "Never", "Never", "Never",
  "Ever", "Ever", "Ever", "Ever", "Ever", "Ever", "Ever", "Ever"
)
res_n <- data.frame(res_n)
res_e <- data.frame(res_e)
Method <- c(
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger",
  "MR-IVW", "MVMR-IVW", "MR-Egger", "MVMR-Egger"
)
OR <- c(
  result_nmr_never[3, "or"], mod_n_or[1, 7],
  result_nmr_never[1, "or"], mod_ME_n_or[2, 7],
  result_cpd_never[3, "or"], mod_n_or[2, 7],
  result_cpd_never[1, "or"], mod_ME_n_2_or[3, 7],
  result_nmr_ever[3, "or"], mod_e_or[1, 7],
  result_nmr_ever[1, "or"], mod_ME_e_or[2, 7],
  result_cpd_ever[3, "or"], mod_e_or[2, 7],
  result_cpd_ever[1, "or"], mod_ME_e_2_or[3, 7]
)

LCI <- c(
  result_nmr_never[3, "or_lci95"], mod_n_or[1, 8],
  result_nmr_never[1, "or_lci95"], mod_ME_n_or[2, 8],
  result_cpd_never[3, "or_lci95"], mod_n_or[2, 8],
  result_cpd_never[1, "or_lci95"], mod_ME_n_2_or[3, 8],
  result_nmr_ever[3, "or_lci95"], mod_e_or[1, 8],
  result_nmr_ever[1, "or_lci95"], mod_ME_e_or[2, 8],
  result_cpd_ever[3, "or_lci95"], mod_e_or[2, 8],
  result_cpd_ever[1, "or_lci95"], mod_ME_e_2_or[3, 8]
)

UCI <- c(
  result_nmr_never[3, "or_uci95"], mod_n_or[1, 9],
  result_nmr_never[1, "or_uci95"], mod_ME_n_or[2, 9],
  result_cpd_never[3, "or_uci95"], mod_n_or[2, 9],
  result_cpd_never[1, "or_uci95"], mod_ME_n_2_or[3, 9],
  result_nmr_ever[3, "or_uci95"], mod_e_or[1, 9],
  result_nmr_ever[1, "or_uci95"], mod_ME_e_or[2, 9],
  result_cpd_ever[3, "or_uci95"], mod_e_or[2, 9],
  result_cpd_ever[1, "or_uci95"], mod_ME_e_2_or[3, 9]
)

p <- c(
  result_nmr_never[3, "pval"], mod_n_or[1, 4],
  result_nmr_never[1, "pval"], mod_ME_n_or[2, 4],
  result_cpd_never[3, "pval"], mod_n_or[2, 4],
  result_cpd_never[1, "pval"], mod_ME_n_2_or[3, 4],
  result_nmr_ever[3, "pval"], mod_e_or[1, 4],
  result_nmr_ever[1, "pval"], mod_ME_e_or[2, 4],
  result_cpd_ever[3, "pval"], mod_e_or[2, 4],
  result_cpd_ever[1, "pval"], mod_ME_e_2_or[3, 4]
)

I2Gx <- c(
  ".", ".", ISQ[1, 1], ".",
  ".", ".", ISQ[3, 1], ".",
  ".", ".", ISQ[2, 1], ".",
  ".", ".", ISQ[4, 1], "."
)
Q <- c(
  ptr_n_nmr[1, 1], ptr_n[["Qstat"]], ptr_n_nmr[2, 1], ptr_n[["Qstat"]],
  ptr_n_cpd[1, 1], ptr_n[["Qstat"]], ptr_n_cpd[2, 1], ptr_n[["Qstat"]],
  ptr_e_nmr[1, 1], ptr_e[["Qstat"]], ptr_e_nmr[2, 1], ptr_e[["Qstat"]],
  ptr_e_cpd[1, 1], ptr_e[["Qstat"]], ptr_e_cpd[2, 1], ptr_e[["Qstat"]]
)

EggerI <- c(
  ".", ".", egger_n_nmr[["b_i"]], mod_ME_n_or[1, 7],
  ".", ".", egger_n_cpd[["b_i"]], mod_ME_n_2_or[1, 7],
  ".", ".", egger_e_nmr[["b_i"]], mod_ME_e_or[1, 7],
  ".", ".", egger_e_cpd[["b_i"]], mod_ME_e_2_or[1, 7]
)

EggerIp <- c(
  ".", ".", egger_n_nmr[["pval_i"]], mod_ME_n_or[1, 4],
  ".", ".", egger_n_cpd[["pval_i"]], mod_ME_n_2_or[1, 4],
  ".", ".", egger_e_nmr[["pval_i"]], mod_ME_e_or[1, 4],
  ".", ".", egger_e_cpd[["pval_i"]], mod_ME_e_2_or[1, 4]
)

F_stat <- c(
  ptr_n_nmr[3, 1], Fstat_n[1, 1], ptr_n_nmr[3, 1], Fstat_n[1, 1],
  ptr_n_cpd[3, 1], Fstat_n[1, 2], ptr_n_cpd[3, 1], Fstat_n[1, 2],
  ptr_e_nmr[3, 1], Fstat_e[1, 1], ptr_e_nmr[3, 1], Fstat_e[1, 1],
  ptr_e_cpd[3, 1], Fstat_e[1, 2], ptr_e_cpd[3, 1], Fstat_e[1, 2]
)

all_results <- data.frame(
  Exposure, Smoking, Method,
  OR, LCI, UCI, p,
  I2Gx, Q, EggerI, EggerIp, F_stat
)

write.csv(all_results,
  "Results_NMR_CPD_lung.csv",
  row.names = FALSE
)

MVMRproxies_nmr_n <- subset(SNPs.proxies4, target %in% NMR_SNPlist$SNP)
MVMRproxies_cpd_n <- subset(SNPs.proxies4, target %in% CPD_SNPlist$SNP)
MVMRproxies_nmr_e <- subset(SNPs.proxies5, target %in% NMR_SNPlist$SNP)
MVMRproxies_cpd_e <- subset(SNPs.proxies5, target %in% CPD_SNPlist$SNP)

nSNPs_NMRn <- length(intersect(
  NMR_SNPlist$SNP,
  mvdat_never1$SNP
)) + length(intersect(
  MVMRproxies_nmr_n$RS_Number,
  mvdat_never1$SNP
))
nSNPs_CPDn <- length(intersect(
  CPD_SNPlist$SNP, mvdat_never1$SNP
)) + length(intersect(
  MVMRproxies_cpd_n$RS_Number, mvdat_never1$SNP
))
nSNPs_NMRe <- length(intersect(
  NMR_SNPlist$SNP, mvdat_ever1$SNP
)) + length(intersect(
  MVMRproxies_nmr_e$RS_Number, mvdat_ever1$SNP
))
nSNPs_CPDe <- length(intersect(
  CPD_SNPlist$SNP, mvdat_ever1$SNP
)) + length(intersect(
  MVMRproxies_cpd_e$RS_Number, mvdat_ever1$SNP
))

nSNPs_NMRn_mr <- mr_n_nmr$nsnp
nSNPs_CPDn_mr <- mr_n_cpd$nsnp
nSNPs_NMRe_mr <- mr_e_nmr$nsnp
nSNPs_CPDe_mr <- mr_e_cpd$nsnp

##### Betas #####
# Create dataframe of the effect estimates, lower and upper confidence intervals
##### Ever smokers #####

setwd("")

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
  round(all_results$OR, digits = 2),
  " (",
  round(all_results$LCI, digits = 2),
  ", ",
  round(all_results$UCI, digits = 2),
  ")"
)

tabletext_ever <- cbind(
  c("Exposure", "NMR", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
  c(
    "N SNPs", nSNPs_NMRe_mr, nSNPs_NMRe, nSNPs_NMRe_mr, nSNPs_NMRe,
    nSNPs_CPDe_mr, nSNPs_CPDe, nSNPs_CPDe_mr, nSNPs_CPDe
  ),
  c(
    "Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[9:16]),
  c("P value", round(all_results$p[9:16], digits = 3))
)
tabletext_ever[, 5][tabletext_ever[, 5] == 0] <- "<0.001"

##### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_ever.pdf", width = 15, height = 15)
forestplot(tabletext_ever,
  cochrane_from_rmeta_ever,
  new_page = TRUE,
  is.summary = c(TRUE, rep(FALSE, 9)),
  lineheight = unit(1, "cm"),
  graphwidth = unit(10, "cm"),
  boxsize = 0.15,
  clip = c(0.1, 10),
  xticks = c(0.5, 1, 1.5, 3, 5, 10),
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
setwd("")

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
  c("Exposure", "NMR", NA, NA, NA, "Smoking Heaviness", NA, NA, NA),
  c(
    "N SNPs", nSNPs_NMRn_mr, nSNPs_NMRn, nSNPs_NMRn_mr, nSNPs_NMRn,
    nSNPs_CPDn_mr, nSNPs_CPDn, nSNPs_CPDn_mr, nSNPs_CPDn
  ),
  c(
    "Method", "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger",
    "MR IVW", "MVMR IVW", "MR-Egger", "MVMR-Egger"
  ),
  c("OR (95% CI)", OR_LCI_UCI[1:8]),
  c("P value", round(all_results$p[1:8], digits = 3))
)
tabletext_never[, 5][tabletext_never[, 5] == 0] <- "<0.001"

###### Create graphs and write to file #####
pdf.options(reset = TRUE, onefile = FALSE)
pdf("MVMR_never.pdf", width = 15, height = 15)
forestplot(tabletext_never,
  cochrane_from_rmeta_never,
  new_page = TRUE,
  is.summary = c(TRUE, rep(FALSE, 9)),
  lineheight = unit(1, "cm"),
  graphwidth = unit(10, "cm"),
  boxsize = 0.15,
  clip = c(0.1, 3),
  xticks = c(0.5, 1, 1.5, 3),
  xlog = TRUE, zero = 1,
  ci.vertices = TRUE,
  colgap = unit(4, "mm"),
  cex = 0.5,
  txt_gp = fpTxtGp(ticks = gpar(cex = 1)),
  col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue")
)
dev.off()
