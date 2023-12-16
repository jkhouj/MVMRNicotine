# Script to run MR Egger SIMEX. Created by Jasmine Khouja 25/07/2022.

simexegger<-c()

################################################################################
################################ Cplus3HC Never ################################
################################################################################

################################################################################
##### Format data #####
################################################################################

BetaXG <- dat_Cplus3HC_never$beta.exposure
BetaYG <- dat_Cplus3HC_never$beta.outcome
seBetaXG <- dat_Cplus3HC_never$se.exposure
seBetaYG <- dat_Cplus3HC_never$se.outcome

BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG)         

################################################################################
##### Egger Analyses #####
################################################################################

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

################################################################################
##### Simulation #####
################################################################################

# Simulation extrapolation 
mod.sim1 <- simex(
  Fit1,B=1000, 
  measurement.error = seBetaXG, 
  SIMEXvariable="BXG",
  fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(
  Fit2,B=1000, 
  measurement.error = seBetaXG, 
  SIMEXvariable="BXG",
  fitting.method ="quad",
  asymptotic="FALSE") 

mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

# extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

# convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("exposure", "b", "se", "pval") 
# (following the MRBase naming convention)
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

# extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

################################################################################
##### Save results #####
################################################################################

output<-cbind(
  "Cplus3HC_N", 
  beta1, 
  lci1, 
  uci1, 
  p1, 
  or1, 
  lcior1, 
  ucior1, 
  beta2, 
  lci2, 
  uci2, 
  p2, 
  or2, 
  lcior2, 
  ucior2)
simexegger_Cplus3HCn <-rbind(simexegger, output)
colnames(simexegger_Cplus3HCn) <- c(
  "Exposure",
  "beta_weighted",
  "lowerCI_weighted",
  "upperCI_weighted", 
  "p_weighted", 
  "OR_weighted",
  "lCIOR_weighted", 
  "uCIOR_weighted", 
  "beta_unweighted", 
  "lowerCI_unweighted", 
  "upperCI_unweighted", 
  "p_unweighted",
  "OR_unweighted",
  "lCIOR_unweighted", 
  "uCIOR_unweighted")

################################################################################
################################ Cplus3HC Ever #################################
################################################################################

################################################################################
##### Format data #####
################################################################################

BetaXG <- dat_Cplus3HC_ever$beta.exposure
BetaYG <- dat_Cplus3HC_ever$beta.outcome
seBetaXG <- dat_Cplus3HC_ever$se.exposure
seBetaYG <- dat_Cplus3HC_ever$se.outcome

BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG)         

################################################################################
##### Egger Analyses #####
################################################################################

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

################################################################################
##### Simulation #####
################################################################################

# Simulation extrapolation 
mod.sim1 <- simex(
  Fit1,B=1000, 
  measurement.error = seBetaXG, 
  SIMEXvariable="BXG",
  fitting.method ="quad",
  asymptotic="FALSE") 
mod.sim2 <- simex(
  Fit2,B=1000, 
  measurement.error = seBetaXG,
  SIMEXvariable="BXG",
  fitting.method ="quad",
  asymptotic="FALSE") 

mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

# extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

# convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("exposure", "b", "se", "pval") 
# (following the MRBase naming convention)
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

# extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

################################################################################
##### Save results #####
################################################################################

output<-cbind(
  "Cplus3HC_E",
  beta1,
  lci1, 
  uci1, 
  p1, 
  or1, 
  lcior1, 
  ucior1, 
  beta2, 
  lci2, 
  uci2, 
  p2, 
  or2, 
  lcior2, 
  ucior2)
simexegger_Cplus3HCe <-rbind(simexegger, output)
colnames(simexegger_Cplus3HCe) <- c(
  "Exposure", 
  "beta_weighted", 
  "lowerCI_weighted",
  "upperCI_weighted", 
  "p_weighted", 
  "OR_weighted", 
  "lCIOR_weighted", 
  "uCIOR_weighted", 
  "beta_unweighted",
  "lowerCI_unweighted",
  "upperCI_unweighted",
  "p_unweighted", 
  "OR_unweighted",
  "lCIOR_unweighted", 
  "uCIOR_unweighted")

################################################################################
################################## Cplus3HC Never ##############################
################################################################################
# Note: Will not run with only 2 available SNPs - refer to 5e6 results

################################################################################
################################ Cplus3HC 5e6 Never ############################
################################################################################

################################################################################
##### Format data #####
################################################################################

BetaXG <- dat_Cplus3HC_never_5e6$beta.exposure
BetaYG <- dat_Cplus3HC_never_5e6$beta.outcome
seBetaXG <- dat_Cplus3HC_never_5e6$se.exposure
seBetaYG <- dat_Cplus3HC_never_5e6$se.outcome

BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG)         

################################################################################
##### Egger Analyses #####
################################################################################

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

################################################################################
##### Simulation #####
################################################################################

# Simulation extrapolation 
mod.sim1 <- simex(
  Fit1,B=1000, 
  measurement.error = seBetaXG, 
  SIMEXvariable="BXG",
  fitting.method ="quad",
  asymptotic="FALSE") 
mod.sim2 <- simex(
  Fit2,B=1000, 
  measurement.error = seBetaXG,
  SIMEXvariable="BXG",
  fitting.method ="quad",
  asymptotic="FALSE") 

mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

# extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

# convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("exposure", "b", "se", "pval") 
# (following the MRBase naming convention)
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

# extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

################################################################################
##### Save results #####
################################################################################

output<-cbind(
  "Cplus3HC_N_5e6", 
  beta1, 
  lci1,
  uci1, 
  p1, 
  or1, 
  lcior1,
  ucior1, 
  beta2, 
  lci2, 
  uci2, 
  p2, 
  or2, 
  lcior2, 
  ucior2)
simexegger_Cplus3HCn5e6 <-rbind(simexegger, output)
colnames(simexegger_Cplus3HCn5e6) <- c(
  "Exposure", 
  "beta_weighted", 
  "lowerCI_weighted",
  "upperCI_weighted", 
  "p_weighted", 
  "OR_weighted", 
  "lCIOR_weighted", 
  "uCIOR_weighted", 
  "beta_unweighted",
  "lowerCI_unweighted",
  "upperCI_unweighted",
  "p_unweighted", 
  "OR_unweighted", 
  "lCIOR_unweighted", 
  "uCIOR_unweighted")

################################################################################
################################ Cplus3HC 5e6 Ever #############################
################################################################################

################################################################################
##### Format data #####
################################################################################

BetaXG <- dat_Cplus3HC_ever_5e6$beta.exposure
BetaYG <- dat_Cplus3HC_ever_5e6$beta.outcome
seBetaXG <- dat_Cplus3HC_ever_5e6$se.exposure
seBetaYG <- dat_Cplus3HC_ever_5e6$se.outcome

BYG <- BetaYG*sign(BetaXG)
BXG <- abs(BetaXG)         

################################################################################
##### Egger Analyses #####
################################################################################

# MR-Egger regression (weighted) 
Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

################################################################################
##### Simulation #####
################################################################################

# Simulation extrapolation 
mod.sim1 <- simex(
  Fit1,B=1000, 
  measurement.error = seBetaXG, 
  SIMEXvariable="BXG",
  fitting.method ="quad",
  asymptotic="FALSE") 
mod.sim2 <- simex(
  Fit2,B=1000,
  measurement.error = seBetaXG,
  SIMEXvariable="BXG",
  fitting.method ="quad",
  asymptotic="FALSE") 

mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)

# extract results in beta format
beta1<-mod1$coefficients$jackknife[2,1]
se1<-mod1$coefficients$jackknife[2,2]
p1<-mod1$coefficients$jackknife[2,4]
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

# convert to odds ratios for categorical outcomes
results1<-cbind("weighted", beta1, se1, p1)
results2<-cbind("unweighted", beta2, se2, p2)
results<-rbind(results1, results2)
colnames(results) <- c("exposure", "b", "se", "pval") 
# (following the MRBase naming convention)
results<-data.frame(results)
results$b<-as.numeric(as.character(results$b))
results$se<-as.numeric(as.character(results$se))
results$pval<-as.numeric(as.character(results$pval))

or<-generate_odds_ratios(results)

# extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci1<-or[1,5]
uci1<-or[1,6]
or2<-or[2,7]
lcior2<-or[2,8]
ucior2<-or[2,9]
lci2<-or[2,5]
uci2<-or[2,6]

################################################################################
##### Save results #####
################################################################################

output<-cbind(
  "Cplus3HC_E_5e6",
  beta1, 
  lci1,
  uci1, 
  p1, 
  or1,
  lcior1,
  ucior1,
  beta2, 
  lci2, 
  uci2, 
  p2, 
  or2,
  lcior2,
  ucior2)
simexegger_Cplus3HCe5e6 <-rbind(simexegger, output)
colnames(simexegger_Cplus3HCe5e6) <- c(
  "Exposure",
  "beta_weighted", 
  "lowerCI_weighted", 
  "upperCI_weighted",
  "p_weighted", 
  "OR_weighted", 
  "lCIOR_weighted", 
  "uCIOR_weighted", 
  "beta_unweighted",
  "lowerCI_unweighted", 
  "upperCI_unweighted", 
  "p_unweighted",
  "OR_unweighted", 
  "lCIOR_unweighted",
  "uCIOR_unweighted")

