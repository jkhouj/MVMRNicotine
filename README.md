# MVMRNicotine
This repository contains code relating to a project in which we use multivariable Mendelian randomisation (MVMR) to explore the health effects of nicotine use

Read me file for scripts relating to the project: MVMR of nicotine and non-nicotine constituents.

1. Scripts included

Script "1. MVMR Cancer.R" conducts univariable and multivariable MR exploring the effects of nicotine a
and non-nicotine constituents of tobacco smoke (measured by nicotine metabolite ratio [NMR] and 
cigarettes per day [CPD] on lung cancer.

Script "2. MVMR Cancer_cotinine.r" conducts univariable and multivariable MR exploring the effects of
nicotine and non-nicotine constituents of tobacco smoke (measured by cotinine [cot] and cigarettes per 
day [CPD] on lung cancer.

Script "3. MVMR Cancer_3HCplusCOT.r" conducts univariable and multivariable MR exploring the effects of
nicotine and non-nicotine constituents of tobacco smoke (measured by 3HC plus cotinine [Cplus3HC] and 
cigarettes per day [CPD] on lung cancer.

Script "1. MVMR_NMR_healt.R" conducts univariable and multivariable MR exploring the effects of nicotine 
and non-nicotine constituents of tobacco smoke (measured by nicotine metabolite ratio [NMR] and cigarettes
per day [CPD]) on health outcomes (chronic obstructive pulmonary disease [COPD], coronary heart disease 
[CHD], forced expiratory volume [FEV], forced vital capacity [FVC], heart rate [HR], and body mass index 
[BMI])

Script "2. MVMR_cot_Health.R" conducts univariable and multivariable MR exploring the effects 
of nicotine and non-nicotine constituents of tobacco smoke (measured by cotinine [cot] and cigarettes 
per day [CPD] on  health outcomes (chronic obstructive pulmonary disease [COPD], coronary heart disease 
[CHD], forced expiratory volume [FEV], forced vital capacity [FVC], heart rate [HR], and body mass
index [BMI])

Script "3. MVMR_Cotplus3HC_health.R" conducts univariable and multivariable MR exploring the effects of
nicotine and non-nicotine constituents of tobacco smoke (measured by 3HC plus cotinine [Cplus3HC] and 
cigarettes per day [CPD]) on health outcomes (chronic obstructive pulmonary disease [COPD], coronary 
heart disease [CHD], forced expiratory volume [FEV], forced vital capacity [FVC], heart rate [HR],
and body mass index [BMI])

SIMEX_Cot.R conducts SIMEX for script "2. MVMR Cancer_cotinine.r"

SIMEX_Cplus3HC.R conducts SIMEX for script "3. MVMR Cancer_3HCpluCOT.r" 

SIMEX_Cot_Rep.R conducts SIMEX for script "2. MVMR_cot_Health.R"

SIMEX_Cotplus3HC.R conducts SIMEX for script "3. MVMR_Cotplus3HC_health.R" 

Proxy searching files are not included.

The scripts were run using R version 4.1.1. 

2. Data availability
GSCAN CPD data are available at https://doi.org/10.13020/3b1n-ff32 
(note the datasets with UK Biobank removed are used for analysis with outcome data derived in UK Biobank)
NMR and cotinine plus 3HC data are available on request from the authors https://doi.org/10.1038/s41380-020-0702-z
Cotinine data are available at https://doi.org/10.5523/bris.182rhz19hg3lz1172a7yfcap9v
SNP lists are available in this Github folder.

3. Main packages used, version numbers, dependencies and suggestions

MVMR package
Version: 0.3
Dependencies:R>=3.5, boot
Suggested: captioner, kableExtra, knitr, MendelianRandomization,\nrmarkdown, stringr

TwoSampleMR package
Version: 0.5.6 
Dependencies:R>=3.6, ieugwasr (>= 0.1.5), ggplot2, gridExtra, cowplot, plyr,\nreshape2, stringr, 
	knitr, markdown, gtable, rmarkdown,\nMendelianRandomization, dplyr, mr.raps, psych, 
	magrittr, car,\nrandomForest, meta, data.table, MRPRESSO, MRInstruments,\nRadialMR, 
	MRMix, glmnet, lattice, pbapply, MASS, Cairo
Suggested: covr, testthat

Mendelian Randomisation package
Version: 0.5.1
Dependencies:R>= 3.0.1, knitr, rmarkdown, plotly (>= 3.6.0), ggplot2 (>= 1.0.1),\nrobustbase (>= 0.92-6), 
	Matrix (>= 1.2), iterpc (>= 0.3),\nquantreg (>= 5.01), rjson, glmnet
Suggested:NA

