# Health-and-Wellbeing-MR
The scripts in this repository are for conducting the analysis in Wootton, Lawn, Millard, Davies, Taylor, Munafò, Timpson, Davis, Davey Smith & Haworth, Evaluating the causal effects between subjective wellbeing and cardiometabolic health using Mendelian randomisation

The scripts relate to two separate analyses:
1. A bi-directional two-sample Mendelian randomisation between physical health and subjective wellbeing
2. A one-sample Mendelian randomisation to replicate the effects of BMI on subjective wellbeing. 

Run the scripts in the following order to replicate analysis…

**Two-sample MR**
- Two-sample MR Base script.R
- Two-sample MR Base script - Revisions.R (the addition of waist-to-hip ratio, body fat, waist circumference and blood pressure)

**One-sample MR**
- One-sample MR Script - 1 Create phenotype file.R
- Plink Guide - 2a Creating polygenic risk scores.txt
	* Script1_extractBMISNPs.sh
	* Script2_create_recoded_genotypes_BMI.sh
	* Script3_creating_profiles.sh
	* BMI_betas.txt
	* BMISNPs.txt
- One-sample MR Script - 2b Format PRS.R
- One-sample MR Script - 3 Run one-sample.R
