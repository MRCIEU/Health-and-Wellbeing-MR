#This scipt uses scripts from the MRCIEU GitHub (https://github.com/MRCIEU/TwoSampleMR) - using version TwoSampleMR_0.4.11
#Adapted to explore the association between wellbeing and physical health using two-sample MR
#This script contains the additional analysis conducted for the BMJ revisions
#Robyn May 2018
######################################################################################
#Contents
#1. Load Packages

##ADIPOSITY
#2. Read subjective wellbeing exposure data (5 x 10-8)
#3. Read in adiposity outcome data
#4. Harmonise data
#5. Run 2 sample MR
#6. Read subjective wellbeing exposure data (5 x 10-5)
#7. Read in adiposity outcome data
#8. Harmonise data
#9. Run 2 sample MR
#10. MR-PRESSO
#11. SIMEX Correction
#12. Adiposity as the exposure
#13. Regression dilution 
#14. SIMEX Correction

###BLOOD PRESSURE
#15. Read subjective wellbeing exposure data (5 x 10-8)
#16. Read in BP outcome data
#17. Harmonise data
#18. Run 2 sample MR
#19. Read subjective wellbeing exposure data (5 x 10-5)
#20. Read in BP outcome data
#21. Harmonise data
#22. Run 2 sample MR
#23. Simex correction 
#24. Blood Pressure as the exposure
#25. Regression dilution 
#26. SIMEX Correction

######################################################################################
#1. Load Packages
######################################################################################
#First make sure you have the biomaRt and devtools packages installed:
rm(list=ls(all=TRUE))

#Load packages
#install.packages("devtools")
library(devtools)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("knitr")
library(knitr)
#install_github("MRCIEU/MRInstruments")
#The current results were achieved with using version TwoSampleMR_0.4.11
library(MRInstruments)


# #Preparing blood pressure for later
# #DBP
# load("/pathname/Blood Pressure/DBP_summary_meta_ALL_June2015_FILTERED_NEFF60MAF1pc_rsid.rdata")
# str(temp)

# #currently doesn't have any column names so read it in to add them - NOTE: Louise says that the effects are in ref to the alt allele
# colnames(temp) <- c("drop", "other_allele", "effect_allele", "eaf", "beta", "se", "pval", "N", "CHR", "POS", "drop2", "drop3", "drop4", "SNP", "drop6")
# str(temp)

# #add a phenotype column
# temp$Phenotype<-"DBP"
# dbp<-temp[,c(2:10, 14)]
# str(dbp) 

# #Save dataset
# setwd("/pathname/")
# write.table(dbp, "DBP_GWAS_MRBASE.txt", quote=FALSE, row.names=FALSE)

# #Save the genomw-wide sig SNPs
# sig_dbp<-subset(dbp, dbp$pval<0.00000005)
# str(sig_dbp)

# #Save dataset
# setwd("/pathname")
# write.table(sig_dbp, "DBP_GWAS_MRBASE_sig.txt", quote=FALSE, row.names=FALSE)

# #SBP
# load("pathname/Blood Pressure/SBP_summary_meta_ALL_June2015_FILTERED_NEFF60MAF1pc_rsid.rdata")
# str(temp)

# #currently doesn't have any column names so read it in to add them - NOTE: Louise says that the effects are in ref to the alt allele
# colnames(temp) <- c("drop", "other_allele", "effect_allele", "eaf", "beta", "se", "pval", "N", "CHR", "POS", "drop2", "drop3", "drop4", "SNP", "drop6")
# str(temp)

# #add a phenotype column
# temp$Phenotype<-"SBP"
# sbp<-temp[,c(2:10, 14)]
# str(sbp) 

# #Save dataset
# setwd("pathname")
# write.table(sbp, "SBP_GWAS_MRBASE.txt", quote=FALSE, row.names=FALSE)

# #Save the genome-wide sig SNPs
# sig_sbp<-subset(sbp, sbp$pval<0.00000005)
# str(sig_sbp)

# #Save dataset
# setwd("pathname")
# write.table(sig_sbp, "SBP_GWAS_MRBASE_sig.txt", quote=FALSE, row.names=FALSE)

######################################################################################################################################
#2. Read subjective wellbeing exposure data (5 x 10-8)
######################################################################################################################################
# # Read in the exposure data with 23andMe included 
# setwd("pathname")
# wb<-read.table("SSGAC-SignificantHits.txt", header=TRUE)
# str(wb)

# # #rename the columns correctly
# library("reshape")
# wb1<-rename(wb, c(SNPID="SNP", Effect_allele ="effect_allele", Other_allele="other_allele", Beta="beta", SE="se", EAF="eaf", P.value="pval"))
# str(wb1)
# wb1$Phenotype<-"SWB"

# # #Save dataset
# setwd("pathname")
# write.table(wb1, "SSGAC-SignificantHits_MRBASE.txt", quote=F, row.names=F)

#Read in using MR Base
wb_exp_dat <- read_exposure_data(filename = "/Users/Robyn/Google Drive/PhD/Research/SSGAC/Data/SSGAC-SignificantHits_MRBASE.txt")
str(wb_exp_dat) #N=3

######################################################################################################################################
#3. Read in ADIPOSITY outcome data
######################################################################################################################################
#Search for ADIPOSITY in the available outcomes
ao <- available_outcomes()
head(ao)

whr<-subset(ao, ao$consortium=="GIANT") #id=79 - wthR id=61 - waist circum
fat <- subset(ao, grepl("fat", trait)) #id=999

outcome_dat <- extract_outcome_data(wb_exp_dat$SNP, c(79, 61, 999))

######################################################################################################################################
#4. Harmonise data
######################################################################################################################################
# harminising data (NOTE: action 2 means check palindromes using EAF)
dat1 <- harmonise_data( 
	exposure_dat = wb_exp_dat,
	outcome_dat = outcome_dat,
	action = 2
)

######################################################################################################################################
#5. Run 2 sample MR
######################################################################################################################################
# heterogeneity measures
mr_het1 <- mr_heterogeneity(dat1)

# primary methods
res1 <- mr(dat1, method_list=c("mr_ivw", "mr_weighted_median")) #Cannot conduct Egger or mode with only 3 SNPs

# single SNP analyses
res1_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_weighted_median"))

# leave one out analyses
res1_loo <- mr_leaveoneout(dat1)

# Main results
res1
mr_het1
mr_ruck
res1_single
res1_loo
  
######################################################################################################################################
#6. Read subjective wellbeing exposure data (5 x 10-5)
######################################################################################################################################

ao <- available_outcomes()

#Always check that the column numbers are the same because they change as new GWAS are added
#find wellbeing
str(ao)
ex<-subset(ao, ao$consortium=="SSGAC") # id=1009, c=13
head(ex)

#extract wellbeing snps
exposure_dat <- extract_instruments(outcomes=1009, p1=5e-5)
exposure_dat_5 <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001) #84 SNPs

######################################################################################################################################
#7. Read in adiposity outcome data
######################################################################################################################################
#Search for adiposity in the available outcomes
ao <- available_outcomes()
head(ao)

whr<-subset(ao, ao$consortium=="GIANT") #id=79
fat <- subset(ao, grepl("fat", trait)) #id=999

outcome_dat <- extract_outcome_data(exposure_dat_5$SNP, c(79, 61, 999))

######################################################################################################################################
#8. Harmonise data
######################################################################################################################################
# harminising data (NOTE: action 2 means check palindromes using EAF)
dat1 <- harmonise_data( 
	exposure_dat = exposure_dat_5 ,
	outcome_dat = outcome_dat,
	action = 2
)

######################################################################################################################################
#9. Run 2 sample MR
######################################################################################################################################
# heterogeneity measures
mr_het1 <- mr_heterogeneity(dat1)

# primary methods
res1 <- mr(dat1, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

# Egger intercept
mr_egger_int1 <- mr_pleiotropy_test(dat1)

# single SNP analyses
res1_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

# leave one out analyses
res1_loo <- mr_leaveoneout(dat1)

# Main results
res1
mr_het1
mr_egger_int1
res1_single
res1_loo

######################################################################################################################################
#10. MR-PRESSO
######################################################################################################################################
#MR-PRESSO
press<-run_mr_presso(dat1, NbDistribution = 5000, SignifThreshold = 0.05)

#install.packages("rlist")
library("rlist")
# #Save results
setwd("pathname")
list.save(press, 'mrpresso_wb-health.rdata')

#################################################################################
#11. Simex correction 
#################################################################################
#Run the simex correction for each phenotype
#install.packages("simex")
#load package
library(simex) 

#get wellbeing instrument
exposure_dat <- extract_instruments(outcome=1009, p1 = 5e-05)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#create empty dataframe to store output
simexegger<-c()

cols<-c(79, 61, 999)

#run simex in a loop
for(i in 1:length(cols)){
	Vars<-cols[i]
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(Vars), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
str(dat)

#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)

#extract results in beta format
beta2<-mod2$coefficients$jackknife[2,1]
se2<-mod2$coefficients$jackknife[2,2]
p2<-mod2$coefficients$jackknife[2,4]

#convert to odds ratios for categorical outcomes
results2<-cbind("unweighted", beta2, se2, p2)
colnames(results2) <- c("exposure", "b", "se", "pval") #following the MRBase naming convention
results2<-data.frame(results2)
results2$b<-as.numeric(as.character(results2$b))
results2$se<-as.numeric(as.character(results2$se))
results2$pval<-as.numeric(as.character(results2$pval))

or<-generate_odds_ratios(results2)

#extract confidence intervals and odds ratios
or1<-or[1,7]
lcior1<-or[1,8]
ucior1<-or[1,9]
lci2<-or[1,5]
uci2<-or[1,6]

#Save results
output<-cbind(Vars, beta2, lci2, uci2, p2, or1, lcior1, ucior1)
simexegger<-rbind(simexegger, output)
colnames(simexegger) <- c("Exposure", "beta_unweighted", "lowerCI_unweighted", "upperCI_unweighted", "p_unweighted", "OR_unweighted", "lCIOR_unweighted", "uCIOR_unweighted")
setwd("pathname")
write.csv(simexegger, file="mreggersimex_wb-health.csv", row.names = FALSE)

}


######################################################################################
#12. Adiposity as the exposure
######################################################################################
ao <- available_outcomes()

whr<-subset(ao, ao$consortium=="GIANT") #id=79, c=878 #id =61 c=679
fat <- subset(ao, grepl("fat", trait)) #id=999, c=1094

#extract exposure snps
exposure_dat <- extract_instruments(outcomes=c(79,61,999))
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

# Egger intercept
mr_egger_int1 <- mr_pleiotropy_test(dat)

#heterogeneity
mr_het1 <- mr_heterogeneity(dat)

#run MR
mr_results <- mr(dat)
mr_results

#Save results
setwd("pathname")
write.csv(mr_results, "MRBase_health_exposure.csv", row.names=F, quote=F)

#MR-PRESSO
press<-run_mr_presso(dat, NbDistribution = 5000, SignifThreshold = 0.05)

######################################################################################
#13. Regression dilution I2 GX
######################################################################################
# I-squared function
Isq <- function(y,s){
k          = length(y)
w          = 1/s^2; sum.w  = sum(w)
mu.hat     = sum(y*w)/sum.w  
Q          = sum(w*(y-mu.hat)^2)
Isq        = (Q - (k-1))/Q
Isq        = max(0,Isq)
return(Isq)
}

#calculate isq one at a atime for each exposure
##Waist-to-Hip ratio
#extract exposure snps
exposure_dat <- extract_instruments(outcome=79)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted = 0.575
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

##Waist circum
#extract exposure snps
exposure_dat <- extract_instruments(outcome=61)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted = 0.882
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted


##Body fat
#extract exposure snps
exposure_dat <- extract_instruments(outcome=999)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted = 0.565
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted


#################################################################################
#14. Simex correction 
#################################################################################
#Run the simex correction for each phenotype
#install.packages("simex")
#load package
library(simex) 

##Waist-to-Hip ratio
#extract exposure snps
exposure_dat <- extract_instruments(outcome=79)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation  
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)

##Waist Circumference
#extract exposure snps
exposure_dat <- extract_instruments(outcome=61)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation  
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)


##Body fat
#extract exposure snps
exposure_dat <- extract_instruments(outcome=999)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation  
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)


######################################################################################################################################
#15. Read subjective wellbeing exposure data (5 x 10-8)
######################################################################################################################################
rm(list=ls(all=TRUE))

# # Read in the exposure data with 23andMe included 
# setwd("/Users/Robyn/Google Drive/PhD/Research/SSGAC/Data")
# wb<-read.table("SSGAC-SignificantHits.txt", header=TRUE)
# str(wb)

# # #rename the columns correctly
# library("reshape")
# wb1<-rename(wb, c(SNPID="SNP", Effect_allele ="effect_allele", Other_allele="other_allele", Beta="beta", SE="se", EAF="eaf", P.value="pval"))
# str(wb1)
# wb1$Phenotype<-"SWB"

# # #Save dataset
# setwd("/Users/Robyn/Google Drive/PhD/Research/SSGAC/Data")
# write.table(wb1, "SSGAC-SignificantHits_MRBASE.txt", quote=F, row.names=F)

#Read in using MR Base
wb_exp_dat <- read_exposure_data(filename = "pathname/SSGAC-SignificantHits_MRBASE.txt")
str(wb_exp_dat) #N=3

######################################################################################################################################
#16. Read in BP outcome data
######################################################################################################################################
#Read in DBP data
dbp_out_dat <- read_outcome_data("pathname/Blood Pressure/DBP_GWAS_MRBASE.txt", snps=wb_exp_dat$SNP)
str(dbp_out_dat)

#Read in sbp data
sbp_out_dat <- read_outcome_data("pathname/Blood Pressure/SBP_GWAS_MRBASE.txt", snps=wb_exp_dat$SNP)
str(sbp_out_dat)

dbp_out_dat$outcome<-"DBP"
sbp_out_dat$outcome<-"SBP"

#join the two together
out_dat<-rbind(dbp_out_dat, sbp_out_dat)
str(out_dat)

######################################################################################################################################
#17. Harmonise data
######################################################################################################################################
# harminising data (NOTE: action 2 means check palindromes using EAF)
dat1 <- harmonise_data( 
	exposure_dat = wb_exp_dat,
	outcome_dat = out_dat,
	action = 2
)


######################################################################################################################################
#18. Run 2 sample MR
######################################################################################################################################
# heterogeneity measures
mr_het1 <- mr_heterogeneity(dat1)

# primary methods
res1 <- mr(dat1, method_list=c("mr_ivw", "mr_weighted_median")) #Cannot conduct Egger with only 3 SNPs

# single SNP analyses
res1_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_weighted_median", "mr_weighted_mode"))

# leave one out analyses
res1_loo <- mr_leaveoneout(dat1)

# Main results
res1
mr_het1
res1_single
res1_loo

  
######################################################################################################################################
#19. Read subjective wellbeing exposure data (5 x 10-5)
######################################################################################################################################
ao <- available_outcomes()

#Always check that the column numbers are the same because they change as new GWAS are added
#find wellbeing
str(ao)
ex<-subset(ao, ao$consortium=="SSGAC") # id=1009, c=13
head(ex)

#extract wellbeing snps
exposure_dat <- extract_instruments(outcome=1009, p1=5e-5)
exposure_dat_5 <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001) #84 SNPs

######################################################################################################################################
#20. Read in BP outcome data
######################################################################################################################################
#Read in DBP data
dbp_out_dat <- read_outcome_data("pathname/Blood Pressure/DBP_GWAS_MRBASE.txt", snps= exposure_dat_5 $SNP)
str(dbp_out_dat)

#Read in sbp data
sbp_out_dat <- read_outcome_data("pathname/Blood Pressure/SBP_GWAS_MRBASE.txt", snps= exposure_dat_5 $SNP)
str(sbp_out_dat)

dbp_out_dat$outcome<-"DBP"
sbp_out_dat$outcome<-"SBP"

#join the two together
out_dat<-rbind(dbp_out_dat, sbp_out_dat)
str(out_dat)

######################################################################################################################################
#21. Harmonise data
######################################################################################################################################
# harminising data (NOTE: action 2 means check palindromes using EAF)
dat1 <- harmonise_data( 
	exposure_dat = exposure_dat_5 ,
	outcome_dat = out_dat,
	action = 2
)

######################################################################################################################################
#22. Run 2 sample MR
######################################################################################################################################
# heterogeneity measures
mr_het1 <- mr_heterogeneity(dat1)

# primary methods
res1 <- mr(dat1, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

# Egger intercept
mr_egger_int1 <- mr_pleiotropy_test(dat1)

# single SNP analyses
res1_single <- mr_singlesnp(dat1, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))

# leave one out analyses
res1_loo <- mr_leaveoneout(dat1)

# Main results
res1
mr_het1
mr_egger_int1
res1_single
res1_loo

#################################################################################
#23. Simex correction 
#################################################################################
#Run the simex correction for each phenotype
#install.packages("simex")
#load package
library(simex) 

##DBP

dat <- harmonise_data( 
	exposure_dat = exposure_dat_5 ,
	outcome_dat = dbp_out_dat,
	action = 2
)

#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation  
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)


##sbp
dat1 <- harmonise_data( 
	exposure_dat = exposure_dat_5 ,
	outcome_dat = sbp_out_dat,
	action = 2
)

#Rename required columns
dat1$BetaXG<-dat1$beta.exposure
dat1$seBetaXG<-dat1$se.exposure
dat1$BetaYG<-dat1$beta.outcome
dat1$seBetaYG<-dat1$se.outcome
BetaXG <- dat1$BetaXG
BetaYG <- dat1$BetaYG
seBetaXG <- dat1$seBetaXG
seBetaYG <- dat1$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation  
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)


######################################################################################
#24. BP as the exposure
######################################################################################
##DBP
#Read in DBP data
dbp_exp_dat <- read_exposure_data("pathname/Blood Pressure/DBP_GWAS_MRBASE_sig.txt") #excludes indels
str(dbp_exp_dat)

dbp_exp_dat$exposure<-"DBP"
dbp_exp_dat <- clump_data(dbp_exp_dat, clump_kb = 10000, clump_r2 = 0.001) #63 SNPs
str(dbp_exp_dat)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(dbp_exp_dat $SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(dbp_exp_dat, outcome_dat, action = 2)

#run MR
mr_results <- mr(dat)
mr_results

#Egger intercept
mr_egger_int1 <- mr_pleiotropy_test(dat)

# heterogeneity measures
mr_het1 <- mr_heterogeneity(dat)

#MR-PRESSO
press<-run_mr_presso(dat, NbDistribution = 5000, SignifThreshold = 0.05)

##SBP
#Read in DBP data
sbp_exp_dat <- read_exposure_data("pathname/Blood Pressure/SBP_GWAS_MRBASE_sig.txt") #excludes indels
str(sbp_exp_dat)

sbp_exp_dat $exposure<-"SBP"
sbp_exp_dat <- clump_data(sbp_exp_dat, clump_kb = 10000, clump_r2 = 0.001) #64 SNPs
str(sbp_exp_dat)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(sbp_exp_dat $SNP, c(1009), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat1 <- harmonise_data(sbp_exp_dat, outcome_dat, action = 2)

#run MR
mr_results <- mr(dat1)
mr_results

#Egger intercept
mr_egger_int1 <- mr_pleiotropy_test(dat1)

# heterogeneity measures
mr_het1 <- mr_heterogeneity(dat1)

#MR-PRESSO
press<-run_mr_presso(dat1, NbDistribution = 5000, SignifThreshold = 0.05)


######################################################################################
#25. Regression dilution I2 GX
######################################################################################
# I-squared function
Isq <- function(y,s){
k          = length(y)
w          = 1/s^2; sum.w  = sum(w)
mu.hat     = sum(y*w)/sum.w  
Q          = sum(w*(y-mu.hat)^2)
Isq        = (Q - (k-1))/Q
Isq        = max(0,Isq)
return(Isq)
}

#calculate isq one at a atime for each exposure
##DBP
#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted 
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

##SBP
#Rename required columns
dat1 $BetaXG<-dat1 $beta.exposure
dat1$seBetaXG<-dat1 $se.exposure
BetaXG   = dat1 $BetaXG
seBetaXG = dat1 $seBetaXG 
seBetaYG<-dat1 $se.outcome

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted 
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted


#################################################################################
#26. Simex correction 
#################################################################################
#Run the simex correction for each phenotype
#install.packages("simex")
#load package
library(simex) 

##DBP
#Rename required columns
dat$BetaXG<-dat$beta.exposure
dat$seBetaXG<-dat$se.exposure
dat$BetaYG<-dat$beta.outcome
dat$seBetaYG<-dat$se.outcome
BetaXG <- dat$BetaXG
BetaYG <- dat$BetaYG
seBetaXG <- dat$seBetaXG
seBetaYG <- dat$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation  
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)


##sbp
#Rename required columns
dat1$BetaXG<-dat1$beta.exposure
dat1$seBetaXG<-dat1$se.exposure
dat1$BetaYG<-dat1$beta.outcome
dat1$seBetaYG<-dat1$se.outcome
BetaXG <- dat1$BetaXG
BetaYG <- dat1$BetaYG
seBetaXG <- dat1$seBetaXG
seBetaYG <- dat1$seBetaYG


BYG <- BetaYG*sign(BetaXG)# Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)         

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation  
mod.sim2 <- simex(Fit2,B=1000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod2<-summary(mod.sim2)
