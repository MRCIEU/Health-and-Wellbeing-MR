#This scipt uses scripts from the MRCIEU GitHub (https://github.com/MRCIEU/TwoSampleMR) - Using version TwoSampleMR_0.4.11
#Adapted to explore the association between wellbeing and physical health using two-sample MR
#Robyn June 2017
######################################################################################
#Contents
#1. Load Packages
#2. Subjective wellbeing on physical health (p<5x10-8)
#3. Subjective wellbeing on physical health (p<5x10-5)
#4. MR-PRESSO
#5. Physical health on subjective wellbeing
#6. MR-PRESSO
#7. Plot the effects of BMI on subjective wellbeing
#8. Regression dilution I2 - wellbeing
#9. Regression dilution I2 - health
#10. SIMEX Correction

######################################################################################
#1. Load Packages
######################################################################################
#First make sure you have the biomaRt and devtools packages installed:
#install.packages("devtools")
#install.packages("biomaRt")

#Then to install:
library(devtools)
library(biomaRt)
#install_github("MRCIEU/TwoSampleMR") 
#To update the package just run the install_github("MRCIEU/TwoSampleMR") command again.
#The current results were achieved with using version TwoSampleMR_0.4.11
library(TwoSampleMR)

#Note: if these do not run then they might be working on updating the website

######################################################################################
#2. Subjective wellbeing on physical health (p<5x10-8)
######################################################################################
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
wb_exp_dat <- read_exposure_data(filename = "pathname/SSGAC-SignificantHits_MRBASE.txt")
str(wb_exp_dat) #N=3

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(wb_exp_dat$SNP, c(2,7,299,300,301,798), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(wb_exp_dat, outcome_dat, action = 2)

#run MR
mr_results <- mr(dat, method_list=c("mr_ivw", "mr_weighted_median"))
mr_results

or<-generate_odds_ratios(mr_results)
or

# #Save results
# setwd("pathname")
# write.csv(mr_results, "MRBase_RScript_wellbeing_exposure_health_outcome.csv", row.names=F, quote=F)

######################################################################################
#3. Subjective wellbeing on physical health (p<5x10-5)
######################################################################################
ao <- available_outcomes()

#Always check that the column numbers are the same because they change as new GWAS are added
#find wellbeing
str(ao)
ex<-subset(ao, ao$consortium=="SSGAC") #use the SWB full set (including 23andme) for the exposure (id=1009, c=13) Note. column might change but ID remains the same
str(ex)

#extract wellbeing snps at the relaxed threshold
exposure_dat <- extract_instruments(outcomes=1009, p1 = 5e-05)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#extract outcomes
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c(2,7,299,300,301,798), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#run MR
mr_results <- mr(dat)
mr_results

#generate odds ratios 
or<-generate_odds_ratios(mr_results)
or

# #Save results
# setwd("pathname")
# write.csv(or, "MRBase_RScript_wellbeing_exposure_health_outcome5e-5.csv", row.names=F, quote=F)

# Egger intercept
mr_egger_int1 <- mr_pleiotropy_test(dat)

######################################################################################
#4. MR-PRESSO
######################################################################################
#MR-PRESSO
press<-run_mr_presso(dat, NbDistribution = 5000, SignifThreshold = 0.05)

#install.packages("rlist")
library("rlist")
# #Save results
setwd("pathname")
list.save(press, 'mrpresso_wb-health.rdata')

######################################################################################
#5. Physical health on wellbeing
######################################################################################
ao <- available_outcomes()

#Identify exposures
table(ao$consortium)
bmi<-subset(ao, ao$consortium=="GIANT") #(id = 2, c = 249) - make sure col number is correct
chol<-subset(ao, ao$consortium=="GLGC") #HDL = (id = 299, c = 344), LDL = (id=300, c=347), Total = (id=301, c=348) 
heart<-subset(ao, ao$consortium=="CARDIoGRAMplusC4D") #CAD = (id=7, c=778), MI = (id=798, c=887) 
exposure_dat <- extract_instruments(outcomes=c(2,7,299,300,301,798))
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=1009, proxies = T, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)
mr_results_health <- mr(dat)
mr_results_health

# #Save results
# setwd("pathname")
# write.csv(mr_results_health, "MRBase_RScript_health_exposure_wellbeing_outcome.csv", row.names=F, quote=F)

# Egger intercept
mr_egger_int1 <- mr_pleiotropy_test(dat)

# #Save results
# setwd("pathname")
# write.csv(mr_egger_int1, "MRBase_health_exposure_wellbeing_outcome_intercept.csv", row.names=F, quote=F)

######################################################################################
#6. MR-PRESSO
######################################################################################
#MR-PRESSO
press<-run_mr_presso(dat, NbDistribution = 5000, SignifThreshold = 0.05)

#install.packages("rlist")
library("rlist")
# #Save results
setwd("pathname")
list.save(press, 'mrpresso_wb-health.rdata')

######################################################################################
#7. Plot the effects of Physical health on wellbeing
######################################################################################
#Just get BMI->wellbeing
bmi_exposure_dat <- extract_instruments(outcomes=c(2))
bmi_exposure_dat <- clump_data(bmi_exposure_dat, clump_kb = 10000, clump_r2 = 0.001)
ao <- available_outcomes()
outcome_dat <- extract_outcome_data(snps= bmi_exposure_dat$SNP, outcomes=1009, proxies = T, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

dat <- harmonise_data(bmi_exposure_dat, outcome_dat, action = 2)
mr_results_health <- mr(dat)

p <- mr_scatter_plot(mr_results_health, dat)
p[[1]]

#Leave one out analysis
l <- mr_leaveoneout(dat)
p <- mr_leaveoneout_plot(l)
p[[1]]

#Forest plot
s <- mr_singlesnp(dat)
p <- mr_forest_plot(s)
p[[1]]

#Funnel plot
p <- mr_funnel_plot(s)
p[[1]]

######################################################################################
#8. Regression dilution I2 - wellbeing
######################################################################################
I2<-c()

# I-squared function
Isq = function(y,s){
k          = length(y)
w          = 1/s^2; sum.w  = sum(w)
mu.hat     = sum(y*w)/sum.w  
Q          = sum(w*(y-mu.hat)^2)
Isq        = (Q - (k-1))/Q
Isq        = max(0,Isq)
return(Isq)
}

#Read in expsoure data
wellbeing <- extract_instruments(outcomes=c(1009), p1 = 5e-05)
wellbeing <- clump_data(wellbeing, clump_kb = 10000, clump_r2 = 0.001)

#Rename required columns
wellbeing$BetaXG<-wellbeing$beta.exposure
wellbeing$seBetaXG<-wellbeing$se.exposure
BetaXG   = wellbeing$BetaXG
seBetaXG = wellbeing$seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)

#unweightedIsq
unIsq = Isq(BXG,seBetaXG) #unweighted

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance

#Save results
output<-cbind("wellbeing", mF, unIsq, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq Unweighted",  "Isq_weighted")
setwd("pathname")
write.csv(I2, file="regression_dilution_isq.csv", row.names = FALSE)

######################################################################################
#9. Regression dilution I2 - health
######################################################################################

##BMI
#Read in expsoure data
bmi <- extract_instruments(outcomes=c(2))
bmi <- clump_data(bmi)

#Rename required columns
bmi $BetaXG<-bmi $beta.exposure
bmi $seBetaXG<-bmi $se.exposure
BetaXG   = bmi $BetaXG
seBetaXG = bmi $seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
unIsq = Isq(BXG,seBetaXG)

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance

#Save results
output<-cbind("BMI", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq")
setwd("pathname")
write.csv(I2, file="regression_dilution_isq.csv", row.names = FALSE)


##CHD
#Read in expsoure data
chd <- extract_instruments(outcomes=c(7))
chd <- clump_data(chd)

#Rename required columns
chd $BetaXG<-chd $beta.exposure
chd $seBetaXG<-chd $se.exposure
BetaXG   = chd $BetaXG
seBetaXG = chd $seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
unIsq = Isq(BXG,seBetaXG)

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance

#Save results
output<-cbind("CHD", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq")
setwd("pathname")
write.csv(I2, file="regression_dilution_isq.csv", row.names = FALSE)


##MI
#Read in expsoure data
mi <- extract_instruments(outcomes=c(798))
mi <- clump_data(mi)

#Rename required columns
mi$BetaXG<-mi$beta.exposure
mi$seBetaXG<-mi$se.exposure
BetaXG   = mi$BetaXG
seBetaXG = mi$seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
unIsq = Isq(BXG,seBetaXG)

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance

#Save results
output<-cbind("MI", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq")
setwd("pathname")
write.csv(I2, file="regression_dilution_isq.csv", row.names = FALSE)


##LDL cholesterol
#Read in expsoure data
ldl <- extract_instruments(outcomes=c(300))
ldl <- clump_data(ldl)

#Rename required columns
ldl $BetaXG<-ldl $beta.exposure
ldl $seBetaXG<-ldl $se.exposure
BetaXG   = ldl $BetaXG
seBetaXG = ldl $seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
unIsq = Isq(BXG,seBetaXG)

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance

#Save results
output<-cbind("LDL", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq")
setwd("pathname")
write.csv(I2, file="regression_dilution_isq.csv", row.names = FALSE)

##HDL cholesterol
#Read in expsoure data
hdl <- extract_instruments(outcomes=c(299))
hdl <- clump_data(hdl)

#Rename required columns
hdl$BetaXG<-hdl$beta.exposure
hdl$seBetaXG<-hdl$se.exposure
BetaXG   = hdl$BetaXG
seBetaXG = hdl$seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
unIsq = Isq(BXG,seBetaXG)

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance

#Save results
output<-cbind("HDL", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq")
setwd("pathname")
write.csv(I2, file="regression_dilution_isq.csv", row.names = FALSE)


##Total cholesterol
#Read in expsoure data
chol <- extract_instruments(outcomes=c(301))
chol <- clump_data(chol)

#Rename required columns
chol $BetaXG<-chol $beta.exposure
chol $seBetaXG<-chol $se.exposure
BetaXG   = chol $BetaXG
seBetaXG = chol $seBetaXG 

BXG             = abs(BetaXG)         # gene--exposure estimates are positive  

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

F   = BXG^2/seBetaXG^2
mF  = mean(F)
unIsq = Isq(BXG,seBetaXG)

mF  # Want mF to be high for good IVW performance
unIsq # Want Isq to be close to 1 for good MR-Egger performance

#Save results
output<-cbind("Total Chol", mF, unIsq)
I2<-rbind(I2, output)
colnames(I2) <- c("Exposure", "mF", "Isq")
setwd("pathname")
write.csv(I2, file="regression_dilution_isq.csv", row.names = FALSE)


#################################################################################
#10. SIMEX Correction
#################################################################################
#Run the simex correction for each phenotype
#install.packages("simex")
#load package
library(simex) 

#get wellbeing instrument
exposure_dat <- extract_instruments(outcomes=c(1009), p1 = 5e-05)
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001)

#create empty dataframe to store output
simexegger<-c()

cols<-c(2,7,299,300,301,798)

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











