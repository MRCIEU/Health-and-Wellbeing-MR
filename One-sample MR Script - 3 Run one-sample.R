#Prep for running One-sample MR using UK Biobank data
#x->BMI
#y->Wellbeing
######################################################################################
#Contents
#1. Load BMI polygenic scores
#2. Load phenotype data
#3. Link the genetic and phenotypic data
#4. Exclude all individuals without any genetic data (because then they haven't met the new QC exclusions)
#5. Test Assosiation of genetic score with BMI
#6. Descriptive stats - 500k
#7. Creating Histograms
#8. Test Assosiation of BMI with baseline confounders
#9. Test Assosiation of genetic score with baseline confounders
#10. Phenotypic association between BMI and wellbeing
#11. Merge in principal components
#12. One sample MR - using ivreg
#13. Sensitivity analysis with the 350k (and no relateds) 
#14. Sensitivity analysis - splitting down age
#15. Sensitivity analysis - split by sex
######################################################################################
#1. Load BMI genotypes
######################################################################################
rm(list=ls())

#Readin genotype file
setwd("pathname")
scores<-read.csv("UKBiobank_BMI_PRS.csv", header=T) #Check N=337,115

######################################################################################
#2. Load phenotype data
######################################################################################
#Readin genotype file
setwd("pathname")
phen<-read.csv("UKBiobank_BMIandWellbeing_Phenotypes.csv", header=T) 
str(phen)#Check N=502,647

######################################################################################
#3. Link the genetic and phenotypic data
######################################################################################
#read in the linkage file
setwd("pathname")
link<-read.csv("data.7341.csv", header=T)
str(link)

#merge the link file with the phenotype file
phen.link<-merge(phen, link, by.x="eid", by.y="app16729", all.x=T)
str(phen.link)

#Merge scores and phenotype data using the genetic identifier
all.equal(phen.link$app8786, scores$IID) #only 337,115 of the 500,000 match - this is what we expected because we have removed the rest with our exclusion criteria

data<-merge(phen.link, scores, by.x="app8786", by.y="IID", all.x=T)
str(data)
summary(data)

######################################################################################
#4. Exclude all individuals without any genetic data (because then they haven't met the new QC exclusions)
######################################################################################
str(data)
ex<-subset(data, (complete.cases(data$SCORE500)|complete.cases(data$SCORE150)|complete.cases(data$SCORE350) ))
str(ex)
#N=337,112

#And any individuals who have withdrawn consent
setwd("pathname")
exclude<-read.csv("w914_20170726.csv", header=F)
str(exclude)
library(reshape)
exclude<-rename(exclude, c(V1="IID"))
exclude$exclude<-1
excluded<-merge(ex, exclude, by.x="app8786", by.y="IID", all.x=T)
str(excluded)
summary(excluded)
#None of the withdrawn individuals are within this dataset anyway

library(psych)
describe(excluded)
######################################################################################
#5. Test Assosiation of genetic score with BMI
######################################################################################
#install.packages("foreign")
#install.packages("plyr")
#install.packages("AER")
# load packages into this R session
library(foreign)
library(plyr)
library(AER)

#rename dataset
data<-ex
str(data)

# Association of IV (genes) with exposure (BMI)
##### 150k
model1<-lm(data$BMI1~data$SCORE150) 
summary(model1)
plot(model1) 

##### 500k
model1<-lm(data$BMI1~data$SCORE500) 
summary(model1)

##### 350k
model1<-lm(data$BMI1~data$SCORE350) 
summary(model1)


######################################################################################
#6. Descriptive stats - 500k (Time point 1)
######################################################################################
#BMI
mean(data$BMI1, na.rm=T)
sd(data$BMI1, na.rm=T)
#27.39 (4.75)

#happ
mean(data$happ1, na.rm=T)
sd(data$happ1, na.rm=T)
#4.45 (0.70)

#work
mean(data$work1, na.rm=T)
sd(data$work1, na.rm=T)
#4.40 (0.87)

#health
mean(data$health1, na.rm=T)
sd(data$health1, na.rm=T)
#4.25 (0.87)

#finances
mean(data$finance1, na.rm=T)
sd(data$finance1, na.rm=T)
#4.31 (0.94)

#friends
mean(data$friends1, na.rm=T)
sd(data$friends1, na.rm=T)
#4.76 (0.74)

#family
mean(data$family1, na.rm=T)
sd(data$family1, na.rm=T)
#4.79 (0.90)

#Mean baseline confounders
#Sex
table(data$sex) #0=female (N=181363) 1=male (N=155749) 54% female

#Age at baseline
mean(data$Age)
sd(data$Age) 
#56.87 (8.00)years

#age range
summary(data$Age)

#SES - Townsend deprivation index
mean(data$SES, na.rm=T) 
sd(data$SES, na.rm=T) 
#-1.58 (2.93)

#get skew and kurtosis
library(psych)
wellbeing<-subset(data, select=c("happ1", "work1", "health1", "finance1", "friends1", "family1"))
describe(wellbeing)

######################################################################################
#7. Creating Histograms
######################################################################################
gen<-data
str(gen)

setwd("pathname")
pdf("histogram_BMI.pdf",width=6,height=4,paper='special') 
hist(gen$BMI1, main="", xlab="BMI", col="lightblue3")
dev.off()
pdf("histogram_happ.pdf",width=6,height=4,paper='special') 
hist(gen$happ1, main="", xlab="Subjective Happiness", col="lightblue3", breaks=6)
dev.off()
pdf("histogram_work.pdf",width=6,height=4,paper='special') 
hist(gen$work1, main="", xlab="Satisfaction with Work", col="lightblue3", breaks=6)
dev.off()
pdf("histogram_health.pdf",width=6,height=4,paper='special') 
hist(gen$health1, main="", xlab="Satisfaction with Health", col="lightblue3", breaks=6)
dev.off()
pdf("histogram_finance.pdf",width=6,height=4,paper='special') 
hist(gen$finance1, main="", xlab="Satisfaction with Finances", col="lightblue3", breaks=6)
dev.off()
pdf("histogram_friends.pdf",width=6,height=4,paper='special') 
hist(gen$friends1, main="", xlab="Satisfaction with Friends", col="lightblue3", breaks=6)
dev.off()
pdf("histogram_family.pdf",width=6,height=4,paper='special') 
hist(gen$family1, main="", xlab="Satisfaction with Family", col="lightblue3", breaks=6)
dev.off()

######################################################################################
#8. Test Assosiation of BMI with baseline confounders
######################################################################################
#I have got sex, age at recruitment and ses as baseline confounders. 
#I need to check that the genetic scores for BMI do not predict any of these. 

#Sex
#Conduct independent t-test to test th association between sex and BMI 
t.test(gen$BMI1~gen$sex)

#Age
#Conduct linear regression to test th association between age and BMI 
model2<-lm(gen$BMI1~gen$Age)
summary(model2)

#SES
#Conduct linear regression to test th association between SES and BMI 
model3<-lm(gen$BMI1~ gen$SES)
summary(model3)


######################################################################################
#9. Test Assosiation of genetic score with baseline confounders
######################################################################################
#I have got sex, age at recruitment and ses as baseline confounders. 
#I need to check that the genetic scores for BMI do not predict any of these. 

#Sex
#Conduct independent t-test to test th association between sex and BMI score
t.test(data$SCORE500~data$sex)

#Age
#Conduct linear regression to test th association between age and BMI score
model5<-lm(data$SCORE500~data$Age)
summary(model5)


#SES
#Conduct linear regression to test th association between SES and BMI score
model6<-lm(data$SCORE500~data$SES)
summary(model6)


#Merge in education
setwd("pathname")
edu<-read.csv("Oct2017_smoking_mortality_SNP_tidied.csv", header=T)
str(edu)
str(data)

#merge in genetic IDs
#Read in linker file
setwd("pathname")
library(readstata13)
link2<-read.dta13("20160211_Mapping_phenotype_genetic.dta")
str(link2)

#merge the link file with the phenotype file
age.gen.link<-merge(edu, link2, by.x="ID", by.y="appid9142", all.x=T)
str(age.gen.link)

all.equal(age.gen.link $appid8786, data$app8786)

data.edu<-merge(data, age.gen.link, by.x="app8786", by.y="appid8786", all.x=T)
str(data.edu)

#make values of higher smoking higher
table(data.edu$eversmoking)
data.edu$eversmoking<-as.character(data.edu$eversmoking)
data.edu$eversmoking[data.edu$eversmoking=="Ever"]<-1
data.edu$eversmoking[data.edu$eversmoking=="Never"]<-0
data.edu$eversmoking<-as.numeric(data.edu$eversmoking)

table(data.edu$smoking_status)
data.edu$smoking_status<-as.character(data.edu$smoking_status)
data.edu$smoking_status[data.edu$smoking_status =="Current"]<-2
data.edu$smoking_status[data.edu$smoking_status =="Previous"]<-1
data.edu$smoking_status[data.edu$smoking_status =="Never"]<-0
data.edu$smoking_status <-as.numeric(data.edu$smoking_status)

summary(data.edu$SES)
data.edu$SES<-data.edu$SES*-1

#Save the data as .dta and create the bias plots in STATA using bias_plots.do
library("foreign")
#save stata data set
setwd("pathname")
write.dta(data.edu, "Nov2017_BMI_wellbeing.dta")


######################################################################################
#10. Phenotypic association between BMI and wellbeing
######################################################################################
str(gen) #only using people with genetic data
library(psych)
library(Hmisc) 

###Rerun as regressions
str(gen)
#Happiness
m1<-lm(gen$happ1~gen$BMI1+gen$sex+gen$Age+gen$SES)
summary(m1) 
#beta= -0.0013937  SE=0.0004423  p=0.00163 **

#Work
m1<-lm(gen$work1~gen$BMI1+gen$sex+gen$Age+gen$SES)
summary(m1) 
#beta= -0.0011064  SE=0.0006666  p=0.097

#Health
m1<-lm(gen$health1~gen$BMI1+gen$sex+gen$Age+gen$SES)
summary(m1) 
#beta= -0.0480265  SE=0.0005280 p< 2e-16 ***

#Finance
m1<-lm(gen$finance1~gen$BMI1+gen$sex+gen$Age+gen$SES)
summary(m1) 
#beta= -0.0207047 SE= 0.0005816 p< 2e-16 ***

#Friends
m1<-lm(gen$friends1~gen$BMI1+gen$sex+gen$Age+gen$SES)
summary(m1) 
#beta= 0.0016238  SE=0.0004683   p=0.000525 ***

#Family
m1<-lm(gen$family1~gen$BMI1+gen$sex+gen$Age+gen$SES)
summary(m1) 
#beta= 4.454e-05  SE=5.711e-04   p=0.938


######################################################################################
#11. Merge in principal components
######################################################################################
## Read in principal components
PC <- read.table ("/pathname/data.pca1-10.field_22009.txt", colClasses = c(rep("character", 12)))
str(PC)

#merge with the genetic id
str(gen)
gen_pc<-merge(gen, PC, by.x="app8786", by.y="V1", all.x=T)
str(gen_pc)

#rename the PC columns
library(reshape)
gen_pc<-rename(gen_pc, c(V2="ID", V3="PC1", V4="PC2", V5="PC3", V6="PC4", V7="PC5", V8="PC6", V9="PC7", V10="PC8", V11="PC9", V12="PC10"))
str(gen_pc)

#Make PCs numeric
gen_pc$PC1<-as.numeric(gen_pc$PC1)
gen_pc$PC2<-as.numeric(gen_pc$PC2)
gen_pc$PC3<-as.numeric(gen_pc$PC3)
gen_pc$PC4<-as.numeric(gen_pc$PC4)
gen_pc$PC5<-as.numeric(gen_pc$PC5)
gen_pc$PC6<-as.numeric(gen_pc$PC6)
gen_pc$PC7<-as.numeric(gen_pc$PC7)
gen_pc$PC8<-as.numeric(gen_pc$PC8)
gen_pc$PC9<-as.numeric(gen_pc$PC9)
gen_pc$PC10<-as.numeric(gen_pc$PC10)

######################################################################################
#12. One sample MR - using ivreg
######################################################################################
# perform iv regression
iv_mod <- ivreg(happ1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(work1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(health1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(finance1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)
  
# perform iv regression
  iv_mod <- ivreg(friends1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)
  

# perform iv regression
iv_mod <- ivreg(family1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
summary(iv_mod, diagnostics = TRUE)

  
######################################################################################
#13. Sensitivity analysis with the 350k (and no relateds)
######################################################################################

iv_mod <- ivreg(happ1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE350 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(work1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE350 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(health1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE350 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(finance1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE350 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)
  
# perform iv regression
  iv_mod <- ivreg(friends1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE350 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(family1 ~ BMI1 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE350 + Age +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = gen_pc)
  summary(iv_mod, diagnostics = TRUE)


######################################################################################
#14. Sensitivity analysis - splitting down age
######################################################################################
#find the median age
median(gen_pc$Age)

#split into median and younger or older than median
young<-subset(gen_pc, gen_pc$Age<=58)
table(young$Age)
str(young)
old<-subset(gen_pc, gen_pc$Age>58)
table(old$Age)
str(old)

# perform iv regression
  iv_mod <- ivreg(health1 ~ BMI1 +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = young)
  # check diagnostics
  summary(iv_mod, diagnostics = TRUE)

# perform iv regression
  iv_mod <- ivreg(health1 ~ BMI1 +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 +sex +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = old)
  # check diagnostics
  summary(iv_mod, diagnostics = TRUE)

######################################################################################
#15. Sensitivity analysis - split by sex
######################################################################################
#Mean baseline confounders
#Sex
table(data$sex) #0=female (N=181363) 1=male (N=155749) 54% female

#split into median and younger or older than median
female<-subset(gen_pc, gen_pc$sex==0)
table(female$sex)
str(female)
male<-subset(gen_pc, gen_pc$sex==1)
table(male$sex)
str(male)

# perform iv regression
  iv_mod <- ivreg(health1 ~ BMI1 +Age +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 +Age +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = female)
  # check diagnostics
  summary(iv_mod, diagnostics = TRUE) 

# perform iv regression
  iv_mod <- ivreg(health1 ~ BMI1 +Age +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10 | SCORE500 +Age +PC1 +PC2 +PC3 +PC4 +PC5 +PC6 +PC7 +PC8 +PC9+PC10, data = male)
  # check diagnostics
  summary(iv_mod, diagnostics = TRUE)

