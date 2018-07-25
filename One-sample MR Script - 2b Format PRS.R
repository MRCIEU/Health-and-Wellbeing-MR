#This script takes the PLINK output from "Plink Guide - 2a Creating polygenic risk scores.txt"
#This output has already had the standard IEU exclusions applied (DOI: 10.5523/bris.3074krb6t2frj29yh2b03x3wxj)
#In this script, I create 3 different subsets - those who were in the original 150k, those who weren't (350K), and everyone from the full release
#Robyn Wootton
######################################################################################
#To create BMI genotype profiles, see guide: Plink Guide - 2a Creating polygenic risk scores.txt
rm(list = ls())

#Genetic scores from 500k
setwd("pathname")
score500<-read.table("20180725_scores_BMI_97proxy.profile", header=T)
str(score500)
#Check N=337,115
#This data already has recommended exclusions, non-white british, highly relateds and minimal relateds removed

#Genetic scores from 150k
setwd("pathname")
score150<-read.table("20170530_scores_BMI_Louise.profile", header=T)
str(score150)
#Check N=152249
#This data has NOT yet had exlcusions applied

#Genetic scores from 350k
merge.scores<-merge(score500, score150, by="FID", all.x=T)
str(merge.scores)
score350<-subset(merge.scores, is.na(merge.scores$SCORE.y))
str(score350)
#Check N=242219
#And as this is made from the 500k sample, this should have all exclusions and relateds removed

#Create one dataset with just SCORES needed
score<-merge(score500, score150, by="FID", all.x=T)
str(score)
score<-merge(score, score350, by="FID", all.x=T)
str(score)
#Subset the variables we want
scores<-subset(score, select=c("FID", "IID.x.x", "SCORE.x.x", "SCORE.y.x", "SCORE.x.y"))
str(scores)
library(reshape)
scores<-rename(scores, c(IID.x.x="IID", SCORE.x.x="SCORE500", SCORE.y.x="SCORE150", SCORE.x.y="SCORE350"))
summary(scores) #Check Ns remain the same
#Note. SCORE150 becomes N-94896 because exclusions are now applied

#Save genetic data on RDSF
#Connect to RDSF on server
setwd("pathname")
write.csv(scores, "UKBiobank_BMI_PRS_97.csv", quote=F, row.names=F)