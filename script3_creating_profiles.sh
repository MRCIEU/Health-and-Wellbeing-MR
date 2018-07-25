#program to generate genomewide score profiles based on SNP beta values 
#Plink v1.09 is needed to run the extractsnps.sh script called in this program 


#!bin/bash

#################################################
#only needed for qsub'd runs
#set time 
#PBS -l nodes=1:ppn=1,walltime=8:00:00

#Define working directory 
export WORK_DIR=$HOME/path/to/my/directory

#Change into working directory 
cd $WORK_DIR

#store the date
now=$(date +"%Y%m%d")

#################################################

#generate scores
 plink --file add_snps_02Oct2017 \
 		--noweb \
 		--allow-extra-chr \
 		--score BMI_betas.txt \
 		no-mean-imputation \
 		--out ${now}_scores_BMI


