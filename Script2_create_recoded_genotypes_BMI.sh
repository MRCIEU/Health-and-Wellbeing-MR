#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -me
#PBS -l nodes=1:ppn=8
#
#---------------------------------------------
#
cd "${HOME}/path/to/my/directory"
#
#---------------------------------------------

export DIR1="path/to/UKBiobank/genetic/data"
export DIR2="path/to/my/directory"


# This scripts creates a file with genotypes coded as 0/1/2 for each of the snps 

plink --gen $DIR2/snps-out.gen --sample $DIR1/sample-stats/data.chr01.sample  --remove $DIR2/biobank_exclusion_list_unique_ids.txt --allow-extra-chr --extract $DIR2/BMISNPs.txt --recode compound-genotypes --out $DIR2/add_snps_02Oct2017


 
