#This script reads in the UK Biobank Phenotype data and subsests/formats the required phenotypes
#Robyn Wootton 
######################################################################################
#Contents
#1. Load Phenotype Data
#2. Tidy data
#3. Make missing data NA
#4. Reverse scores (high score = high wellbeing)
#5. Create mean scores across the three time points
#6. Save the formatted phenotype file

######################################################################################
#1. Load phenotype data
######################################################################################
#Connext to IEU server using finder-go-connect to server
setwd("pathname")
#install.packages("data.table")
library(data.table)
read <- fread("data.6564.csv", select = c("eid", "31-0.0", "21022-0.0", "189-0.0", "21000-0.0", "1558-0.0", "4526-0.0", "4526-1.0", "4526-2.0", "4537-0.0", "4537-1.0", "4537-2.0", "4548-0.0", "4548-1.0", "4548-2.0", "4581-0.0", "4581-1.0", "4581-2.0", "4570-0.0", "4570-1.0", "4570-2.0", "4559-0.0", "4559-1.0", "4559-2.0", "21001-0.0", "21001-1.0", "21001-2.0"), sep=",")
phen<-read
str(phen)

######################################################################################
#2. Tidy data
######################################################################################
#Rename cols so they are not numbers
library(reshape)
phen<-rename(phen, c("31-0.0"="sex", "21022-0.0"="Age", "189-0.0"="SES", "21000-0.0"="ethnic", "1558-0.0"="alcohol", "4526-0.0"="happ1", "4526-1.0"="happ2", "4526-2.0"="happ3", "4537-0.0"="work1", "4537-1.0"="work2", "4537-2.0"="work3", "4548-0.0"="health1", "4548-1.0"="health2", "4548-2.0"="health3", "4581-0.0"="finance1", "4581-1.0"="finance2", "4581-2.0"="finance3", "4570-0.0"="friends1", "4570-1.0"="friends2", "4570-2.0"="friends3", "4559-0.0"="family1", "4559-1.0"="family2", "4559-2.0"="family3", "21001-0.0"="BMI1", "21001-1.0"="BMI2", "21001-2.0"="BMI3"))
str(phen)


#Make all variables numeric
phen<-data.frame(phen)
cols <- c("eid", "sex", "Age", "SES", "ethnic", "alcohol", "happ1", "happ2", "happ3", "work1", "work2", "work3", "health1", "health2", "health3", "finance1", "finance2", "finance3", "friends1", "friends2", "friends3", "family1", "family2", "family3", "BMI1", "BMI2", "BMI3");    
phen[,cols] <- apply(phen[,cols], 2, function(x) as.numeric(as.character(x)));
str(phen)


######################################################################################
#3. Make missing data NA
######################################################################################
#Make all missing data NA
table(phen[,5])
phen[phen == ""] <-NA
table(is.na(phen[,5]))


#make 'I don't know' and 'Prefer not to say' into NA
phen$happ1[phen$happ1 == -1] <-NA
phen$happ1[phen$happ1 == -3] <-NA
phen$happ2[phen$happ2 == -1] <-NA
phen$happ2[phen$happ2 == -3] <-NA
phen$happ3[phen$happ3 == -1] <-NA
phen$happ3[phen$happ3 == -3] <-NA
phen$work1[phen$work1 == -1] <-NA
phen$work1[phen$work1 == -3] <-NA
phen$work1[phen$work1 == 7] <-NA #This is response 'I am not employed'
phen$work2[phen$work2 == -1] <-NA
phen$work2[phen$work2 == -3] <-NA
phen$work2[phen$work2 == 7] <-NA 
phen$work3[phen$work3 == -1] <-NA
phen$work3[phen$work3 == -3] <-NA
phen$work3[phen$work3 == 7] <-NA 
phen$health1[phen$health1 == -1] <-NA
phen$health1[phen$health1 == -3] <-NA
phen$health2[phen$health2 == -1] <-NA
phen$health2[phen$health2 == -3] <-NA
phen$health3[phen$health3 == -1] <-NA
phen$health3[phen$health3 == -3] <-NA
phen$finance1[phen$finance1 == -1] <-NA
phen$finance1[phen$finance1 == -3] <-NA
phen$finance2[phen$finance2 == -1] <-NA
phen$finance2[phen$finance2 == -3] <-NA
phen$finance3[phen$finance3 == -1] <-NA
phen$finance3[phen$finance3 == -3] <-NA
phen$friends1[phen$friends1 == -1] <-NA
phen$friends1[phen$friends1 == -3] <-NA
phen$friends2[phen$friends2 == -1] <-NA
phen$friends2[phen$friends2 == -3] <-NA
phen$friends3[phen$friends3 == -1] <-NA
phen$friends3[phen$friends3 == -3] <-NA
phen$family1[phen$family1 == -1] <-NA
phen$family1[phen$family1 == -3] <-NA
phen$family2[phen$family2 == -1] <-NA
phen$family2[phen$family2 == -3] <-NA
phen$family3[phen$family3 == -1] <-NA
phen$family3[phen$family3 == -3] <-NA
summary(phen)

######################################################################################
#4. Reverse scores (high score = high wellbeing)
######################################################################################
#reverse score all of the wellbeing variables so that 6 is high wellbeing and 1 is low
phen$happ1<-(7-phen$happ1)
phen$happ2<-(7-phen$happ2)
phen$happ3<-(7-phen$happ3)
phen$work1<-(7-phen$work1)
phen$work2<-(7-phen$work2)
phen$work3<-(7-phen$work3)
phen$health1<-(7-phen$health1)
phen$health2<-(7-phen$health2)
phen$health3<-(7-phen$health3)
phen$finance1<-(7-phen$finance1)
phen$finance2<-(7-phen$finance2)
phen$finance3<-(7-phen$finance3)
phen$friends1<-(7-phen$friends1)
phen$friends2<-(7-phen$friends2)
phen$friends3<-(7-phen$friends3)
phen$family1<-(7-phen$family1)
phen$family2<-(7-phen$family2)
phen$family3<-(7-phen$family3)
summary(phen) #check means have flipped

######################################################################################
#5. Create mean scores across the three time points
######################################################################################
#Create mean scores of each of the variables
str(phen)
#assign the column names
happ123<-c("happ1", "happ2", "happ3")
work123<-c("work1", "work2", "work3")
health123<-c("health1", "health2", "health3")
finance123<-c("finance1", "finance2", "finance3")
friends123<-c("friends1", "friends2", "friends3")
family123<-c("family1", "family2", "family3")
BMI123<-c("BMI1", "BMI2", "BMI3")

#create the row means
phen$happ <- apply(phen[, happ123],1,mean, na.rm=T)
phen$work <- apply(phen[, work123],1,mean, na.rm=T)
phen$health <- apply(phen[, health123],1,mean, na.rm=T)
phen$finance <- apply(phen[, finance123],1,mean, na.rm=T)
phen$friends <- apply(phen[, friends123],1,mean, na.rm=T)
phen$family <- apply(phen[, family123],1,mean, na.rm=T)
phen$BMI <- apply(phen[, BMI123],1,mean, na.rm=T)
head(phen, n=50)

######################################################################################
#6. Save the formatted phenotype file
######################################################################################
#Save phenotypes to RDSF
#Connect to RDSF on server
setwd("pathname")
write.csv(phen, "UKBiobank_BMIandWellbeing_Phenotypes.csv", quote=F, row.names=F)
