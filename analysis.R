##############################################################################
# Script name: analysis.R
# Purpose:
# Author: 
# Date Created: 2021-08-24
# Copyright (c) Phillip Hungerford, 2020
# Email: phillip.hungerford@gmail.com
##############################################################################
# Notes:
#
#
##############################################################################
# set working directory for Mac and PC
#setwd("~/Google Drive/")
#setwd("C:/Users/tim/Google Drive/")

##############################################################################
options(scipen = 6, digits = 4) 
#memory.limit(30000000)

##############################################################################
# Dependencies:
library(dplyr) # for filtering and data manipulation
library(ggplot2) # for data visualisation
library(survival) # for survival analysis
##############################################################################
# LOAD DATA

# SOURCE: https://healthdatainsight.org.uk/project/the-simulacrum/
# Load patient level data
sim_av_patient <- read.csv("simulacrum_release_v1.2.0.2017/data/sim_av_patient.csv", na.strings = c(""))

# Check variables match Dictionary
names(sim_av_patient)

# preview data
head(sim_av_patient)

# examine data
str(sim_av_patient)

# format
sim_av_patient %>% 
  mutate(
    VITALSTATUSDATE = as.Date(VITALSTATUSDATE , format = "%Y-%m-%d")
    )

# clean vital status
zvitalstatus_lookup <- read.csv('simulacrum_release_v1.2.0.2017/data_dictionary_files/zvitalstatus_lookup.csv')

# create event variable
sim_av_patient$event <- sim_av_patient$NEWVITALSTATUS

for (i in 1:nrow(zvitalstatus_lookup)){
  sim_av_patient$NEWVITALSTATUS[sim_av_patient$NEWVITALSTATUS == zvitalstatus_lookup$ZVITALSTATUSID[i]] <- zvitalstatus_lookup$SHORTDESC[i]
}

table(sim_av_patient$NEWVITALSTATUS)

# clean age
zsex_lookup <- read.csv('simulacrum_release_v1.2.0.2017/data_dictionary_files/zsex_lookup.csv')

for (i in 1:nrow(zsex_lookup)){
  sim_av_patient$SEX[sim_av_patient$SEX == zsex_lookup$ZSEXID[i]] <- zsex_lookup$SHORTDESC[i]
}

# clean ethnicity
zethnicity_lookup <- read.csv('simulacrum_release_v1.2.0.2017/data_dictionary_files/zethnicity_lookup.csv')

for (i in 1:nrow(zethnicity_lookup)){
  sim_av_patient$ETHNICITY[sim_av_patient$ETHNICITY == zethnicity_lookup$ZETHNICITYID[i]] <- zethnicity_lookup$SHORTDESC[i]
}

# make factors
sim_av_patient <- sim_av_patient %>% mutate_at(c('SEX',
                                                 'ETHNICITY',
                                                 'NEWVITALSTATUS'), as.factor)


summary(sim_av_patient)

#=============================================================================
# Load tumour level data
sim_av_tumour <- read.csv("simulacrum_release_v1.2.0.2017/data/sim_av_tumour.csv")

# Check variables match Dictionary
names(sim_av_tumour)

# preview data
head(sim_av_tumour)

# examine data
str(sim_av_tumour)

# Format data
sim_av_tumour %>% 
  mutate(
    DIAGNOSISDATEBEST  = as.Date(DIAGNOSISDATEBEST , format = "%Y-%m-%d")
  )
# Format the data
## CODE

##############################################################################
# Results
#=============================================================================
# 1. Analysis

#-----------------------------------------------------------------------------
# a. What is the sex distribution of patients?
sim_av_patient %>%
  count(SEX) %>%
  # count creates a column called 'n'
  mutate(percent = n / sum(n) * 100)

#-----------------------------------------------------------------------------
# b. Prepare a dataset with all lung cancer patients (SITE_ICD10_O2_3CHAR= C34) and their cause of deaths.

# get earliest diagnosis date for patients
first_diagnosis <- sim_av_tumour %>% 
  select(PATIENTID, DIAGNOSISDATEBEST) %>% 
  group_by(PATIENTID) %>%
  filter(DIAGNOSISDATEBEST == min(DIAGNOSISDATEBEST)) %>% distinct()

# there seem to be more patients

# lets find the duplicated patient data
duplicate_ids <- first_diagnosis$PATIENTID[duplicated(first_diagnosis$PATIENTID)]

# view data
duplicate_df <- first_diagnosis[first_diagnosis$PATIENTID %in% t, ]

# inspect
head(t2)

# we just remove exact duplicates
first_diagnosis <- first_diagnosis[!duplicated(first_diagnosis),] # N = 220062

# merge patient diagnosis data to patient level data by patient id
analysis_data <- left_join(sim_av_patient, first_diagnosis,by="PATIENTID")

analysis_data <- analysis_data %>% mutate(
  DIAGNOSISDATEBEST = as.Date(DIAGNOSISDATEBEST, format = "%Y-%m-%d"),
  VITALSTATUSDATE = as.Date(VITALSTATUSDATE, format = "%Y-%m-%d"),
  time = as.numeric(difftime(VITALSTATUSDATE, DIAGNOSISDATEBEST, units = "days"))
)
# Extract lung cancer patients
target <- "C34"

# find patients with target from tumour data
lung_cancer_tumours <- sim_av_tumour %>% filter(SITE_ICD10_O2_3CHAR == target)

# get patient Ids
lung_cancer_patient_ids <- unique(lung_cancer_tumours$PATIENTID) # n = 167233

# extract these patients
lung_cancer_patients <- analysis_data[analysis_data$PATIENTID %in% lung_cancer_patient_ids, ]

# check missing for key variables
na_count <-sapply(lung_cancer_patients, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)
print(na_count)

# for survival analysis we need time

#-----------------------------------------------------------------------------
# c. Calculate the mean survival time in days for the lung cancer patients who were diagnosed in 2013.

# filter patients that were diagnosed in 2013 only as per condition.
lung_cancer_patients_2013 <- lung_cancer_patients %>% filter(DIAGNOSISDATEBEST >= "2013-01-01" & DIAGNOSISDATEBEST <= "2013-12-31") # N = 37081

table(lung_cancer_patients_2013$NEWVITALSTATUS)

# Alive = 9358        
# Died = 27707
# Exit posting = 16

lung_cancer_patients_2013$status <- NA
lung_cancer_patients_2013$status[lung_cancer_patients_2013$NEWVITALSTATUS == "Alive"] <- 0
lung_cancer_patients_2013$status[lung_cancer_patients_2013$NEWVITALSTATUS == "Dead"] <- 1
lung_cancer_patients_2013$status[lung_cancer_patients_2013$NEWVITALSTATUS == "Exit Posting"] <- 0 #!!

# Survival from diagnosis -> death

# Kaplan-Meier estimate can be obtained as numbers by using survfit 
om.fit <- survfit(Surv(time, status) ~ 1, data = lung_cancer_patients_2013) 
summary(om.fit)

## Report how many deaths were there during the follow-up and what was the median follow-up time until death.
# get the mean and median survival times
print(om.fit, print.rmean = TRUE) #1976

# plotting Nelson-Aalen estimate
om.fit <- survfit( Surv(time, status) ~ 1, data = lung_cancer_patients_2013) 
summary(om.fit)
plot(om.fit, fun="cumhaz", main = "Nelson-Aalen estimate of cumulative hazard function")


# Plot the Kaplan-Meier survival curves by all categorical variables in the data.
om.fit3 <- survfit(Surv(time, status) ~ SEX, data = lung_cancer_patients_2013) 

# plotting Kaplan-Meier estimate
plot(om.fit3, col=c('blue', 'red', 'green'), main = 'Kaplan-Meier estimate for Sex')

# Carry out the log-rank tests for differences in survival experience by the categorical variables and interpret the results.
# 2-sample log-rank test
om.fit2b <- coxph( Surv(time, status) ~ SEX, data = lung_cancer_patients_2013) 
summary(om.fit2b)

# fitting a simple Cox model for sex
cox.rx <- coxph( Surv(time, status) ~ SEX, data = lung_cancer_patients_2013) 
summary(cox.rx) #sig p=0.003


#=============================================================================
# 2. Visualisation
#-----------------------------------------------------------------------------
# a. Visualise the proportion of all ethnic groups among non-British male patients. You can find ethnic information of patients in the data dictionary provided in the simulacrum dataset.
#-----------------------------------------------------------------------------
# b. Explore the data to find an interesting or insightful aspect that you can then communicate visually. Be creative!
##############################################################################
#################################### END #####################################
##############################################################################