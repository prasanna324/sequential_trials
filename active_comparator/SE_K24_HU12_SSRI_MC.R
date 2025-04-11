# Load required packages

if (!require("tableone")) install.packages("tableone")
library(tableone)
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
if (!require("boot")) install.packages("boot")
library(boot)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("data.table")) install.packages("data.table")
library(data.table)
if (!require("splitstackshape")) install.packages("splitstackshape")
library("splitstackshape")
if (!require("Hmisc")) install.packages("Hmisc")
library(Hmisc)
if (!require("speedglm")) install.packages("speedglm")
library(speedglm)
if (!require("tidyverse")) install.packages("tidyverse")
library(tidyverse)
if (!require("arrow")) install.packages("arrow")
library(arrow)
if (!require("ggpubr")) install.packages("ggpubr")
library(ggpubr)
if (!require("zoo")) install.packages("zoo")
library(zoo)
if (!require("lubridate")) install.packages("lubridate")
library(lubridate)

setwd('~/Medications')

# Set healthcare utililization (HU) censoring time, event name, and time point for estimation of risks in months
HU_time <- 12
event_name <- 'MC'
K <- 24

# Import long format data file with covariate and outcome data
med_obs <- read_feather(paste0('ETT/SSRI/Processed_files/SSRI_longdata_', event_name,'.ft'))
med_obs$id <- med_obs$Pid
length(unique(med_obs$Pid))

# Set eligibility and encounter dates based on censoring time for HU
med_obs$elig <- med_obs$elig12
med_obs$encounter_date <- med_obs$encounter_date12
med_obs$encounter_e_date <- med_obs$encounter_e_date12

med_obs <- med_obs %>%
  dplyr::select(1:10, elig, 11:ncol(med_obs))
med_obs <- med_obs %>%
  dplyr::select(1:3, encounter_date, 4:ncol(med_obs))

# Create a time varying continuous variable for indication
med_obs$ind_SSRI_con <- as.numeric(difftime(med_obs$current_date, med_obs$indication_date, units = "days")+1)/365.25

# Replace NAs with 0 in the ind_SSRI_con column
med_obs$ind_SSRI_con[is.na(med_obs$ind_SSRI_con)] <- 0

# Baseline indication time
med_obs <- med_obs %>%
  arrange(Pid, cal_time) %>%
  group_by(Pid) %>%
  mutate(ind_SSRI_con_sb = if_else(cal_time == 0, ind_SSRI_con, NA_integer_)) %>%
  fill(ind_SSRI_con_sb, .direction = "down") %>%
  ungroup()

# Create a replicate of encounter date but with no NAs. This will be used for computing time since last HU
# Fill NAs in encounter_date column with the last available non-NA value
med_obs <- med_obs %>%
  group_by(Pid) %>%
  mutate(enc_date_nomiss = zoo::na.locf(encounter_date))

# Calculate time difference in months between current_date and encounter_date
med_obs$time_lastHU <- round(as.numeric(interval(med_obs$enc_date_nomiss, med_obs$current_date) / months(1)), 2)

med_obs$time_lastHU_days <- round(as.numeric(interval(med_obs$enc_date_nomiss, med_obs$current_date) / days(1)), 2)

med_obs <- med_obs %>%
  dplyr::select(1:3, time_lastHU, 4:ncol(med_obs))

med_obs <- med_obs %>%
  dplyr::select(1:4, time_lastHU_days, 5:ncol(med_obs))

### Preprocess variables (create dummies for categorical covariates) ###

### Calendar Year

# Create the calendar Year variable
med_obs <- med_obs %>%
  mutate(date_admission = as.Date(date_admission),
         Year = as.numeric(format(date_admission, "%Y")))

# Define the breaks for the categories and create a factor variable
breaks <- c(2005, 2008, 2011, 2014, 2017)
med_obs$Year_cat <- cut(med_obs$Year, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$Year_cat <- as.factor(med_obs$Year_cat)

# Create dummy variables
dummy_variables <- model.matrix(~ Year_cat - 1, data = med_obs)
colnames(dummy_variables) <- c("Year_2006-2008", "Year_2009-2011", "Year_2012-2014", "Year_2015-2017")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

### Age

# Define the breaks for the categories and create a factor variable
breaks <- c(64, 70, 75, 80, 85, Inf)
med_obs$Age_cat <- cut(med_obs$Age, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$Age_cat <- as.factor(med_obs$Age_cat)

# Create dummy variables
dummy_variables <- model.matrix(~ Age_cat - 1, data = med_obs)
colnames(dummy_variables) <- c("Age_65-70","Age_71-75","Age_76-80", "Age_81-85", "Age_86+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

### Sex

# Re-factoring the Sex variable
med_obs <- med_obs %>%
  mutate(Sex = ifelse(Sex == 1, 0, ifelse(Sex == 2, 1, Sex)))

### Education

# Recode the levels of the educ variable
med_obs$educ <- ifelse(med_obs$educ %in% c(0, 1), 0, 
                       ifelse(med_obs$educ == 2, 1, 
                              ifelse(med_obs$educ == 3, 2, med_obs$educ)))

med_obs$educ <- as.factor(med_obs$educ)
# Create dummy variables
dummy_variables <- model.matrix(~ educ - 1, data = med_obs)

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

### Med indication

# Define the breaks for the categories and create a factor variable
breaks_med <- c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, Inf)
med_obs$ind_cat <- cut(med_obs$ind_SSRI_con, breaks = breaks_med, include.lowest = FALSE, right = TRUE)
med_obs$ind_cat <- as.factor(med_obs$ind_cat)

# Create dummy variables
dummy_variables_med <- model.matrix(~ ind_cat - 1, data = med_obs)
colnames(dummy_variables_med) <- c("ind_SSRI_0", "ind_SSRI_1", "ind_SSRI_2", "ind_SSRI_3", "ind_SSRI_4", "ind_SSRI_5",
                                   "ind_SSRI_6", "ind_SSRI_7", "ind_SSRI_8", "ind_SSRI_9+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables_med)
rm(dummy_variables_med)

# Baseline
med_obs$ind_cat_sb <- cut(med_obs$ind_SSRI_con_sb, breaks = breaks_med, include.lowest = FALSE, right = TRUE)
med_obs$ind_cat_sb <- as.factor(med_obs$ind_cat_sb)

# Create dummy variables
dummy_variables_med <- model.matrix(~ ind_cat_sb - 1, data = med_obs)
colnames(dummy_variables_med) <- c("ind_SSRI_sb_0", "ind_SSRI_sb_1", "ind_SSRI_sb_2", "ind_SSRI_sb_3", "ind_SSRI_sb_4", "ind_SSRI_sb_5",
                                   "ind_SSRI_sb_6", "ind_SSRI_sb_7", "ind_SSRI_sb_8", "ind_SSRI_sb_9+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables_med)
rm(dummy_variables_med)

### CCI

# Define the breaks for the categories and create a factor variable
breaks <- c(-1, 0, 2, 4, Inf)
med_obs$CCIW_cat <- cut(med_obs$CCIW, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$CCIW_cat <- as.factor(med_obs$CCIW_cat)

# Create dummy variables
dummy_variables <- model.matrix(~ CCIW_cat - 1, data = med_obs)
colnames(dummy_variables) <- c("CCI_0","CCI_1-2","CCI_3-4","CCI_5+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

# Baseline
med_obs$CCIW_cat_sb <- cut(med_obs$CCIW_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$CCIW_cat_sb <- as.factor(med_obs$CCIW_cat_sb)

# Create dummy variables
dummy_variables <- model.matrix(~ CCIW_cat_sb - 1, data = med_obs)
colnames(dummy_variables) <- c("CCI_sb_0","CCI_sb_1-2","CCI_sb_3-4","CCI_sb_5+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)


### N endo
breaks <- c(-1, 0, Inf)
med_obs$N.endo_cat <- cut(med_obs$N.endo, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.endo_cat <- as.factor(med_obs$N.endo_cat)

# Create dummy variables
dummy_variables <- model.matrix(~ N.endo_cat - 1, data = med_obs)
colnames(dummy_variables) <- c("N.endo_0","N.endo_1+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

# Baseline
med_obs$N.endo_cat_sb <- cut(med_obs$N.endo_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.endo_cat_sb <- as.factor(med_obs$N.endo_cat_sb)

# Create dummy variables
dummy_variables <- model.matrix(~ N.endo_cat_sb - 1, data = med_obs)
colnames(dummy_variables) <- c("N.endo_sb_0","N.endo_sb_1+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

### N GI
breaks <- c(-1, 0, Inf)
med_obs$N.GI_cat <- cut(med_obs$N.GI, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.GI_cat <- as.factor(med_obs$N.GI_cat)

# Create dummy variables
dummy_variables <- model.matrix(~ N.GI_cat - 1, data = med_obs)
colnames(dummy_variables) <- c("N.GI_0","N.GI_1+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

# Baseline
med_obs$N.GI_cat_sb <- cut(med_obs$N.GI_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.GI_cat_sb <- as.factor(med_obs$N.GI_cat_sb)

# Create dummy variables
dummy_variables <- model.matrix(~ N.GI_cat_sb - 1, data = med_obs)
colnames(dummy_variables) <- c("N.GI_sb_0","N.GI_sb_1+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)


### N visits
breaks <- c(-1, 0, 5, 10, Inf)
med_obs$N.visits_cat <- cut(med_obs$N.visits, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.visits_cat <- as.factor(med_obs$N.visits_cat)

# Create dummy variables
dummy_variables <- model.matrix(~ N.visits_cat - 1, data = med_obs)
colnames(dummy_variables) <- c("N.visits_0","N.visits_1-5","N.visits_6-10", "N.visits_11+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

# Baseline
med_obs$N.visits_cat_sb <- cut(med_obs$N.visits_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.visits_cat_sb <- as.factor(med_obs$N.visits_cat_sb)

# Create dummy variables
dummy_variables <- model.matrix(~ N.visits_cat_sb - 1, data = med_obs)
colnames(dummy_variables) <- c("N.visits_sb_0","N.visits_sb_1-5","N.visits_sb_6-10", "N.visits_sb_11+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)


### N drugs
breaks <- c(-1, 0, 4, 8, 12, 16, Inf)
med_obs$N.drugs_cat <- cut(med_obs$N.drugs, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.drugs_cat <- as.factor(med_obs$N.drugs_cat)

# Create dummy variables
dummy_variables <- model.matrix(~ N.drugs_cat - 1, data = med_obs)
colnames(dummy_variables) <- c("N.drugs_0","N.drugs_1-4", "N.drugs_5-8", "N.drugs_9-12", "N.drugs_13-16", "N.drugs_17+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

# Basline
med_obs$N.drugs_cat_sb <- cut(med_obs$N.drugs_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.drugs_cat_sb <- as.factor(med_obs$N.drugs_cat_sb)

# Create dummy variables
dummy_variables <- model.matrix(~ N.drugs_cat_sb - 1, data = med_obs)
colnames(dummy_variables) <- c("N.drugs_sb_0","N.drugs_sb_1-4", "N.drugs_sb_5-8", "N.drugs_sb_9-12", "N.drugs_sb_13-16", "N.drugs_sb_17+")

# Append the dummy variables to the original dataframe
med_obs <- cbind(med_obs, dummy_variables)
rm(dummy_variables)

##############################################
##############################################

# Create HU censoring variable
med_obs <- med_obs %>%
  group_by(id) %>%
  mutate(censor_hu = ifelse(is.na(encounter_date), 1, 0)) %>%
  ungroup()

# Update 'elig' to 0 where 'censor_hu' is equal to 1
med_obs$elig[med_obs$censor_hu == 1] <- 0

# Add new eligibility criteria of HU within the past month
# med_obs$elig[med_obs$time_lastHU_days > 31] <- 0

# Set elig to 0 for rows where encounter_date is not between current_date and current_date_end
med_obs$elig <- ifelse(!is.na(med_obs$encounter_e_date) &
                         (med_obs$encounter_e_date >= med_obs$current_date &
                            med_obs$encounter_e_date <= med_obs$current_date_end),
                       med_obs$elig, 0)

# Set elig to 0 when no visits or prescriptions are present in the past 12 months
med_obs$elig <- ifelse(med_obs$N.visits == 0 & med_obs$N.drugs == 0, 0, med_obs$elig)

# Update censoring due to loss of follow-up variable
med_obs <- med_obs %>%
  group_by(id) %>%
  mutate(censor_enc = ifelse(N.visits == 0 & N.drugs == 0, 1, 0)) %>%
  ungroup()

# For active comparator cohorts, "elig" can only be 1 in the first trial
med_obs$elig <- ifelse(med_obs$cal_time == 0, med_obs$elig, 0)

# Create indicator for baseline treatment group
med_obs$treat_b <- med_obs$med_type=="SSRI"
med_obs$treat_b <- as.numeric(med_obs$treat_b)


#################################################################################
#################################################################################

### Perform sequential emulation ###

# Write a function to create sequentially emulated dataset
seq.em <- function(data, trial.num) {
  datalist <- list()
  for(i in (1:trial.num)){
    # Extract list of eligible IDs for trial i
    ids.trial <- data$id[which(data$elig==1 & data$cal_time==(i-1))]
    # Create person-time data for trial i
    trial <- data[data$id %in% ids.trial, ]
    # Create follow-up time and trial number variables
    trial$trial_num <- i - 1
    trial$time <- trial$cal_time - (i - 1)
    trial$timesqr <- trial$time * trial$time
    # Create concatenated unique ID number
    trial$id_new <- paste(trial$id, "-", i, sep = "")
    # Delete person-time rows before the start of trial i
    # + administratively censor at 24 weeks after baseline
    trial <- trial[(trial$time >= 0 & trial$time <= K-1), ]
    # Add trial i to the list
    datalist[[i]] <- trial
  }
  # Combine all sequentially emulated trials into a single dataset
  all <- rbindlist(datalist)
  return(all)
}

# Define the start and end dates
start_date <- as.Date("2006-01-01")
end_date <- as.Date("2017-12-31")

# Calculate the difference in months
numberofmonths <- as.numeric(interval(start_date, end_date) %/% months(1))

# Sequentially emulating trials
trials_all <- seq.em(data = med_obs, trial.num = numberofmonths)

trials_all <- trials_all %>%
  dplyr::select(1:10, id_new, 11:ncol(trials_all))


#################################################################################
#################################################################################

##### Create censoring variables ######

# Censoring due to no healthcare utilization in the past 24 months
trials_all <- trials_all %>%
  group_by(id_new) %>%
  mutate(censor_hu = ifelse(cumsum(censor_hu)>=1, 1, 0)) %>%
  ungroup()

trials_all <- trials_all %>%
  dplyr::select(1:17, censor_hu, 18:ncol(trials_all))

# Censoring due to loss of follow-up
trials_all <- trials_all %>%
  group_by(id_new) %>%
  mutate(censor_enc = ifelse(cumsum(censor_enc)>=1, 1, 0)) %>%
  ungroup()


#### Counts ####

# Count the number of unique ids where censor_hu is equal to 1
unique_hu_censor <- trials_all %>%
  filter(censor_hu == 1) %>%
  summarise(unique_ids = n_distinct(id))

cat('Censoring due to no health care utilization:', unique_hu_censor$unique_ids, "\n")
rm(unique_hu_censor)

# Count the number of unique ids where IBD equal to 1 in sequential emulation
unique_ibd <- trials_all %>%
  filter(current_date > First_IBD_DATE) %>%
  summarise(unique_ids = n_distinct(id))

cat('Censoring due to IBD:', unique_ibd$unique_ids, "\n")
rm(unique_ibd)

#################################################################################
#################################################################################

### Create an indicator for ever being censored due to lack of healthcare utilization
trials_all <- trials_all %>%
  dplyr::group_by(id_new) %>%
  dplyr::mutate(
    censor_hu_any = ifelse(is.na(censor_hu), NA, max(censor_hu, na.rm = T))) %>%
  dplyr::ungroup()

### Describe number of individuals lost to follow-up ###

# Count & proportion of unique individuals censored
count_censor_hu <- trials_all[which(trials_all$censor_hu_any==1),]
no_dups_censor_hu <- count_censor_hu[!duplicated(count_censor_hu$id),]
no_dups_censor_hu %>% dplyr::summarise(
  n_censored_hu = length(censor_hu_any[censor_hu_any == 1]),
  prop_censored_hu = length(censor_hu_any[censor_hu_any == 1])/length(trials_all[!duplicated(trials_all$id),]$id))

rm(count_censor_hu)

#################################################################################
#################################################################################

#### Exclude censored rows ####

trials_all <- trials_all %>% filter(!(censor_enc==1 | censor_hu==1 | death==1))

trials_all <- trials_all %>%
  group_by(id_new) %>%
  mutate(month_num = row_number() - 1) %>%
  ungroup()

# Create dataset restricted to baseline (follow-up time = 0)
trials_all_base <- trials_all %>%
  dplyr::filter(time == 0)

print('Adding baseline data to person-time data')

# Renaming the long dataset to be consistent with other cohorts
med_elig <- trials_all

rm(trials_all)
gc()

# Create indicator of any MC occurrence over the 24-week follow-up period
med_elig <- med_elig %>%
  dplyr::group_by(id_new) %>%
  dplyr::mutate(
    MC_any = ifelse(is.na(eventMC), NA, max(eventMC, na.rm = T))) %>%
  dplyr::ungroup()

table(med_elig$MC_any[which(med_elig$time==0)], useNA = "always")


### Create a descriptive table ###

# List of variables to be included in the table

myVars <- c("Age", "Age_65-70","Age_71-75","Age_76-80", "Age_81-85", "Age_86+",
            "Sex", "nordic",
            "Year", "Year_2006-2008", "Year_2009-2011", "Year_2012-2014", "Year_2015-2017",
            "educ0", "educ1", "educ2", 
            "N.GI_sb_1+", "N.endo_sb_1+",
            "N.visits_sb_0","N.visits_sb_1-5","N.visits_sb_6-10", "N.visits_sb_11+",
            "N.drugs_sb_0","N.drugs_sb_1-4", "N.drugs_sb_5-8", "N.drugs_sb_9-12", "N.drugs_sb_13-16", "N.drugs_sb_17+",
            "CCI_sb_0","CCI_sb_1-2","CCI_sb_3-4","CCI_sb_5+",
            "ind_SSRI_sb_0", "ind_SSRI_sb_1", "ind_SSRI_sb_2", "ind_SSRI_sb_3", "ind_SSRI_sb_4", "ind_SSRI_sb_5",
            "ind_SSRI_sb_6", "ind_SSRI_sb_7", "ind_SSRI_sb_8", "ind_SSRI_sb_9+")


# List of categorical variables
catVars <- c("Age_65-70","Age_71-75","Age_76-80", "Age_81-85", "Age_86+",
             "Sex", "nordic",
             "Year", "Year_2006-2008", "Year_2009-2011", "Year_2012-2014", "Year_2015-2017",
             "educ0", "educ1", "educ2", 
             "N.GI_sb_1+", "N.endo_sb_1+",
             "N.visits_sb_0","N.visits_sb_1-5","N.visits_sb_6-10", "N.visits_sb_11+",
             "N.drugs_sb_0","N.drugs_sb_1-4", "N.drugs_sb_5-8", "N.drugs_sb_9-12", "N.drugs_sb_13-16", "N.drugs_sb_17+",
             "CCI_sb_0","CCI_sb_1-2","CCI_sb_3-4","CCI_sb_5+",
             "ind_SSRI_sb_0", "ind_SSRI_sb_1", "ind_SSRI_sb_2", "ind_SSRI_sb_3", "ind_SSRI_sb_4", "ind_SSRI_sb_5",
             "ind_SSRI_sb_6", "ind_SSRI_sb_7", "ind_SSRI_sb_8", "ind_SSRI_sb_9+")

# List of continuous variables which should be displayed as median (IQR)
medVars <- c("Age")

# Create table
tab1 <- CreateTableOne(vars = myVars, # set descriptive variables
                       strata = "treat_b", # define stratifying variable
                       data = trials_all_base, # baseline
                       factorVars = catVars) # define categorical variables

# Print table
print(tab1,
      nonnormal = medVars,
      formatOptions = list(big.mark = ","),
      test = FALSE)


#################################################################################
#################################################################################

### Describe the number of events ###

# Create a baseline dataset without duplicate IDs
no_dups <- med_elig[!duplicated(med_elig$id),]

# Number of unique individuals
length(no_dups$id)

# Number of unique events
MC_count <- med_elig[which(med_elig$MC_any==1),]
no_dups_MC <- MC_count[!duplicated(MC_count$id),]
table(no_dups_MC$MC_any)
rm(MC_count)

# Compute average number of trials each individual contributed to
avg_trials <- trials_all_base %>% dplyr::group_by(id) %>%
  dplyr::summarise(trials_count = length(id))

# Print average number of trials
summary(avg_trials$trials_count)

#################################################################################
#################################################################################

### Fit a pooled logistic model that predicts the probability of remaining uncensored
# (due to lack of HU) at each time ###

no_dups <- med_elig[!duplicated(med_elig$id),]
med_obs_elig_ids <- no_dups$id
med_obs_elig <- med_obs[med_obs$id %in% med_obs_elig_ids, ]

med_obs_elig_subset <- med_obs_elig %>%
  filter(time_lastHU >= HU_time & time_lastHU < HU_time+1)

# Fit the pooled logistic model

Knots_year <- Hmisc::rcspline.eval(med_obs_elig_subset$Year, knots = quantile(med_obs_elig_subset$Year, c(0.05, .35, .65, .95)),
                                   nk = 4, knots.only = TRUE)

psc.denom <- speedglm(censor_hu == 0 ~ treat_b + cal_time  + cal_timesqr + Hmisc::rcspline.eval(Year,  knots = Knots_year, nk = 4) +
                        as.factor(Age_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) +
                        as.factor(N.endo_cat_sb) + as.factor(N.GI_cat_sb) +  as.factor(N.visits_cat_sb) +
                        as.factor(N.drugs_cat_sb) + as.factor(CCIW_cat_sb) + as.factor(ind_cat_sb) + as.factor(CCIW_cat),
                      family=binomial(link="logit"),
                      data=med_obs_elig_subset)

summary(psc.denom)
gc()

# Obtain predicted probabilities of being uncensored in the original data
med_obs_elig$psc.denom <- predict(psc.denom, med_obs_elig, type = "response")
gc()

summary(med_obs_elig$psc.denom)
quantile(med_obs_elig$psc.denom, probs = 0.95)

#################################################################################
#################################################################################

###  Fit model for the numerator of the stabilized weights ###
# Using original data prior to sequential emulation
psc.num <- speedglm(censor_hu==0 ~ treat_b + cal_time + cal_timesqr,
                    family=binomial(link="logit"),
                    data=med_obs_elig_subset)
summary(psc.num)

# Obtain predicted probabilities of being uncensored in the original data
med_obs_elig$psc.num <- predict(psc.num, med_obs_elig, type = "response")

summary(med_obs_elig$psc.num)

### Estimate stabilized inverse probability weights for censoring ###

# Add probabilities to sequentially emulated data

# Convert data frames to data tables
setDT(med_elig)
setDT(med_obs_elig)

# Merge the data tables
med_elig <- merge(med_elig, med_obs_elig[, c("id", "cal_time", "psc.num", "psc.denom"), with = FALSE],
                  by = c("id", "cal_time"), all.x = TRUE)

# Assign the probability of being uncensored as 1 if a healthcare encounter occurred within the past 12 months
med_elig <- med_elig %>%
  mutate(psc.num = if_else(time_lastHU < HU_time, 1, psc.num),
         psc.denom = if_else(time_lastHU < HU_time, 1, psc.denom))

# Arrange data
med_elig <- med_elig %>% arrange(id_new, time)

# Take cumulative products starting at baseline of a given sequential target trial
med_elig <- med_elig %>%
  group_by(id_new) %>%
  mutate(sw_c = cumprod(psc.num)/cumprod(psc.denom)) %>%
  ungroup() %>%
  mutate(sw_c = ifelse(is.na(psc.num), 1, sw_c))

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: stabilized weights ###

print('Summary of censoring weights:')
summary(med_elig$sw_c)
cat("Standard deviation of sw_c:", sd(med_elig$sw_c), "\n")
quantile(med_elig$sw_c, probs = c(0.975, 0.99, 0.999, 0.9999))
gc()

# Assign censoring weights to be the main weights
med_elig$sw <- med_elig$sw_c

print('Summary of final weights before truncation:')
summary(med_elig$sw)
quantile(med_elig$sw, probs = c(0.975, 0.99, 0.999, 0.9999))

### Truncate final stabilized weight at the 99th percentile (no extremes, so doesn't need truncation) ###
# threshold_99 <- quantile(med_elig$sw, 0.9999)
# med_elig$sw_99 <- med_elig$sw
# med_elig$sw_99[med_elig$sw_99 > threshold_99] <- threshold_99
med_elig$sw_99 <- med_elig$sw

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: truncated weights ###
print('Summary of final weights after truncation:')
summary(med_elig$sw_99)
sd(med_elig$sw_99)

rm(med_obs_elig, trials_all_base, med_obs_elig_subset)
gc()

#################################################################################
#################################################################################

### Fit weighted pooled logistic regression with final stabilized weights ###

if (K %in% c(24,60)) {
  write_feather(med_elig, sprintf("ETT/SSRI/Processed_files/SSRI_med_elig_%d_%s.ft", K, event_name))
}

gc()
Sys.sleep(60)

# Include product terms between time and treatment
fit.pool1 <- speedglm(formula = eventMC==1 ~ treat_b + time + timesqr + Year_cat +
                        I(treat_b*time) +  I(treat_b*timesqr) +
                        as.factor(Age_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) + 
                        as.factor(N.endo_cat_sb) + as.factor(N.GI_cat_sb) + as.factor(N.visits_cat_sb) + 
                        as.factor(N.drugs_cat_sb) + as.factor(CCIW_cat_sb) + as.factor(ind_cat_sb),
                      family = binomial(link = 'logit'),
                      data = med_elig,
                      weights = med_elig$sw_99)

# Recall that the standard variance estimates and CIs reported here are invalid;

# Print results 
summary(fit.pool1)
gc()

### Transform estimates to risks at each time point in each group ###

# Create dataset with all time points for each individual under each treatment level
treat0 <- expandRows(med_elig[which(med_elig$time==0),], count=K, count.is.col=F) 

print('Expanded Rows')

treat0$time <- rep(seq(0, K-1), nrow(med_elig[which(med_elig$time==0),]))
treat0$timesqr <- treat0$time^2
# Under no baseline vaccination
treat0$treat_b <- 0
# Under CROWN vaccination
treat1 <- treat0
treat1$treat_b <- 1

# Extract predicted values from pooled logistic regression model for each person-time row
# Predicted values correspond to discrete-time hazards
treat0$p.event0 <- predict(fit.pool1, treat0, type = "response")
treat1$p.event1 <- predict(fit.pool1, treat1, type = "response")

# Obtain predicted survival probabilities from discrete-time hazards
treat0.surv <- treat0 %>% group_by(id_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
treat1.surv <- treat1 %>% group_by(id_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
treat0.surv$risk0 <- 1 - treat0.surv$surv0
treat1.surv$risk1 <- 1 - treat1.surv$surv1

# Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
risk0 <- aggregate(treat0.surv, by=list(treat0.surv$time), FUN=mean)[c("treat_b", "time", "risk0")]
risk1 <- aggregate(treat1.surv, by=list(treat1.surv$time), FUN=mean)[c("treat_b", "time", "risk1")]

# Prepare data

# Convert data frames to data tables
setDT(risk0)
setDT(risk1)

# Merge the data tables
graph.pred <- merge(risk0, risk1, by = "time")

# Edit data frame to reflect that risks are estimated at the END of each interval
graph.pred$time_0 <- graph.pred$time + 1
zero <- data.frame(cbind(0,0,0,1,0,0))
zero <- setNames(zero,names(graph.pred))
graph <- rbind(zero, graph.pred)

### Use pooled logistic regression estimates to compute causal estimates ###

num_iters <- ifelse(K == 12, 1, ifelse(K == 24, 2, ifelse(K == 60, 3, NA)))

for (i in 1:num_iters) {
  
  time_index <- ifelse(i == 1, 12, ifelse(i == 2, 24, ifelse(i == 3, 60, NA)))
  
  print(paste('Risk estimates for', time_index, 'months', sep = " "))
  
  # risk in no vaccine group
  risk0.plr <- graph$risk0[which(graph$time==time_index-1)]
  cat(risk0.plr, "\n")
  
  # risk in CROWN vaccine group
  risk1.plr <- graph$risk1[which(graph$time==time_index-1)]
  cat(risk1.plr, "\n")
  
  # risk difference
  rd.plr <- risk1.plr - risk0.plr
  cat(rd.plr, "\n")
  
  # risk ratio
  rr.plr <- risk1.plr / risk0.plr
  cat(rr.plr, "\n\n")
}

rm(treat0, treat1)
gc()
Sys.sleep(60)


#################################################################################
#################################################################################

### Construct marginal parametric cumulative incidence (risk) curves ###

if (event_name == "MC") {
  CI_title <- "Microscopic Colitis"
  if (K %in% c(12, 24)) {
    CI_x_by <- 2
    CI_y_max <- 0.004
    CI_y_by <- 0.001
  } else if (K == 60) {
    CI_x_by <- 5
    CI_y_max <- 0.006
    CI_y_by <- 0.0015
  }
} else if (event_name == "NM") {
  CI_title <- "Normal Mucosa"
  if (K %in% c(12, 24)) {
    CI_x_by <- 2
    CI_y_max <- 0.012
    CI_y_by <- 0.003
  } else if (K == 60) {
    CI_x_by <- 5
    CI_y_max <- 0.020
    CI_y_by <- 0.005
  }
}

# Create plot (without CIs)
plot.plr <- ggplot(graph, 
                   aes(x=time_0, y=risk)) + # set x and y axes
  geom_line(aes(y = risk1, # create line for vaccine group
                color = "SSRI"),
            linewidth = 1.5) + 
  geom_line(aes(y = risk0, # create line for no vaccine group
                color = "Mirtazapine"),
            linewidth = 1.5) +
  xlab("Months") + # label x axis
  scale_x_continuous(limits = c(0, K), # format x axis
                     breaks=seq(0, K, by = CI_x_by)) + 
  ylab("Cumulative Incidence (%)") + # label y axis
  scale_y_continuous(limits=c(0, CI_y_max), # format y axis
                     breaks=seq(0, CI_y_max, by = CI_y_by),
                     labels=sprintf("%.2f%%", seq(0, CI_y_max, by = CI_y_by) * 100)) + 
  labs(title = CI_title) + # add title
  theme_minimal()+ # set plot theme elements
  theme(axis.text = element_text(size=14), legend.position = c(0.2, 0.8),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  font("xlab",size=14)+
  font("ylab",size=14)+
  font("legend.text",size=10)+
  scale_color_manual(values=c("#E7B800","#2E9FDF"), # set colors
                     breaks=c('Mirtazapine', 'SSRI'))  

# Increase the white margin around the image
plot.plr <- plot.plr + theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))

# Plot
plot.plr

ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave(paste0("ETT/SSRI/Output/CI_K", K, "_HU", HU_time, "_", event_name, ".png"), plot = plot.plr, width = 8, height = 8)
while (!is.null(dev.list()))  dev.off()

#################################################################################
#################################################################################

### Construct a covariate balance plot for the weighted population ###

# Create subsets of data, according to treat_b status
treat_b0 <- subset(med_elig,treat_b==0)
treat_b1 <- subset(med_elig,treat_b==1)

dim(treat_b0)
dim(treat_b1)

rm(med_elig)
gc()

# List of variables to compare
varlist <- c("Age_65-70","Age_71-75","Age_76-80", "Age_81-85", "Age_86+",
             "Sex", "nordic",
             "Year_2006-2008", "Year_2009-2011", "Year_2012-2014", "Year_2015-2017",
             "educ0", "educ1", "educ2", 
             "N.GI_sb_1+", "N.endo_sb_1+",
             "N.visits_sb_0","N.visits_sb_1-5","N.visits_sb_6-10", "N.visits_sb_11+",
             "N.drugs_sb_0","N.drugs_sb_1-4", "N.drugs_sb_5-8", "N.drugs_sb_9-12", "N.drugs_sb_13-16", "N.drugs_sb_17+",
             "CCI_sb_0","CCI_sb_1-2","CCI_sb_3-4","CCI_sb_5+",
             "ind_SSRI_sb_0", "ind_SSRI_sb_1", "ind_SSRI_sb_2", "ind_SSRI_sb_3", "ind_SSRI_sb_4", "ind_SSRI_sb_5",
             "ind_SSRI_sb_6", "ind_SSRI_sb_7", "ind_SSRI_sb_8", "ind_SSRI_sb_9+")

# Convert treat_b0 and treat_b1 to data.tables
setDT(treat_b0)
setDT(treat_b1)

# Create function to take mean difference, or SMD for age
meanfctn <- function(x){
  if (x %in% c("Age")){
    t0 <- treat_b0[[x]]
    t1 <- treat_b1[[x]]
    md <- (mean(t1) - mean(t0))/sd(t1)
  }else{
    t0 <- treat_b0[[x]]
    t1 <- treat_b1[[x]]
    md <- mean(t1) - mean(t0)}
  return(c(var = x, md = md))
}

# Calculate mean differences for covariates (SMD for age)
wmean_fctn <- function(x){
  if(x %in% c("Age")){
    
    print(x)
    print(length(treat_b1[[x]]))
    print(length(treat_b1$sw_99))
    print(sum(is.na(treat_b1[[x]])))
    print(sum(is.na(treat_b1$sw_99)))
    
    md <- (weighted.mean(treat_b1[[x]], treat_b1$sw_99) - weighted.mean(treat_b0[[x]], treat_b0$sw_99))/sqrt(wtd.var(treat_b1[[x]],                                                                                                               treat_b1$w_a))
  }else{
    t0 <- weighted.mean(treat_b0[[x]], treat_b0$sw_99)
    t1 <- weighted.mean(treat_b1[[x]], treat_b1$sw_99)
    md <- t1-t0}
  return(c(t0 = round(t0*100, 3), t1 = round(t1*100, 3), var = x, md = round(md, 4)))
  
}

# Function to calculate weighted statistics for Age
weighted_stats <- function(data, variable, weight) {
  wtd_mean <- wtd.mean(data[[variable]], data[[weight]])
  wtd_sd <- sqrt(wtd.var(data[[variable]], data[[weight]], normwt=TRUE))
  wtd_median <- wtd.quantile(data[[variable]], data[[weight]], probs = 0.5, normwt=TRUE)
  wtd_q25 <- wtd.quantile(data[[variable]], data[[weight]], probs = 0.25, normwt=TRUE)
  wtd_q75 <- wtd.quantile(data[[variable]], data[[weight]], probs = 0.75, normwt=TRUE)
  
  return(data.frame(
    Mean = wtd_mean,
    SD = wtd_sd,
    Median = wtd_median,
    Q25 = wtd_q25,
    Q75 = wtd_q75
  ))
}

print('Computing weighted means')

# Calculate weighted statistics for treat_b0
stats_treat_b0 <- weighted_stats(treat_b0, "Age", "sw_99")
print("Weighted Statistics for treat_b0:")
print(stats_treat_b0)

# Calculate weighted statistics for treat_b1
stats_treat_b1 <- weighted_stats(treat_b1, "Age", "sw_99")
print("Weighted Statistics for treat_b1:")
print(stats_treat_b1)
print('\n\n')

# Create the covariate plot
covplot_w <- lapply(varlist, wmean_fctn) %>% do.call(rbind,.) %>% as.data.frame()
covplot_w$md <- as.numeric(covplot_w$md)

# Order the dataframe by the 'var' column
covplot_w$var <- factor(covplot_w$var, levels = varlist)
covplot_w <- covplot_w[order(covplot_w$var), ]
covplot_w$var <- gsub("_sb", "", covplot_w$var)

print(covplot_w)

# Plot it

covplot_weighted <- ggplot(data = covplot_w) +
  geom_point(aes(x = md, y = var), color = "steelblue") + scale_x_continuous(limits = c(-0.25, 0.25)) +
  geom_vline(xintercept = 0) +
  labs(y = "Covariates", x = "Mean Difference", title = "Covariate Balance Plot")

# Increase the white margin around the image
covplot_weighted <- covplot_weighted + theme(plot.margin = margin(0.9,0.9,0.9,0.9, "cm"))

ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave(paste0("ETT/SSRI/Output/covplot_K", K, "_HU", HU_time, "_", event_name, ".png"), plot = covplot_weighted, width = 8, height = 9)
while (!is.null(dev.list()))  dev.off()

