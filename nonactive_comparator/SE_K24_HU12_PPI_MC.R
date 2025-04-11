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
if (!require("disk.frame")) install.packages("disk.frame")
library(disk.frame)
if (!require("scales")) install.packages("scales")
library(scales)
if (!require("future")) install.packages("future")
library(future)
if (!require("furrr")) install.packages("furrr")
library(furrr)
if (!require("rlang")) install.packages("rlang")
library(rlang)

setwd('~/Medications/')

# Set healthcare utililization (HU) censoring time, event name, and time point for estimation of risks in months
HU_time <- 12
event_name <- 'MC'
K <- 24

# Import long format data file with covariate and outcome data
med_obs <- read_feather(sprintf("ETT/PPI/Processed_files/PPI_longdata_GI_%s.ft", event_name))

med_obs$id <- med_obs$Lopnr
length(unique(med_obs$Lopnr))

# Set eligibility and encounter dates based on censoring time for HU
med_obs$elig <- med_obs$elig_gi12
med_obs$encounter_date <- med_obs$encounter_date12
med_obs$encounter_e_date <- med_obs$encounter_e_date12

med_obs <- med_obs %>%
  dplyr::select(1:10, elig, 11:ncol(med_obs))
med_obs <- med_obs %>%
  dplyr::select(1:3, encounter_date, 4:ncol(med_obs))

# Create a time varying continuous variable for indication
med_obs$ind_PPI_con <- as.numeric(difftime(med_obs$current_date, med_obs$indication_date, units = "days")+1)/365.25

# Replace NAs with 0 for time since indication
med_obs$ind_PPI_con[is.na(med_obs$ind_PPI_con)] <- 0

# Create a replicate of encounter date but with no NAs. This will be used for computing time since last HU
# Fill NAs in encounter_date column with the last available non-NA value
med_obs <- med_obs %>%
  group_by(Lopnr) %>%
  mutate(enc_date_nomiss = zoo::na.locf(encounter_date))

# Calculate time difference in months between current_date and encounter_date
med_obs$time_lastHU <- round(as.numeric(interval(med_obs$enc_date_nomiss, med_obs$current_date) / months(1)), 2)
med_obs$time_lastHU_days <- round(as.numeric(interval(med_obs$enc_date_nomiss, med_obs$current_date) / days(1)), 2)

# Set column order
med_obs <- med_obs %>%
  dplyr::select(1:3, time_lastHU, 4:ncol(med_obs))

med_obs <- med_obs %>%
  dplyr::select(1:4, time_lastHU_days, 5:ncol(med_obs))

### Preprocess variables ###

### Calendar Year

med_obs <- med_obs %>%
  mutate(date_admission = as.Date(date_admission),
         Year = as.numeric(format(date_admission, "%Y")))

# Define the breaks for the categories and create a factor variable
breaks <- c(2005, 2008, 2011, 2014, 2017)
med_obs$Year_cat <- cut(med_obs$Year, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$Year_cat <- as.factor(med_obs$Year_cat)

### Age

# Define the breaks for the categories and create a factor variable
breaks <- c(64, 70, 75, 80, 85, Inf)
med_obs$Age_cat <- cut(med_obs$Age, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$Age_cat <- as.factor(med_obs$Age_cat)

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

### Med indication

# Define the breaks for the categories and create a factor variable
breaks_ind <- c(-1, 2, 5, Inf)
med_obs$ind_cat <- cut(med_obs$ind_PPI_con, breaks = breaks_ind, include.lowest = FALSE, right = TRUE)
med_obs$ind_cat <- as.factor(med_obs$ind_cat)

### CCI

# Define the breaks for the categories and create a factor variable
breaks <- c(-1, 2, 4, Inf)
med_obs$CCIW_cat <- cut(med_obs$CCIW, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$CCIW_cat <- as.factor(med_obs$CCIW_cat)

# Baseline
med_obs$CCIW_cat_sb <- cut(med_obs$CCIW_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$CCIW_cat_sb <- as.factor(med_obs$CCIW_cat_sb)

### N endo
breaks <- c(-1, 0, Inf)
med_obs$N.endo_cat <- cut(med_obs$N.endo, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.endo_cat <- as.factor(med_obs$N.endo_cat)

# Baseline
med_obs$N.endo_cat_sb <- cut(med_obs$N.endo_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.endo_cat_sb <- as.factor(med_obs$N.endo_cat_sb)

### N GI
breaks <- c(-1, 0, Inf)
med_obs$N.GI_cat <- cut(med_obs$N.GI, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.GI_cat <- as.factor(med_obs$N.GI_cat)

# Baseline
med_obs$N.GI_cat_sb <- cut(med_obs$N.GI_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.GI_cat_sb <- as.factor(med_obs$N.GI_cat_sb)

### N visits
breaks <- c(-1, 5, 10, Inf)
med_obs$N.visits_cat <- cut(med_obs$N.visits, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.visits_cat <- as.factor(med_obs$N.visits_cat)

# Baseline
med_obs$N.visits_cat_sb <- cut(med_obs$N.visits_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.visits_cat_sb <- as.factor(med_obs$N.visits_cat_sb)

### N drugs
breaks <- c(-1, 1, 4, 7, 10, 13, 16, 19, Inf)
med_obs$N.drugs_cat <- cut(med_obs$N.drugs, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.drugs_cat <- as.factor(med_obs$N.drugs_cat)

# Basline
med_obs$N.drugs_cat_sb <- cut(med_obs$N.drugs_sb, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_obs$N.drugs_cat_sb <- as.factor(med_obs$N.drugs_cat_sb)

### Create a variable for cumulative of treatment variable
med_obs <- med_obs %>%
  group_by(Lopnr) %>%
  mutate(treat_cum = cumsum(treat)) %>%
  ungroup()

### Create healthcare utilization (HU) censoring variable
med_obs <- med_obs %>%
  group_by(id) %>%
  mutate(censor_hu = ifelse(is.na(encounter_date), 1, 0)) %>%
  ungroup()

### Update the eligibility variable

# Update 'elig' to 0 where 'censor_hu' is equal to 1
med_obs$elig[med_obs$censor_hu == 1] <- 0

# Set elig to 0 for rows where encounter_date is not between current_date and current_date_end
med_obs$elig <- ifelse(!is.na(med_obs$encounter_e_date) &
                         (med_obs$encounter_e_date >= med_obs$current_date & med_obs$encounter_e_date <= med_obs$current_date_end),
                       med_obs$elig, 0)

# Set elig to 0 when no visits or prescriptions are present in the past 12 months
med_obs$elig <- ifelse(med_obs$N.visits == 0 & med_obs$N.drugs == 0, 0, med_obs$elig)

### Update censoring due to loss of follow-up variable
med_obs <- med_obs %>%
  group_by(id) %>%
  mutate(censor_enc = ifelse(N.visits == 0 & N.drugs == 0, 1, 0)) %>%
  ungroup()

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

dim(trials_all)

#################################################################################
#################################################################################

##### Create censoring variables ######

# Censoring due to non-adherence to baseline treatment strategy
trials_all <- trials_all %>%
  group_by(id_new) %>%
  mutate(censor_trt = if_else(treat == 1 & lag(treat, default = 0) == 0, 1, 0)) %>%
  mutate(censor_trt = if_else(any(time == 0 & treat == 1) & treat == 1, 0, censor_trt)) %>%
  mutate(censor_trt = if_else(cumsum(censor_trt)>=1, 1, 0)) %>%
  ungroup()

trials_all <- trials_all %>%
  dplyr::select(1:16, censor_trt, 17:ncol(trials_all))

# Censoring due to no healthcare utilization in the past 12 months
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


#### Summary statistics ####

# Count the number of unique ids where censor_hu is equal to 1
unique_hu_censor <- trials_all %>%
  filter(censor_hu == 1) %>%
  summarise(unique_ids = n_distinct(id))


# Print the result
cat('Censoring due to no health care utilization:', unique_hu_censor$unique_ids, "\n")
rm(unique_hu_censor)

# Count the number of unique ids where censor_trt is equal to 1
unique_trt_censor <- trials_all %>%
  filter(censor_trt == 1) %>%
  summarise(unique_ids = n_distinct(id))

# Print the result
cat('Censoring due to non-adherence to treatment strategy:', unique_trt_censor$unique_ids, "\n")
rm(unique_trt_censor)

# Count the number of unique ids where IBD equal to 1 in sequential emulation
unique_ibd <- trials_all %>%
  filter(current_date > First_IBD_DATE) %>%
  summarise(unique_ids = n_distinct(id))

# Print the result
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

### Describe number of (i) unique and (ii) non-unique individuals lost to follow-up ###

# Count & proportion of unique individuals censored
count_censor_hu <- trials_all[which(trials_all$censor_hu_any==1),]
no_dups_censor_hu <- count_censor_hu[!duplicated(count_censor_hu$id),]
no_dups_censor_hu %>% dplyr::summarise(
  n_censored_hu = length(censor_hu_any[censor_hu_any == 1]),
  prop_censored_hu = length(censor_hu_any[censor_hu_any == 1])/length(trials_all[!duplicated(trials_all$id),]$id))

# Count & proportion of non-unique individuals censored
trials_all %>% filter(time==0) %>% dplyr::summarise(
  n_censored_hu = length(censor_hu_any[censor_hu_any == 1]),
  prop_censored_hu = mean(censor_hu_any))

rm(count_censor_hu)

#################################################################################
#################################################################################

### Exclude censored rows

trials_all <- trials_all %>% filter(!(censor_trt==1 | censor_enc==1 | censor_hu==1 | death==1))

trials_all <- trials_all %>%
  group_by(id_new) %>%
  mutate(month_num = row_number() - 1) %>%
  ungroup()

### Create dataset restricted to baseline (follow-up time = 0)
trials_all_base <- trials_all %>%
  dplyr::filter(time == 0)

# Create indicator for baseline treatment group
trials_all_base$treat_b <- trials_all_base$treat

# Create period variables for baseline week of each sequential target trial
trials_all_base$period <- trials_all_base$cal_time
trials_all_base$periodsqr <- trials_all_base$cal_timesqr

# Create confounder variables for baseline week of each sequential target trial
trials_all_base$CCIW_b <- trials_all_base$CCIW
trials_all_base$N.endo_b <- trials_all_base$N.endo
trials_all_base$N.GI_b <- trials_all_base$N.GI
trials_all_base$N.visits_b <- trials_all_base$N.visits
trials_all_base$N.drugs_b <- trials_all_base$N.drugs
trials_all_base$ind_PPI_con_b <- trials_all_base$ind_PPI_con

# Add baseline treatment group to person-time data

# Convert data.frames to data.tables
setDT(trials_all)
setDT(trials_all_base)

# Set key for fast merging
setkey(trials_all, id, trial_num)
setkey(trials_all_base, id, trial_num)

# Perform the merge and select required columns
med_elig <- merge(trials_all, trials_all_base[,c("id","trial_num","treat_b",
                                                 "period","periodsqr",
                                                 "CCIW_b", "N.endo_b", "N.GI_b", "N.visits_b", "N.drugs_b", "ind_PPI_con_b")],
                  by=c("id","trial_num"))

rm(trials_all)
gc()

# Create indicator of any MC occurrence over the follow-up period
med_elig <- med_elig %>%
  dplyr::group_by(id_new) %>%
  dplyr::mutate(
    MC_any = ifelse(is.na(eventMC), NA, max(eventMC, na.rm = T))) %>%
  dplyr::ungroup()

print('Any MC occurence over the follow-up period')
table(med_elig$MC_any[which(med_elig$time==0)], useNA = "always")

### Convert continuous confounders to factors ###

### Med indication

# Define the breaks for the categories and create a factor variable
breaks_ind <- c(-1, 2, 5, Inf)
med_elig$ind_cat_b <- cut(med_elig$ind_PPI_con_b, breaks = breaks_ind, include.lowest = FALSE, right = TRUE)
med_elig$ind_cat_b <- as.factor(med_elig$ind_cat_b)

### CCI

# Define the breaks for the categories and create a factor variable
breaks <- c(-1, 2, 4, Inf)

# Baseline
med_elig$CCIW_cat_b <- cut(med_elig$CCIW_b, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_elig$CCIW_cat_b <- as.factor(med_elig$CCIW_cat_b)

### N endo
breaks <- c(-1, 0, Inf)

# Baseline
med_elig$N.endo_cat_b <- cut(med_elig$N.endo_b, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_elig$N.endo_cat_b <- as.factor(med_elig$N.endo_cat_b)

### N GI
breaks <- c(-1, 0, Inf)

# Baseline
med_elig$N.GI_cat_b <- cut(med_elig$N.GI_b, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_elig$N.GI_cat_b <- as.factor(med_elig$N.GI_cat_b)

### N visits
breaks <-  c(-1, 5, 10, Inf)

# Baseline
med_elig$N.visits_cat_b <- cut(med_elig$N.visits_b, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_elig$N.visits_cat_b <- as.factor(med_elig$N.visits_cat_b)

### N drugs
breaks <- c(-1, 1, 4, 7, 10, 13, 16, 19, Inf)

# Basline
med_elig$N.drugs_cat_b <- cut(med_elig$N.drugs_b, breaks = breaks, include.lowest = FALSE, right = TRUE)
med_elig$N.drugs_cat_b <- as.factor(med_elig$N.drugs_cat_b)

### Create a descriptive table ###

# Define categorical columns

selected_columns_obs <- c("Lopnr", "dataset", "time_lastHU", "current_date", "treat", 
                          "eventMC", "timeMC", "cal_time", "cal_timesqr", "Year", "Year_cat",
                          "Age", "Sex", "ind_PPI", "ind_PPI_sb", "N.GI", "N.endo", "N.drugs", 
                          "N.visits", "N.GI_sb", "N.endo_sb", "N.drugs_sb", "N.visits_sb", 
                          "educ", "nordic", "CCIW", "CCIW_sb", "ID", "id", "Age_cat", 
                          "ind_cat", "CCIW_cat", "CCIW_cat_sb", "N.endo_cat", "N.endo_cat_sb", 
                          "N.GI_cat", "N.GI_cat_sb", "N.visits_cat", "N.visits_cat_sb", 
                          "N.drugs_cat", "N.drugs_cat_sb", "treat_cum", "censor_hu", "censor_enc")

categorical_columns_obs <- c("censor_hu",
                             "Age_cat",  "Year", "Year_cat", "Sex", "educ", "nordic",
                             "N.endo_cat", "N.GI_cat", "N.visits_cat",
                             "N.endo_cat_sb", "N.GI_cat_sb", "N.visits_cat_sb",
                             "N.drugs_cat_sb", "CCIW_cat_sb", "ind_PPI_sb")

selected_columns_elig <- c("id", "trial_num", "Lopnr", "dataset", "time_lastHU", "id_new", 
                           "censor_hu", "treat", "censor", "eventMC", "timeMC",  "Year", "Year_cat",
                           "cal_time", "cal_timesqr", "Age", "Sex", "ind_PPI", "ind_PPI_sb", "N.GI", 
                           "N.endo", "N.drugs", "N.visits", "N.GI_sb", "N.endo_sb", "N.drugs_sb", 
                           "N.visits_sb", "educ", "nordic", "CCIW", "CCIW_sb", "ID", "ind_PPI_con", 
                           "Age_cat", "ind_cat", "CCIW_cat", "CCIW_cat_sb", "N.endo_cat", "N.endo_cat_sb", 
                           "N.GI_cat", "N.GI_cat_sb", "N.visits_cat", "N.visits_cat_sb", "N.drugs_cat", 
                           "N.drugs_cat_sb", "treat_cum", "time", "timesqr", "treat_b", "period", 
                           "periodsqr", "CCIW_b", "N.endo_b", "N.GI_b", "N.visits_b", "N.drugs_b", 
                           "MC_any", "ind_cat_b", "CCIW_cat_b", "N.endo_cat_b", "N.GI_cat_b", 
                           "N.visits_cat_b", "N.drugs_cat_b")

categorical_columns_elig <- c("censor_hu",
                              "Age_cat", "Year", "Year_cat", "Sex", "educ", "nordic",
                              "N.GI_cat", "N.endo_cat", "N.visits_cat",
                              "N.drugs_cat", "CCIW_cat", "ind_cat",
                              "N.endo_cat_b", "N.GI_cat_b", "N.visits_cat_b",
                              "N.drugs_cat_b", "CCIW_cat_b", "ind_cat_b",
                              "N.endo_cat_sb", "N.GI_cat_sb", "N.visits_cat_sb",
                              "N.drugs_cat_sb", "CCIW_cat_sb", "ind_PPI_sb")

# Initialize disk.frame
setup_disk.frame()
options(future.globals.maxSize = Inf)

disk_directory <- '/tmp/'

# Initialize final disk.frames
med_obs_disk <- as.disk.frame(med_obs, outdir = file.path(disk_directory, paste0("med_obs.df")), overwrite = TRUE)
med_elig_disk <- as.disk.frame(med_elig, outdir = file.path(disk_directory, paste0("med_elig.df")), overwrite = TRUE)

rm(med_obs, med_elig)
gc()


#############################################

# Get number of rows and columns for med_obs_disk
nrows_obs <- nrow(med_obs_disk)
ncols_obs <- ncol(med_obs_disk)

# Get number of rows and columns for med_elig_disk
nrows_elig <- nrow(med_elig_disk)
ncols_elig <- ncol(med_elig_disk)

# Print the dimensions
cat(sprintf("Final med_obs_disk has %d rows and %d columns\n", nrows_obs, ncols_obs))
cat(sprintf("Final med_elig_disk has %d rows and %d columns\n", nrows_elig, ncols_elig))

# Check disk usage after performing H2O operations:
print('Memory and swap usage:')
system("du -sh /tmp")
system("df -h /tmp", intern = TRUE)
system("awk '/MemTotal/ {total=$2} /MemAvailable/ {available=$2} END {used=total-available; print \"Total: \" total/1024/1024 \"GB, Used: \" used/1024/1024 \"GB (\" used/total*100 \"%), Available: \" available/1024/1024 \"GB (\" available/total*100 \"%)\"}' /proc/meminfo")

### Create a descriptive table ###

new_path <- tempfile(fileext = ".df")
med_elig_disk <- rechunk(med_elig_disk, 1, new_path, overwrite = FALSE)


# Number of unique IDs
unique_ids_summary <- med_elig_disk %>%
  filter(time == 0) %>%
  summarise(unique_ids = n_distinct(id)) %>%
  collect()

# Display the summary
print("Unique individuals in the cohort:")
print(unique_ids_summary)

# Number of unique IDs at time == 0 grouped by treat_b
unique_new_ids_n <- med_elig_disk %>%
  filter(time == 0) %>%
  summarise(unique_id_new = n(id_new)) %>%
  collect()

# Display the summary
print("Non-unique individuals in the cohort:")
print(unique_new_ids_n)

# Number of unique IDs at time == 0 grouped by treat_b
unique_new_ids_n_grp <- med_elig_disk %>%
  filter(time == 0) %>%
  group_by(treat_b) %>%
  summarise(unique_id_new = n(id_new)) %>%
  collect()

print('Non-unique individuals at baseline by treat_b')
print(unique_new_ids_n_grp)

# Function to summarize unique id_new and ID occurrences for a given time filter
summarize_eventMC <- function(df, time_filter) {
  id_new_summary <- df %>%
    filter(time < time_filter, eventMC == 1) %>%
    chunk_group_by(id_new) %>%
    chunk_summarise(count = n()) %>%
    chunk_summarise(unique_ids = n()) %>%
    collect()
  
  id_summary <- df %>%
    filter(time < time_filter, eventMC == 1) %>%
    chunk_group_by(ID) %>%
    chunk_summarise(count = n()) %>%
    chunk_summarise(unique_ids = n()) %>%
    collect()
  
  list(id_new_summary = id_new_summary, id_summary = id_summary)
}

# Summarize for different time filters
eventMC_summary_60 <- summarize_eventMC(med_elig_disk, 60)
eventMC_summary_24 <- summarize_eventMC(med_elig_disk, 24)
eventMC_summary_12 <- summarize_eventMC(med_elig_disk, 12)

# Print the results with appropriate titles
cat("Number of unique id_new occurrences for eventMC with time < 60:\n")
print(eventMC_summary_60$id_new_summary)
cat("Number of unique ID occurrences for eventMC with time < 60:\n")
print(eventMC_summary_60$id_summary)

cat("\nNumber of unique id_new occurrences for eventMC with time < 24:\n")
print(eventMC_summary_24$id_new_summary)
cat("Number of unique ID occurrences for eventMC with time < 24:\n")
print(eventMC_summary_24$id_summary)

cat("\nNumber of unique id_new occurrences for eventMC with time < 12:\n")
print(eventMC_summary_12$id_new_summary)
cat("Number of unique ID occurrences for eventMC with time < 12:\n")
print(eventMC_summary_12$id_summary)

## Summary statistics for continuous variable Age
age_summary <- med_elig_disk %>%
  filter(time == 0) %>%
  chunk_group_by(treat_b) %>%
  chunk_summarise(
    mean_age = mean(Age, na.rm = TRUE),
    median_age = median(Age, na.rm = TRUE),
    sd_age = sd(Age, na.rm = TRUE),
    iqr_age = IQR(Age, na.rm = TRUE),
    p25_age = quantile(Age, 0.25, na.rm = TRUE),
    p75_age = quantile(Age, 0.75, na.rm = TRUE)
  ) %>%
  collect()

print(age_summary)

# List of categorical variables
catVars <- c("Age_cat", "Year_cat", "nordic", "Sex", "educ", "N.GI_cat", "N.endo_cat", 
             "N.visits_cat", "N.drugs_cat", "CCIW_cat", "ind_cat", "Year")

# Define the function to calculate sums and percentages for categorical variables
summarize_categorical <- function(df, catVars) {
  results <- list()
  
  for (var in catVars) {
    summary_df <- df %>%
      filter(time == 0) %>%
      chunk_group_by(treat_b, !!sym(var)) %>%
      chunk_summarise(
        count = n()
      ) %>%
      collect() %>%
      group_by(treat_b) %>%
      mutate(percent = 100 * count / sum(count, na.rm = TRUE)) %>%
      arrange(treat_b, desc(count))
    
    results[[var]] <- summary_df
  }
  
  return(results)
}

## Call the function and store the results
summary_results <- summarize_categorical(med_elig_disk, catVars)

print(summary_results)

#################################################################################
#################################################################################

# Create a baseline dataset without duplicate IDs
no_dups <- med_elig_disk %>%
  filter(!duplicated(id))

# Number of unique individuals
cat('\nN, unique individuals:\n')
unique_individuals <- no_dups %>%
  summarise(n = n_distinct(id)) %>%
  collect()
cat(unique_individuals$n)

# Number of unique events
cat('\nN, unique events:\n')
unique_events <- med_elig_disk %>%
  filter(MC_any == 1) %>%
  filter(!duplicated(id)) %>%
  summarise(n = n_distinct(id)) %>%
  collect()
cat(unique_events$n)

# Number of non-unique individuals
cat('\nN, non-unique individuals:\n')
non_unique_individuals <- med_elig_disk %>%
  summarise(n = n()) %>%
  collect()
cat(non_unique_individuals$n)

# Number of non-unique events
cat('\nN, non-unique events:\n')
non_unique_events <- med_elig_disk %>%
  filter(time == 0) %>%
  group_by(MC_any) %>%
  summarise(n = n()) %>%
  collect()
print(table(non_unique_events$MC_any, non_unique_events$n))

cat('N, non-unique events, treat==0:\n')
non_unique_events_treat_0 <- med_elig_disk %>%
  filter(time == 0 & treat_b == 0) %>%
  group_by(MC_any) %>%
  summarise(n = n()) %>%
  collect()
print(table(non_unique_events_treat_0$MC_any, non_unique_events_treat_0$n))

cat('N, non-unique events, treat==1:\n')
non_unique_events_treat_1 <- med_elig_disk %>%
  filter(time == 0 & treat_b == 1) %>%
  group_by(MC_any) %>%
  summarise(n = n()) %>%
  collect()
print(table(non_unique_events_treat_1$MC_any, non_unique_events_treat_1$n))

# Compute average number of trials each individual contributed to
avg_trials <- med_elig_disk %>%
  group_by(id) %>%
  summarise(trials_count = n()) %>%
  collect()

# Print average number of trials
cat('\nAverage number of trials each individual contributed to:\n')
print(summary(avg_trials$trials_count))

#################################################################################
#################################################################################

# Extract unique IDs from med_elig_disk
unique_ids_df <- med_elig_disk %>%
  dplyr::select(id) %>%
  chunk_distinct() %>%
  collect()

# Filter med_obs_disk for these unique IDs
med_obs_disk_elig <- med_obs_disk %>%
  inner_join(unique_ids_df, by = "id", merge_by_chunk_id = FALSE) %>%
  mutate(treat_1 = ifelse(treat_cum >= 1, 1, 0))

# Filter the dataset for a specific condition
med_obs_disk_elig_hu <- med_obs_disk_elig %>%
  filter(time_lastHU >= HU_time & time_lastHU < HU_time + 1)

# Ensure censor_hu is numeric within disk.frame operations
med_obs_disk_elig_hu <- med_obs_disk_elig_hu %>%
  mutate(censor_hu = as.numeric(as.character(censor_hu)))

# Calculate the mean of censor_hu within disk.frame
mean_censor_hu <- med_obs_disk_elig_hu %>%
  summarise(mean_censor_hu = mean(censor_hu, na.rm = TRUE)) %>%
  collect() %>%
  pull(mean_censor_hu) * 100
cat("Mean of censor_hu:", mean_censor_hu, "\n")

# Fit the pooled logistic model using speedglm
psc.denom <- speedglm(censor_hu == 0 ~ treat_1 + cal_time + cal_timesqr + Year_cat +
                        as.factor(Age_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) +
                        as.factor(N.endo_cat_sb) + as.factor(N.visits_cat_sb) +
                        as.factor(N.drugs_cat_sb) + as.factor(CCIW_cat_sb) + as.factor(ind_PPI_sb) +
                        as.factor(CCIW_cat),
                      family = binomial(link = "logit"),
                      data = med_obs_disk_elig_hu)

# Print model summary
summary(psc.denom)
gc()

# Obtain predicted probabilities of being uncensored in the original data
med_obs_disk_elig$psc.denom <- as.vector(predict(psc.denom, med_obs_disk_elig, type = "response"))
med_obs_disk_elig$psc.denom <- as.numeric(med_obs_disk_elig$psc.denom)

cat("\nSummary of predicted probabilities (psc.denom):\n")
summary(med_obs_disk_elig$psc.denom)
gc()

#################################################################################
#################################################################################

###  Fit model for the numerator of the stabilized weights ###
# Using original data prior to sequential emulation
psc.num <- speedglm(censor_hu==0 ~ treat_1 + cal_time  + cal_timesqr + Year_cat +
                      as.factor(Age_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) +
                      as.factor(N.endo_cat_sb) + as.factor(N.visits_cat_sb) +
                      as.factor(N.drugs_cat_sb) + as.factor(CCIW_cat_sb) + as.factor(ind_PPI_sb),
                    family=binomial(link="logit"),
                    data=med_obs_disk_elig_hu)
summary(psc.num)

# Obtain predicted probabilities of being uncensored in the original data
med_obs_disk_elig$psc.num <- as.vector(predict(psc.num, med_obs_disk_elig, type = "response"))
med_obs_disk_elig$psc.num <- as.numeric(med_obs_disk_elig$psc.num)

summary(med_obs_disk_elig$psc.num)
gc()

### Estimate stabilized inverse probability weights for censoring ###

# Select the desired columns using srckeep

# Step 1: Select and collect id and cal_time
df1 <- med_obs_disk_elig %>%
  dplyr::select(id, cal_time) %>%
  collect()
df1 <- as.data.table(df1)

# Step 2: Extract psc.num and psc.denom using $
psc_num <- med_obs_disk_elig$psc.num
psc_denom <- med_obs_disk_elig$psc.denom

# Step 3: Combine into a single DataFrame
med_obs_disk_elig_sub <- df1 %>%
  mutate(psc.num = psc_num, psc.denom = psc_denom)
med_obs_disk_elig_sub <- as.disk.frame(med_obs_disk_elig_sub)

new_path <- tempfile(fileext = ".df")
med_obs_disk_elig_sub <- rechunk(med_obs_disk_elig_sub, 1, new_path, overwrite = TRUE)

rm(df1)
gc()

# Merge estimated probabilities with the primary dataset
med_elig_disk <- merge(med_elig_disk, med_obs_disk_elig_sub,
                       by=c("id","cal_time"), all.x = TRUE, merge_by_chunk_id = TRUE)

rm(med_obs_disk_elig)
gc()

# Assign the probability of being uncensored as 1 if a healthcare encounter occurred within the past 12 months
med_elig_disk <- med_elig_disk %>%
  mutate(psc.num = if_else(time_lastHU < HU_time, 1, psc.num),
         psc.denom = if_else(time_lastHU < HU_time, 1, psc.denom))

# Take cumulative products starting at baseline of a given sequential target trial
med_elig_disk <- med_elig_disk %>%
  chunk_group_by(id, trial_num) %>%
  mutate(cumprod_psc_num = cumprod(psc.num),
         cumprod_psc_denom = cumprod(psc.denom),
         sw_c = cumprod_psc_num/cumprod_psc_denom
  ) %>%
  chunk_ungroup() 

med_elig_disk <- med_elig_disk %>%
  dplyr::mutate(sw_c = ifelse(is.na(psc.num), 1, sw_c))

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: stabilized weights ###

# Summary of censoring weights
summary(med_elig_disk %>% dplyr::select(sw_c) %>% collect())
cat("Standard deviation of sw_c:", sd(med_elig_disk %>% dplyr::select(sw_c) %>% collect() %>% pull(sw_c)), "\n")

med_elig_disk %>%
  dplyr::select(sw_c) %>%
  collect() %>%
  { quantile(.$sw_c, probs = c(0.90, 0.95, 0.975, 0.99, 0.999, 0.9999), na.rm = TRUE) }
gc()

# Check disk usage after performing H2O operations:
print('Memory and swap usage:')
system("du -sh /tmp")
system("df -h /tmp", intern = TRUE)
system("awk '/MemTotal/ {total=$2} /MemAvailable/ {available=$2} END {used=total-available; print \"Total: \" total/1024/1024 \"GB, Used: \" used/1024/1024 \"GB (\" used/total*100 \"%), Available: \" available/1024/1024 \"GB (\" available/total*100 \"%)\"}' /proc/meminfo")

#################################################################################
#################################################################################

### Fit a model to estimate the denominator of the stabilized weights for confounding ###

ipw.denom <- speedglm(treat ~ cal_time + cal_timesqr + as.factor(Year_cat) +
                        as.factor(Age_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) +
                        as.factor(N.endo_cat) + as.factor(N.visits_cat) +
                        as.factor(N.drugs_cat) + as.factor(CCIW_cat) + as.factor(ind_cat) +
                        as.factor(N.endo_cat_b) + as.factor(N.visits_cat_b) +
                        as.factor(N.drugs_cat_b) + as.factor(CCIW_cat_b) + as.factor(ind_cat_b),
                      family=binomial(link="logit"),
                      data=med_elig_disk %>% filter(time > 0), sparse = FALSE)

summary(ipw.denom)

### Fit a model to estimate the numerator of the stabilized weights for confounding ###
ipw.num <- speedglm(treat ~ cal_time + cal_timesqr +
                      as.factor(Sex) + as.factor(educ) + as.factor(nordic) +
                      as.factor(N.endo_cat_b) + as.factor(N.visits_cat_b) +
                      as.factor(N.drugs_cat_b) + as.factor(CCIW_cat_b) + as.factor(ind_cat_b),
                    family=binomial(link="logit"),
                    data=med_elig_disk %>% filter(time > 0), sparse=FALSE)
summary(ipw.num)

#  Estimate denominator probabilities
ipw_denom_v <- as.vector(predict(ipw.denom, med_elig_disk, type = "response"))
ipw_num_v <- as.vector(predict(ipw.num, med_elig_disk, type = "response"))

# Remove variables no longer needed
rm(ipw.num, ipw.denom, psc.num, psc.denom)
rm(psc_num, psc_denom)
delete(med_obs_disk)
delete(med_obs_disk_elig_sub)
gc()
Sys.sleep(10)

# Convert weights into disk frame format
ipw_disk <- as.disk.frame(data.frame(ipw_denom = ipw_denom_v, ipw_num = ipw_num_v), overwrite = TRUE)
new_path <- tempfile(fileext = ".df")
ipw_disk <- rechunk(ipw_disk, 1, new_path, overwrite = FALSE)
rm(ipw_num_v, ipw_denom_v)
gc()

# Assign weights
med_elig_disk <- med_elig_disk %>%
  mutate(ipw_num = ipw_disk %>% dplyr::select(ipw_num) %>% collect() %>% pull(ipw_num),
         ipw_denom = ipw_disk %>% dplyr::select(ipw_denom) %>% collect() %>% pull(ipw_denom))

rm(list=setdiff(ls(), c("K", "HU_time", "event_name", "med_elig_disk")))
gc()
Sys.sleep(10)

# Update weights to follow censoring due to treatment non-adherence
med_elig_disk <- med_elig_disk %>%
  chunk_group_by(id, trial_num) %>%
  mutate(
    ipw_num2 = ifelse(time == 0, 1,
                      ifelse(treat_b == 1 & time > 0, 1,
                             if_else(treat_b == 0 & time > 0, 1 - ipw_num, NA_real_))),
    ipw_denom2 = ifelse(time == 0, 1,
                        ifelse(treat_b == 1 & time > 0, 1,
                               if_else(treat_b == 0 & time > 0, 1 - ipw_denom, NA_real_)))
  ) %>%
  chunk_ungroup()

# Cumulative product of weights
med_elig_disk <- med_elig_disk %>%
  chunk_group_by(id, trial_num) %>%
  mutate(
    sw_a = cumprod(ipw_num2)/cumprod(ipw_denom2)
  ) %>%
  chunk_ungroup()

# Summary of treatment weights
summary(med_elig_disk %>% dplyr::select(sw_a) %>% collect())
cat("Standard deviation of sw_a:", sd(med_elig_disk %>% dplyr::select(sw_a) %>% collect() %>% pull(sw_a)), "\n")

med_elig_disk %>%
  dplyr::select(sw_a) %>%
  collect() %>%
  { quantile(.$sw_a, probs = c(0.90, 0.95, 0.975, 0.99, 0.999, 0.9999)) }
gc()

# Take product of weight for censoring and treatment for each individual
med_elig_disk <- med_elig_disk %>%
  mutate(sw = sw_a * sw_c)

# Summary of final weights before truncation
summary(med_elig_disk %>% dplyr::select(sw) %>% collect())
med_elig_disk %>%
  dplyr::select(sw) %>%
  collect() %>%
  { quantile(.$sw, probs = c(0.90, 0.95, 0.975, 0.98, 0.99, 0.995, 0.999)) }

### Truncate final stabilized weight at the 98th percentile ###
threshold_99 <- med_elig_disk %>%
  dplyr::select(sw) %>%
  collect() %>%
  { quantile(.$sw, probs = c(0.98)) }

med_elig_disk <- med_elig_disk %>%
  mutate(sw_99 = ifelse(sw > threshold_99, threshold_99, sw))

###  Min, 25th percentile, median, mean, SD, 75th percentile, and max: truncated weights ###
print('Summary of final weights after truncation:')
summary(med_elig_disk %>% dplyr::select(sw_99) %>% collect())
cat("Standard deviation of sw_99:", sd(med_elig_disk %>% dplyr::select(sw_99) %>% collect() %>% pull(sw_99)), "\n")

#################################################################################
#################################################################################

rm(list=setdiff(ls(), c("K", "HU_time", "event_name", "med_elig_disk")))
gc()
Sys.sleep(10)

### Construct a covariate balance plot for the weighted population ###

varlist <- c("Age_cat", "Year_cat", "nordic", "Sex", "educ", "N.endo_cat_b",
             "N.visits_cat_b", "N.drugs_cat_b", "CCIW_cat_b", "ind_cat_b")

treat_0 <- med_elig_disk %>% 
  filter(treat_b == 0) %>% 
  dplyr::select(c(varlist, "sw_99", "Age"))

treat_1 <- med_elig_disk %>% 
  filter(treat_b == 1) %>% 
  dplyr::select(c(varlist, "sw_99", "Age"))

# Define the function to calculate weighted statistics for Age
weighted_stats <- function(data, variable, weight) {
  wtd_mean <- data %>%
    chunk_summarize(wtd_mean = Hmisc::wtd.mean(!!sym(variable), !!sym(weight), na.rm = TRUE)) %>%
    collect() %>%
    pull(wtd_mean)
  
  wtd_variance <- data %>%
    chunk_summarize(wtd_variance = Hmisc::wtd.var(!!sym(variable), !!sym(weight), normwt = TRUE, na.rm = TRUE)) %>%
    collect() %>%
    pull(wtd_variance)
  
  wtd_sd <- sqrt(wtd_variance)
  
  wtd_median <- data %>%
    chunk_summarize(wtd_median = Hmisc::wtd.quantile(!!sym(variable), !!sym(weight), probs = 0.5, normwt = TRUE, na.rm = TRUE)) %>%
    collect() %>%
    pull(wtd_median)
  
  wtd_q25 <- data %>%
    chunk_summarize(wtd_q25 = Hmisc::wtd.quantile(!!sym(variable), !!sym(weight), probs = 0.25, normwt = TRUE, na.rm = TRUE)) %>%
    collect() %>%
    pull(wtd_q25)
  
  wtd_q75 <- data %>%
    chunk_summarize(wtd_q75 = Hmisc::wtd.quantile(!!sym(variable), !!sym(weight), probs = 0.75, normwt = TRUE, na.rm = TRUE)) %>%
    collect() %>%
    pull(wtd_q75)
  
  return(data.frame(
    Mean = wtd_mean,
    SD = wtd_sd,
    Median = wtd_median,
    Q25 = wtd_q25,
    Q75 = wtd_q75
  ))
}

# Calculate weighted statistics for treat_b0
stats_treat_b0 <- weighted_stats(treat_0, "Age", "sw_99")
print("Weighted Statistics for treat_b0:")
print(stats_treat_b0)

# Calculate weighted statistics for treat_b1
stats_treat_b1 <- weighted_stats(treat_1, "Age", "sw_99")
print("Weighted Statistics for treat_b1:")
print(stats_treat_b1)

# Function to compute weighted means and differences
wmean_fctn <- function(df_0, df_1, var) {
  # Calculate the overall sum of weights for df_0 and df_1
  overall_sum_0 <- df_0 %>%
    chunk_summarize(total_sw_99 = sum(sw_99, na.rm = TRUE)) %>%
    collect() %>%
    pull(total_sw_99)
  
  overall_sum_1 <- df_1 %>%
    chunk_summarize(total_sw_99 = sum(sw_99, na.rm = TRUE)) %>%
    collect() %>%
    pull(total_sw_99)
  
  # Calculate the weighted proportions for df_0
  df_0_summary <- df_0 %>%
    chunk_group_by(!!sym(var)) %>%
    chunk_summarize(wprop_0 = sum(sw_99, na.rm = TRUE) / overall_sum_0) %>%
    chunk_ungroup() %>%
    collect()
  
  # Calculate the weighted proportions for df_1
  df_1_summary <- df_1 %>%
    chunk_group_by(!!sym(var)) %>%
    chunk_summarize(wprop_1 = sum(sw_99, na.rm = TRUE) / overall_sum_1) %>%
    chunk_ungroup() %>%
    collect()
  
  # Merge summaries and compute differences
  result <- df_0_summary %>%
    left_join(df_1_summary, by = var, suffix = c("_0", "_1")) %>%
    mutate(md = wprop_1 - wprop_0) %>%
    transmute(wprop_0, wprop_1, var = paste0(var, "_", !!sym(var)), md)
  
  return(result)
}

# Initialize an empty data.table to store results
covplot_w <- data.table()

print('Weighted mean differences:')
# Compute weighted means and differences for all variables in the list
covplot_w <- future_map_dfr(varlist, function(var_name) {
  result <- wmean_fctn(treat_0, treat_1, var_name)
  return(result)
}, .progress = TRUE)

# Convert the result to a data frame and ensure 'md' is numeric
covplot_w <- as.data.frame(covplot_w)
covplot_w$md <- as.numeric(covplot_w$md)

# Function to rename variables for the covariate balance plot
rename_lookup <- c(
  "Age_cat_(64,70]" = "Age_65-70",
  "Age_cat_(70,75]" = "Age_71-75",
  "Age_cat_(75,80]" = "Age_76-80",
  "Age_cat_(80,85]" = "Age_81-85",
  "Age_cat_(85,Inf]" = "Age_86+",
  "Year_cat_(2005,2008]" = "Year_2006-2008",
  "Year_cat_(2008,2011]" = "Year_2009-2011",
  "Year_cat_(2011,2014]" = "Year_2012-2014",
  "Year_cat_(2014,2017]" = "Year_2015-2017",
  "nordic_0" = "nordic_0",
  "nordic_1" = "nordic",
  "Sex_0" = "Sex_0",
  "Sex_1" = "Sex",
  "educ_0" = "educ_0",
  "educ_1" = "educ_1",
  "educ_2" = "educ_2",
  "N.endo_cat_b_(-1,0]" = "N.endo_0",
  "N.endo_cat_b_(0,Inf]" = "N.endo_1+",
  "N.visits_cat_b_(-1,5]" = "N.visits_0-5",
  "N.visits_cat_b_(5,10]" = "N.visits_6-10",
  "N.visits_cat_b_(10,Inf]" = "N.visits_11+",
  "N.drugs_cat_b_(-1,1]" = "N.drugs_0-1",
  "N.drugs_cat_b_(1,4]" = "N.drugs_2-4",
  "N.drugs_cat_b_(4,7]" = "N.drugs_5-7",
  "N.drugs_cat_b_(7,10]" = "N.drugs_8-10",
  "N.drugs_cat_b_(10,13]" = "N.drugs_11-13",
  "N.drugs_cat_b_(13,16]" = "N.drugs_14-16",
  "N.drugs_cat_b_(16,19]" = "N.drugs_17-19",
  "N.drugs_cat_b_(19,Inf]" = "N.drugs_20+",
  "CCIW_cat_b_(-1,2]" = "CCIW_0-2",
  "CCIW_cat_b_(2,4]" = "CCIW_3-4",
  "CCIW_cat_b_(4,Inf]" = "CCIW_5+",
  "ind_cat_b_(-1,2]" = "ind_0-2",
  "ind_cat_b_(2,5]" = "ind_3-5",
  "ind_cat_b_(5,Inf]" = "ind_6+"
)

# Rename variables
covplot_w$var <- rename_lookup[covplot_w$var]

# Remove the smallest level for Sex and nordic variables
covplot_w <- covplot_w[!covplot_w$var %in% c("Sex_0", "nordic_0", "N.endo_0"), ]

# Assign order for selected variables
correct_order <- c(
  "Age_65-70", "Age_71-75", "Age_76-80", "Age_81-85", "Age_86+",
  "Year_2006-2008", "Year_2009-2011", "Year_2012-2014", "Year_2015-2017",
  "nordic", "Sex",
  "educ_0", "educ_1", "educ_2",
  "N.endo_1+",
  "N.visits_0-5", "N.visits_6-10", "N.visits_11+",
  "N.drugs_0-1", "N.drugs_2-4", "N.drugs_5-7", 
  "N.drugs_8-10", "N.drugs_11-13", "N.drugs_14-16", "N.drugs_17-19", "N.drugs_20+",
  "CCIW_0-2", "CCIW_3-4", "CCIW_5+",
  "ind_0-2", "ind_3-5", "ind_6+"
)

# Convert the correct_order to a factor with levels in the specified order
covplot_w$var <- factor(covplot_w$var, levels = correct_order)

# Order the dataframe based on the factor levels
covplot_w_ordered <- covplot_w %>%
  arrange((var))

print(covplot_w_ordered)

# Create plot
covplot_weighted <- ggplot(data = covplot_w_ordered) +
  geom_point(aes(x = md, y = var), color = "steelblue") + scale_x_continuous(limits = c(-0.25, 0.25)) +
  geom_vline(xintercept = 0) +
  labs(y = "Covariates", x = "Mean Difference", title = "Covariate Balance Plot") +
  scale_y_discrete(limits = rev(levels(covplot_w$var)))

# Adjust margin around the image
covplot_weighted <- covplot_weighted + theme(plot.margin = margin(0.9,0.9,0.9,0.9, "cm"))

# Save the plot as .png
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave(paste0("ETT/PPI/Output/covplot_K", K, "_HU", HU_time, "_", event_name, '_GI',".png"), plot = covplot_weighted, width = 8, height = 9)
while (!is.null(dev.list()))  dev.off()

#################################################################################
#################################################################################

### Save a copy of the dataset for bootstraps

selected_columns <- c("id", "trial_num", "id_new", "sw", "sw_99", "treat_b", "treat", 
                      "time", "timesqr", "period", "periodsqr", "eventMC",
                      "Age_cat", "Year_cat", "Sex", "educ", "nordic",
                      "N.endo_cat_b", "N.GI_cat_b", "N.visits_cat_b",
                      "N.drugs_cat_b", "CCIW_cat_b", "ind_cat_b")

if (K %in% c(24, 60)) {
  saveRDS(med_elig_disk %>% dplyr::select(selected_columns) %>% collect(), 
          sprintf("ETT/PPI/Processed_files/PPI_BT_GI_%d_%s.rds", K, event_name), compress = "xz")
}

print("Bootstrap data successfully exported")

### Fit weighted pooled logistic regression with final stabilized weights ###

cat("Total events:", sum(med_elig_disk %>% dplyr::select(eventMC) %>% collect() %>% pull(eventMC)), "\n")
cat("Total events, treatb=0:", sum(med_elig_disk %>% filter(treat_b==0) %>% dplyr::select(eventMC) %>% collect() %>% pull(eventMC)), "\n")
cat("Total events, treatb=1:", sum(med_elig_disk %>% filter(treat_b==1) %>%  dplyr::select(eventMC) %>% collect() %>% pull(eventMC)), "\n")
cat("Total events, treat=0:", sum(med_elig_disk %>% filter(treat==0) %>% dplyr::select(eventMC) %>% collect() %>% pull(eventMC)), "\n")
cat("Total events, treat=1:", sum(med_elig_disk %>% filter(treat==1) %>%  dplyr::select(eventMC) %>% collect() %>% pull(eventMC)), "\n")

# Include product terms between time and treatment
fit.pool1 <- speedglm(formula = eventMC==1 ~ treat_b + time + timesqr + period + periodsqr + Year_cat +
                        I(treat_b*time) +  I(treat_b*timesqr) +
                        as.factor(Age_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) + 
                        as.factor(N.endo_cat_b) + as.factor(N.visits_cat_b) + 
                        as.factor(N.drugs_cat_b) + as.factor(CCIW_cat_b) + as.factor(ind_cat_b),
                      family = binomial(link = 'logit'),
                      data = med_elig_disk,
                      weights = med_elig_disk$sw_99, sparse=FALSE)

# Print results 
summary(fit.pool1)
gc()

### Transform estimates to risks at each time point in each group ###

# Create dataset with all time points for each individual under each treatment level

### No treatment arm

# Filter the rows where time is 0
treat0 <- med_elig_disk %>%
  filter(time == 0)

new_path <- tempfile(fileext = ".df")
treat0 <- rechunk(treat0, 8, new_path, overwrite = FALSE, shardby = "id")

# Number of unqiue individuals in treat0
treat0_N <- nrow(treat0)

# Create a dataframe with K rows, each having time = 0
replication_df <- tibble::tibble(time = rep(0, K))

# Save the replication_df as a disk.frame
replication_df <- as.disk.frame(replication_df)

# Perform a cross join to replicate rows
treat0 <- treat0 %>%
  full_join(replication_df, by = "time", merge_by_chunk_id=FALSE)

# Arrange
treat0 <- treat0 %>% chunk_arrange(id, trial_num, time)

treat0 <- treat0 %>%
  mutate(time = rep(seq(0, K-1), treat0_N),
         timesqr = time^2,
         treat_b = 0)

# Extract predicted values from the pooled logistic regression model for each person-time row
# Predicted values correspond to discrete-time hazards
treat0 <- treat0 %>%
  mutate(p.event0 = predict(fit.pool1, ., type = "response"))

# Arrange data
treat0 <- treat0 %>% chunk_arrange(id, trial_num, time)

# Group by id_new and calculate survival probabilities from discrete-time hazards
treat0 <- treat0 %>%
  chunk_group_by(id, trial_num) %>%
  mutate(surv0 = cumprod(1 - p.event0)) %>%
  chunk_ungroup()

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
treat0 <- treat0 %>%
  mutate(risk0 = 1 - surv0)

# Get the mean in each treatment group at each time point
risk_0 <- treat0 %>%
  chunk_group_by(treat_b, time) %>%
  chunk_summarize(risk0 = mean(risk0)) %>%
  dplyr::select(treat_b, time, risk0) %>%
  chunk_ungroup() %>%
  collect()

delete(treat0)

### Treatment arm

# Filter the rows where time is 0
treat1 <- med_elig_disk %>%
  filter(time == 0)

new_path <- tempfile(fileext = ".df")
treat1 <- rechunk(treat1, 8, new_path, overwrite = FALSE, shardby = "id")

# Number of unqiue individuals in treat1
treat1_N <- nrow(treat1)

# Create a dataframe with K rows, each having time = 0
replication_df <- tibble::tibble(time = rep(0, K))

# Save the replication_df as a disk.frame
replication_df <- as.disk.frame(replication_df)

# Perform a cross join to replicate rows
treat1 <- treat1 %>%
  full_join(replication_df, by = "time", merge_by_chunk_id=FALSE)

## Arrange
treat1 <- treat1 %>% chunk_arrange(id, trial_num, time)

treat1 <- treat1 %>%
  mutate(time = rep(seq(0, K-1), treat1_N),
         timesqr = time^2,
         treat_b = 1)

# Extract predicted values from the pooled logistic regression model for each person-time row
# Predicted values correspond to discrete-time hazards
treat1 <- treat1 %>%
  mutate(p.event1 = predict(fit.pool1, ., type = "response"))

# Group by id_new and calculate survival probabilities from discrete-time hazards
treat1 <- treat1 %>%
  chunk_group_by(id, trial_num) %>%
  mutate(surv1 = cumprod(1 - p.event1)) %>%
  chunk_ungroup()

# Estimate risks from survival probabilities
# Risk = 1 - S(t)
treat1 <- treat1 %>%
  mutate(risk1 = 1 - surv1)

# Get the mean in each treatment group at each time point
risk_1 <- treat1 %>%
  chunk_group_by(treat_b, time) %>%
  chunk_summarize(risk1 = mean(risk1)) %>%
  dplyr::select(treat_b, time, risk1) %>%
  chunk_ungroup() %>%
  collect()

delete(treat1)

# Prepare data
graph.pred <- merge(risk_0, risk_1, by=c("time"))

# Edit data frame to reflect that risks are estimated at the end of each interval
graph.pred$time_0 <- graph.pred$time + 1
zero <- data.frame(cbind(0,0,0,1,0,0))
zero <- setNames(zero,names(graph.pred))
graph <- rbind(zero, graph.pred)

if (K==60) {
  saveRDS(graph, sprintf("ETT/PPI/Processed_files/PPI_graph_GI_%s.rds", event_name))
}

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
                color = "PPI"),
            linewidth = 1.5) + 
  geom_line(aes(y = risk0, # create line for no vaccine group
                color = "Non-PPI"),
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
                     breaks=c('Non-PPI', 'PPI'))  

# Increase the white margin around the image
plot.plr <- plot.plr + theme(plot.margin = margin(0.8,0.8,0.8,0.8, "cm"))

# Plot
plot.plr

# Save the plot as .png
ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
ggsave(paste0("ETT/PPI/Output/CI_K", K, "_HU", HU_time, "_", event_name, '_GI', ".png"), plot = plot.plr, width = 8, height = 8)
while (!is.null(dev.list()))  dev.off()

