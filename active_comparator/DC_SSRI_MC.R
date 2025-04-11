# Load required packages

if (!require("arrow")) install.packages("arrow")
library(arrow)
if (!require("dbplyr")) install.packages("dbplyr")
library(dbplyr)
if (!require("data.table")) install.packages("data.table")
library(data.table)
if (!require("rio")) install.packages("rio")
library(rio)
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)
if (!require("lubridate")) install.packages("lubridate")
library(lubridate)
if (!require("tidyr")) install.packages("tidyr")
library(tidyr)
if (!require("zoo")) install.packages("zoo")
library(zoo)
if (!require("survival")) install.packages("survival")
library(survival)
if (!require("parallel")) install.packages("parallel")
library(parallel)

setwd('~/Medications')
set.seed(332)

event_name <- "eventMC"
event_time <- "timeMC"
event_enddate <- "enddateMC"
event_suffix <- 'MC'

# event_name <- "eventNM"
# event_time <- "timeNM"
# event_enddate <- "enddateNM"
# event_suffix <- 'NM'

############### User-defined functions ######################
#############################################################

# Define a function to merge two data.tables by time
merge_by_time <- function(main_df, ref_df, main_date_col, ref_date_col, main_by_col, ref_by_col, roll_value) {
  # Convert data.tables if not already
  main_df <- as.data.table(main_df)
  ref_df <- as.data.table(ref_df)
  
  # Create join_time column
  main_df[, join_time := get(main_date_col)]
  ref_df[, join_time := get(ref_date_col)]
  
  # Set keys
  setkeyv(main_df, c(main_by_col, "join_time"))
  setkeyv(ref_df, c(ref_by_col, "join_time"))
  
  # Perform update in place
  merged_df <- ref_df[main_df, roll = roll_value]
  
  return(merged_df)
}

# Function to randomly select a sample population from treated_data and untreated_data
# Arguments:
#   - treated_data: The dataframe containing the treated group
#   - untreated_data: The dataframe containing the untreated group
#   - percent_treated: The percentage of people to select from treated_data (default is 50)
#   - ratio: The ratio indicating how many people per 1 person in the treated_data should be selected from the untreated_data (default is 1)
# Returns:
#   - A dataframe containing the sampled population

random_sample <- function(treated_data, untreated_data, percent_treated = 50, ratio = 1) {
  
  # selects 1:1 with same sample size equal to that of untreated population when not enough samples are present
  # in the untreated group
  
  treated_count_before <- nrow(treated_data)
  untreated_count_before <- nrow(untreated_data)
  
  treated_data$treated <- TRUE
  untreated_data$treated <- FALSE
  
  treated_sample_size <- round(percent_treated / 100 * treated_count_before)
 
   # considering case when untreated cohort is smaller than treated
  untreated_sample_size <- min((treated_sample_size * ratio), untreated_count_before)
  
  treated_sample_size <- min(untreated_sample_size, treated_sample_size)
  
  treated_sample <- treated_data %>%
    sample_n(treated_sample_size) # no replacement
  
  untreated_sample <- untreated_data %>%
    sample_n(untreated_sample_size) # no replacement
  
  sampled_df <- rbind(treated_sample, untreated_sample)
  
  treated_count_after <- sum(sampled_df$treated)
  untreated_count_after <- nrow(sampled_df) - treated_count_after
  
  cat("Before selection:\n")
  cat("Treated users:", treated_count_before, "\n")
  cat("Untreated non-users:", untreated_count_before, "\n\n")
  
  cat("After selection:\n")
  cat("Treated users:", treated_count_after, "\n")
  cat("Untreated non-users:", untreated_count_after, "\n\n")
  
  return(sampled_df)
}

# Splits survival data into chunks and applies survSplit to reduce memory usage
process_data_in_chunks <- function(data, chunk_size, cut, end, event, id) {
  num_rows <- nrow(data)
  num_chunks <- ceiling(num_rows / chunk_size)
  
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, num_rows)
    
    chunk_data <- data[start_idx:end_idx, ]
    print(nrow(chunk_data))
    
    chunk_split <- survSplit(chunk_data[chunk_data[event_time] > 0, ], cut = cut, end = end, 
                             start = "Tstart", event = event, id = id)
    
    if (i == 1) {
      result <- chunk_split
    } else {
      result <- rbind(result, chunk_split)
      rm(chunk_data, chunk_split)
      gc()
    }
  }
  
  return(result)
}

# Marks the 12-month period following each treatment as "1", skipping overlapping windows
twelvemonth_lag <- function(input_vector) {
  n <- length(input_vector)
  output_vector <- rep(0, n)  # Initialize output vector with zeros
  
  i <- 1
  while (i<n) {
    if (input_vector[i] == 1) {
      # Set the next 11 rows to 1, excluding the current row
      output_vector[(i + 1):(min(i + 11, n))] <- 1
      # Skip the next 11 rows
      i <- i + 12
    }
    
    else {i=i+1}
  }
  
  return(output_vector)
}

# Applies the 12-month lag function within each individual (Pid)
apply_twelvemonth_lag <- function(data) {
  result_df <- data %>%
    group_by(Pid) %>%
    mutate(med_twelvemonths = twelvemonth_lag(treat)) %>%
    ungroup()
  
  return(result_df)
}

##  Function to create rolling counts for endo, drugs and visits
create_rolling_sequence <- function(data_frame) {
  # Create a new column 'month_year' based on the 'event_date' column
  data_frame$month_year <- format(data_frame$event_date, "%m_%Y")
  
  # Count the occurrences of each unique combination of 'Pid' and 'month_year'
  counts_df <- data_frame %>%
    group_by(Pid, month_year) %>%
    summarise(count = n())
  
  # Create a sequence of all unique 'Pid' and 'event_date' combinations
  sequence_df <- expand.grid(
    Pid = unique(data_frame$Pid),
    event_date = seq(as.Date("2005-01-01"), as.Date("2017-12-01"), by = "month")
  )
  
  # Create a new column 'month_year' based on the 'event_date' column
  sequence_df$month_year <- format(sequence_df$event_date, "%m_%Y")
  
  # Arrange the data frame based on 'Pid' and 'event_date'
  sequence_df <- sequence_df %>% arrange(Pid, event_date)
  
  # Left join 'sequence_df' with 'counts_df' on 'Pid' and 'month_year'
  sequence_df <- sequence_df %>% 
    left_join(counts_df, by=c("Pid", "month_year"))
  
  # Replace NA values in the count column with 0
  sequence_df$count[is.na(sequence_df$count)] <- 0
  
  # Create a new column 'N_rolling_lag' representing the rolling sum of 'count' for the previous 12 rows within the same 'Pid'
  sequence_df$N_rolling_lag <- ave(sequence_df$count, sequence_df$Pid, 
                                   FUN = function(x) zoo::rollapply(x, width = 12, align = "right", fill = 0, sum, na.rm = TRUE))
  
  # Group 'sequence_df' by 'Pid' and create a lagged version of 'N_rolling' to represent counts starting from the previous month
  sequence_df <- sequence_df %>%
    group_by(Pid) %>%
    mutate(N_rolling = lag(N_rolling_lag, n = 1, default = NA)) %>%
    dplyr::select(Pid, month_year, N_rolling)
  
  rm(counts_df)
  
  return(sequence_df)
}

##  Function to create rolling counts for drugs
calculate_unique_atc_counts <- function(data_frame) {
  # Load required libraries
  library(dplyr)
  library(zoo)
  
  # Convert 'event_date' to Date type if it's not already
  if (!inherits(data_frame$event_date, "Date")) {
    data_frame$event_date <- as.Date(data_frame$event_date)
  }
  
  # Create a sequence of months from January 2005 to January 2018
  months_sequence <- seq(as.Date("2005-01-01"), as.Date("2018-01-01"), by = "months")
  
  # Initialize an empty list to store results
  result_list <- list()
  
  # Iterate over each month
  for (month in months_sequence) {
    # Calculate the 12-month window
    start_window <- as.Date(month) - months(12)
    end_window <- as.Date(month) - days(1)
    
    # Subset the dataframe to include only rows within the 12-month window
    subset_df <- data_frame %>%
      filter(event_date >= start_window & event_date <= end_window)
    
    # Count unique ATC codes for each Pid
    counts <- subset_df %>%
      group_by(Pid) %>%
      summarize(unique_atc_count = n_distinct(Atc))
    
    # Create a new column for event_month_year
    counts$event_month_year <- format(as.Date(month), "%m_%Y")
    
    # Append the counts to the result list
    result_list[[length(result_list) + 1]] <- counts
  }
  
  # Combine the results into a single dataframe
  final_df <- do.call(rbind, result_list)
  
  # Ensure all combinations of Pid and event_month_year exist
  final_df <- merge(expand.grid(unique(final_df$Pid), unique(final_df$event_month_year)), final_df, by.x = c("Var1", "Var2"), by.y = c("Pid", "event_month_year"), all.x = TRUE)
  
  # Replace missing counts with 0
  final_df[is.na(final_df)] <- 0
  
  # Rename columns
  colnames(final_df) <- c("Pid", "month_year", "N.drugs")
  
  # Order the dataframe by Pid and event_month_year
  final_df <- final_df[order(final_df$Pid, as.Date(final_df$month_year, format = "%m_%Y")), ]
  
  # Return the final dataframe
  return(final_df)
}

# Extracts and optionally saves/loads first event dates based on ICD code
extract_event_dates <- function(icd_code, patients_data, event_name, path = NULL, read_data = FALSE) {
  if (read_data) {
    # Check if the file exists and read if specified
    if (!is.null(path) && file.exists(path)) {
      event_data <- read_feather(path)
      return(event_data)
    } else {
      stop("Specified path does not exist or data is not available for reading.")
    }
  } else {
    # Extract event dates based on the specified event code
    event_subset <- patients_data[grep(icd_code, patients_data$DIA), ]
    
    # Create data frame for the event dates
    event_data <- bind_rows(event_subset) %>%
      group_by(pid) %>%
      filter(row_number(date_admission) == 1) %>%
      ungroup %>%
      rename(!!event_name := date_admission) %>%
      dplyr::select(pid, !!event_name)
    
    # Write event data frame to feather file
    if (!is.null(path)) {
      write_feather(event_data, path)
    }
    
    return(event_data)
  }
}

# Merges CCI or event data by rolling join on Pid and date, and creates an indicator variable
merge_cci <- function(input_df, date_col, result_df, current_date, new_col_name) {
  input_df <- input_df %>% rename(Pid = pid) %>% distinct(across(all_of(c('Pid', date_col))))
  input_df <- as.data.table(input_df)
  input_df[, join_time := input_df[[date_col]]]
  result_df[, join_time := current_date]
  setkey(input_df, Pid, join_time)
  setkey(result_df, Pid, join_time)
  result_df <- input_df[result_df, roll = Inf]
  result_df[, (new_col_name) := as.integer(!is.na(result_df[[date_col]]) & !is.na(result_df$current_date))]
  
  return(result_df)
}

#############################################################
#############################################################

conn <- DBI::dbConnect(
  odbc::odbc(), 
  .connection_string = "Driver=Driver_Name;server=Server_Name;database=DB_Name;schema=Schema_Name;trusted_connection=yes;TrustServerCertificate=yes")

OUTCOME <- read_feather("Data/OUTCOME.ft")

# Convert dates in OUTCOME df to datetime format
OUTCOME <- OUTCOME %>%
  mutate(
    First_IBD_DATE = as.Date(First_IBD_DATE, format = "%Y%m%d"),
    Second_IBD_DATE = as.Date(Second_IBD_DATE, format = "%Y%m%d"),
    MC_D_DATE = as.Date(MC_D_DATE, format = "%Y%m%d"),
    MC_DATE = as.Date(MC_D_DATE, format = "%Y-%m-%d"),
    NM_DATE = as.Date(NM_D_DATE, format = "%Y%m%d")
  )

outcome <- dplyr::select(OUTCOME, PID, MC_D_DATE, MC_DATE, MC_SNOMED, MC_TOPO, MC_Type, NM_D_DATE, NM_DATE, NM_Type) %>% 
  collect %>% filter(!is.na(MC_D_DATE) | !is.na(NM_D_DATE))

# Patient demographic data
COHORT <- read_feather("Data/COHORT.ft")

enddates <- dplyr::select(COHORT,pid,birth_year, birth_month, x_death_date, x_last_imm_date, x_last_emi_date) %>% collect
enddates <- mutate(enddates, migr.date=if_else(as.Date(x_last_emi_date) > as.Date(x_last_imm_date) | 
                                                 (!is.na(x_last_emi_date) & is.na(x_last_imm_date)),x_last_emi_date,NA_character_,missing=NA_character_))

enddates <- inner_join(enddates, outcome,by=c("pid"="PID"))

enddates <- enddates %>%
  mutate(NMdate = as.numeric(gsub("-", "", NM_DATE)),
         MCdate = as.numeric(gsub("-", "", MC_DATE)),
         migrationdate = as.numeric(gsub("-", "", migr.date)),
         deathdate = as.numeric(gsub("-", "", x_death_date)),
         migr.date=as.Date(migr.date),
         x_death_date = as.Date(x_death_date))

enddates <- dplyr::select(enddates, pid, birth_year, MCdate, NMdate,migrationdate, deathdate, MC_DATE, NM_DATE, migr.date, x_death_date)

### Medications file

### SSRI cohort

# Construct the SQL query to extract all SSRI records from the drugs table
query <- "WITH RankedData AS (
              SELECT dm.Pid, dm.Atc, dm.Age, dm.Sex, dm.startdate, dm.enddate, dm.Antal,
                     ROW_NUMBER() OVER (PARTITION BY dm.Pid, dm.startdate ORDER BY dm.enddate DESC) AS RowNum
              FROM SCHEMA_X.DATASET_Y dm
              WHERE dm.Age >= 64 
                AND dm.startdate >= '2005-01-01' 
                AND dm.startdate < '2018-01-01'
                AND dm.startdate <= dm.enddate
                AND dm.Atc LIKE 'N06AB%'
          )
          SELECT Pid, Atc, Age, Sex, startdate, enddate, Antal
          FROM RankedData
          WHERE RowNum = 1;"

# Execute the SQL query
SSRI_cohort <- DBI::dbGetQuery(conn, query)

print(dim(SSRI_cohort))

SSRI_cohort <- SSRI_cohort %>%
  mutate(date_admission=startdate, SSRI.begin=startdate,SSRI.end=enddate, year = year(date_admission)) %>%
  dplyr::select(Pid, Atc, startdate, enddate, Age, Sex, year, date_admission, SSRI.begin, SSRI.end) %>%
  arrange(Pid, SSRI.begin, SSRI.end)

# Keep all medication records for treatment assignment
SSRI_all <- SSRI_cohort  %>%
  rename(SSRI_startdate = startdate, SSRI_endate = enddate) %>% 
  dplyr::select(Pid, SSRI_startdate, SSRI_endate)

### MIR cohort

# Construct the SQL query to extract all MIR records from the drugs table
query <- "WITH RankedData AS (
              SELECT dm.Pid, dm.Atc, dm.Age, dm.Sex, dm.startdate, dm.enddate, dm.Antal,
                     ROW_NUMBER() OVER (PARTITION BY dm.Pid, dm.startdate ORDER BY dm.enddate DESC) AS RowNum
              FROM SCHEMA_X.DATASET_Y dm
              WHERE dm.Age >= 64 
                AND dm.startdate >= '2005-01-01' 
                AND dm.startdate < '2018-01-01'
                AND dm.startdate <= dm.enddate
                AND dm.Atc LIKE 'N06AX11%'
          )
          SELECT Pid, Atc, Age, Sex, startdate, enddate, Antal
          FROM RankedData
          WHERE RowNum = 1;"

# Execute the SQL query
MIR_cohort <- DBI::dbGetQuery(conn, query)

print(dim(MIR_cohort))

MIR_cohort <- MIR_cohort %>%
  mutate(date_admission=startdate, MIR.begin=startdate,MIR.end=enddate, year = year(date_admission)) %>%
  dplyr::select(Pid, Atc, startdate, enddate, Age, Sex, year, date_admission, MIR.begin, MIR.end) %>%
  arrange(Pid, MIR.begin, MIR.end)

# Keep all medication records for treatment assignment
MIR_all <- MIR_cohort  %>%
  rename(MIR_startdate = startdate, MIR_endate = enddate) %>% 
  dplyr::select(Pid, MIR_startdate, MIR_endate)


###### SSRI cohort - exclusions ###### 

# Add first ever MIR use date to SSRI cohort

MIR_cov <- MIR_all %>%
  # Select relevant columns
  select(Pid, MIR_startdate, MIR_endate) %>%
  # Arrange by Pid, and dates
  arrange(Pid, MIR_startdate, MIR_endate) %>%
  # Remove duplicates based on Pid
  distinct(Pid, .keep_all = TRUE) %>%
  rename(comp_startdate = MIR_startdate, comp_enddate = MIR_endate) %>%
  ungroup()

# Merging MIR_cov with SSRI_cohort using Pid with a left join
SSRI_cohort <- left_join(SSRI_cohort, MIR_cov, by = "Pid")

rm(MIR_cov)
gc()

# Preprocess SSRI cohort to have one row per patient

unique_participants <- SSRI_cohort %>%
  summarise(unique_count = n_distinct(Pid))
unique_participants <- unique_participants$unique_count
cat('Number of participants with SSRI use:', unique_participants, '\n')

SSRI_cohort <- filter(SSRI_cohort, year > 2005 & Age > 64)
SSRI_cohort <- subset(SSRI_cohort,select = -year)
exclusion_count <-  unique_participants - length(unique(SSRI_cohort$Pid))
cat('Exclusion, SSRI usage only prior to 2006 or under the age of 65, N:', exclusion_count, '\n')

SSRI_cohort <- SSRI_cohort %>%
  arrange(Pid, startdate) %>%  # Arrange by Pid and startdate
  distinct(Pid, .keep_all = TRUE)  # Drop duplicates, keeping the first occurrence

## Prior MIR
## Prior MC
## Prior IBD
## Prior SSRI
## Prior MAO-I

# Count the rows to exclude where MIR_startdate is on or before date_admission
count_exclusions <- SSRI_cohort %>%
  filter(comp_startdate <= date_admission) %>%
  nrow()

# Print the count of exclusions
cat('Exclusion, prior MIR use, N:', count_exclusions, '\n')

# Exclude rows where MIR_startdate is on or before date_admission
SSRI_cohort <- SSRI_cohort %>%
  filter((comp_startdate > date_admission) | is.na(comp_startdate))

## Create a lag variable for admission date for merging
SSRI_cohort$date_admission_lag <- SSRI_cohort$date_admission - 1

## SSRI usage 6 months prior to baseline

SSRI_cov <- dplyr::select(SSRI_all, Pid, SSRI_startdate, SSRI_endate)

SSRI_cohort <- merge_by_time(main_df = SSRI_cohort, ref_df = SSRI_cov, main_date_col = "date_admission_lag", ref_date_col = "SSRI_startdate", 
                            main_by_col = "Pid", ref_by_col = "Pid", roll_value = Inf)

# Count the number of patients with SSRI before 6 months of index date 
# (checking for SSRI alone is enough since these patients doesn't have a history of MIR usage)
missing_count <- sum(difftime(SSRI_cohort$date_admission, SSRI_cohort$SSRI_endate, units = "days")<6*30.5, na.rm=T)

# Print the number of patients with SSRI before 6 months of index date
cat('Exclusion SSRI within 6 months of index date, N:', missing_count, '\n')

# Exclude patients with SSRI before 6 months of index date (SSRI on the same day as index date are not excluded)
SSRI_cohort <- SSRI_cohort[!(difftime(SSRI_cohort$date_admission, SSRI_cohort$SSRI_endate, units = "days")<6*30.5) | is.na(SSRI_endate)]
SSRI_cohort <- SSRI_cohort[, c('SSRI_startdate', 'SSRI_endate', 'join_time') := NULL]
SSRI_cohort <- as_tibble(SSRI_cohort)

## MAO usage 6 months prior to baseline

# Get the unique Pids from SSRI_cohort data frame
SSRI_pids <- unique(SSRI_cohort$Pid)

# Drop the temporary table if it exists
tryCatch({
  DBI::dbExecute(conn, "BEGIN TRY DROP TABLE #TempPids END TRY BEGIN CATCH END CATCH")
}, error = function(e) {})

# Create a temporary table in the database to store the Pids
DBI::dbExecute(conn, "CREATE TABLE #TempPids (Pid INT)")
DBI::dbWriteTable(conn, "#TempPids", data.frame(Pid = SSRI_pids), overwrite = TRUE)

# Construct the SQL query to join the temporary table with DRUG_MAIN
query <- "SELECT dm.Pid, dm.Atc, dm.Age, dm.Sex, dm.startdate, dm.enddate, dm.Antal
          FROM SCHEMA_X.DATASET_Y dm
          INNER JOIN #TempPids tl ON dm.Pid = tl.Pid
          WHERE dm.Age >= 64 
            AND dm.startdate >= '2005-01-01' AND dm.startdate < '2018-01-01'
            AND dm.startdate <= dm.enddate
            AND dm.Atc LIKE 'N04BD%'"

# Execute the SQL query
mao_all <- DBI::dbGetQuery(conn, query)

# Drop the temporary SQL table
DBI::dbExecute(conn, "DROP TABLE #TempPids")

rm(SSRI_pids)
gc()

# Prepare MAO data
mao_all <- mao_all %>%
  rename(mao_startdate=startdate, mao_enddate=enddate) %>% 
  dplyr::select(Pid, mao_startdate, mao_enddate)  %>% 
  mutate(mao_enddate = as.Date(mao_enddate)) %>%
  distinct(Pid, mao_startdate, .keep_all = TRUE)

mao_cov <- dplyr::select(mao_all, Pid, mao_startdate, mao_enddate)

# Merge MAO data to SSRI cohort
SSRI_cohort <- merge_by_time(main_df = SSRI_cohort, ref_df = mao_cov, main_date_col = "date_admission_lag", ref_date_col = "mao_startdate", 
                             main_by_col = "Pid", ref_by_col = "Pid", roll_value = Inf)

# Count the number of patients with MAO before 6 months of index date (this takes into conderation both start and end dates of MAO)
missing_count <- sum((difftime(SSRI_cohort$date_admission, SSRI_cohort$mao_enddate, units = "days")<6*30.5)  , na.rm=T)

# Print the number of patients with MAO before 6 months of index date
cat('Exclusion, MAO within 6 months of index date, N:', missing_count, '\n')

rm(mao_all, mao_cov)

# Exclude patients with MAO before 6 months of index date (MAO on the same day as index date are not excluded)
SSRI_cohort <- SSRI_cohort[!(difftime(SSRI_cohort$date_admission, SSRI_cohort$mao_enddate, units = "days")<6*30.5) | is.na(mao_enddate)]
SSRI_cohort <- SSRI_cohort[, c('mao_startdate', 'mao_enddate', 'join_time') := NULL]
SSRI_cohort <- as_tibble(SSRI_cohort)

# Add IBD & MC
SSRI_cohort <- left_join(SSRI_cohort,dplyr::select(OUTCOME,PID,First_IBD_DATE, Second_IBD_DATE, 
                                                 IBD_Type,MC_D_DATE, MC_Type, NM_D_DATE),by=c("Pid"="PID"),copy=T)

# Exclude IBD prior to baseline
exclusion_count <- sum(!(as.numeric(SSRI_cohort$date_admission) - as.numeric(SSRI_cohort$First_IBD_DATE) < 0 | is.na(SSRI_cohort$First_IBD_DATE)))
cat('Exclusion, prior IBD, N:', exclusion_count, '\n')
SSRI_cohort <- filter(SSRI_cohort,as.numeric(date_admission) - as.numeric(First_IBD_DATE) < 0 | is.na(First_IBD_DATE))

# Exclude MC prior to baseline
exclusion_count <- sum(!(as.numeric(SSRI_cohort$date_admission) - as.numeric(SSRI_cohort$MC_D_DATE) < 0 | is.na(SSRI_cohort$MC_D_DATE)))
cat('Exclusion, prior MC, N:', exclusion_count, '\n')
SSRI_cohort <- filter(SSRI_cohort,as.numeric(date_admission) - as.numeric(MC_D_DATE) < 0 | is.na(MC_D_DATE))

# Add Death and Migration dates
SSRI_cohort <- left_join(SSRI_cohort, enddates, by=c("Pid"="pid"))

# Exclude Migration prior to baseline (if any)
exclusion_count <- sum(!(as.numeric(SSRI_cohort$date_admission) - as.numeric(SSRI_cohort$migr.date) < 0 | is.na(SSRI_cohort$migr.date)))
cat('Exclusion, migration on or prior to date of admission, N:', exclusion_count, '\n')
SSRI_cohort <- filter(SSRI_cohort,as.numeric(date_admission) - as.numeric(migr.date) < 0 | is.na(migr.date))

# Exclude Death on index date
exclusion_count <- sum(!(as.numeric(SSRI_cohort$date_admission) - as.numeric(SSRI_cohort$x_death_date) < 0 | is.na(SSRI_cohort$x_death_date)))
cat('Exclusion, death on or prior to date of admission, N:', exclusion_count, '\n')
SSRI_cohort <- filter(SSRI_cohort,as.numeric(date_admission) - as.numeric(x_death_date) < 0 | is.na(x_death_date))

###### MIR cohort - exclusions ###### 

# Add first ever SSRI use date to MIR cohort

SSRI_cov <- SSRI_all %>%
  # Select relevant columns
  select(Pid, SSRI_startdate, SSRI_endate) %>%
  # Arrange by Pid, and dates
  arrange(Pid, SSRI_startdate, SSRI_endate) %>%
  # Remove duplicates based on Pid
  distinct(Pid, .keep_all = TRUE) %>%
  rename(comp_startdate = SSRI_startdate, comp_enddate = SSRI_endate) %>%
  ungroup()

# Merging SSRI dates to MIR
MIR_cohort <- left_join(MIR_cohort, SSRI_cov, by = "Pid")

rm(SSRI_cov)

# Preprocess MIR cohort to have one row per patient
unique_participants <- MIR_cohort %>%
  summarise(unique_count = n_distinct(Pid))
unique_participants <- unique_participants$unique_count
cat('Number of participants with MIR use:', unique_participants, '\n')

MIR_cohort <- filter(MIR_cohort, year > 2005 & Age > 64)
MIR_cohort <- subset(MIR_cohort,select = -year)
exclusion_count <-  unique_participants - length(unique(MIR_cohort$Pid))
cat('Exclusion, MIR usage only prior to 2006 or under the age of 65, N:', exclusion_count, '\n')

MIR_cohort <- MIR_cohort %>%
  arrange(Pid, startdate) %>%  # Arrange by Pid and startdate
  distinct(Pid, .keep_all = TRUE)  # Drop duplicates, keeping the first occurrence

## Prior SSRI
## Prior MC
## Prior IBD
## Prior MIR
## Prior MAO

# Count the rows to exclude where MIR_startdate is on or before date_admission
count_exclusions <- MIR_cohort %>%
  filter(comp_startdate <= date_admission) %>%
  nrow()

# Print the count of exclusions
cat('Exclusion, prior SSRI use, N:', count_exclusions, '\n')

# Exclude rows where SSRI_startdate is on or before date_admission
MIR_cohort <- MIR_cohort %>%
  filter((comp_startdate > date_admission) | is.na(comp_startdate))

## Create a lag variable for admission date for merging
MIR_cohort$date_admission_lag <- MIR_cohort$date_admission - 1

## MIR usage 6 months prior to baseline

MIR_cov <- dplyr::select(MIR_all, Pid, MIR_startdate, MIR_endate)

MIR_cohort <- merge_by_time(main_df = MIR_cohort, ref_df = MIR_cov, main_date_col = "date_admission_lag", ref_date_col = "MIR_startdate", 
                            main_by_col = "Pid", ref_by_col = "Pid", roll_value = Inf)

# Count the number of patients with MIR before 6 months of index date
missing_count <- sum(difftime(MIR_cohort$date_admission, MIR_cohort$MIR_endate, units = "days")<6*30.5, na.rm=T)

# Print the number of patients with MIR before 6 months of index date
cat('Exclusion MIR within 6 months of index date, N:', missing_count, '\n')

# Exclude patients with MIR before 6 months of index date (MIR on the same day as index date are not excluded)
MIR_cohort <- MIR_cohort[!(difftime(MIR_cohort$date_admission, MIR_cohort$MIR_endate, units = "days")<6*30.5) | is.na(MIR_endate)]
MIR_cohort <- MIR_cohort[, c('MIR_startdate', 'MIR_endate', 'join_time') := NULL]
MIR_cohort <- as_tibble(MIR_cohort)


## MAO

# Get the unique Pids from MIR_cohort data frame
MIR_pids <- unique(MIR_cohort$Pid)

# Drop the temporary table if it exists
tryCatch({
  DBI::dbExecute(conn, "BEGIN TRY DROP TABLE #TempPids END TRY BEGIN CATCH END CATCH")
}, error = function(e) {})

# Create a temporary table in the database to store the Pids
DBI::dbExecute(conn, "CREATE TABLE #TempPids (Pid INT)")
DBI::dbWriteTable(conn, "#TempPids", data.frame(Pid = MIR_pids), overwrite = TRUE)

# Construct the SQL query to join the temporary table with DRUG_MAIN
query <- "SELECT dm.Pid, dm.Atc, dm.Age, dm.Sex, dm.startdate, dm.enddate, dm.Antal
          FROM SCHEMA_X.DATASET_Y dm
          INNER JOIN #TempPids tl ON dm.Pid = tl.Pid
          WHERE dm.Age >= 64 
            AND dm.startdate >= '2005-01-01' AND dm.startdate < '2018-01-01'
            AND dm.startdate <= dm.enddate
            AND dm.Atc LIKE 'N04BD%'"

# Execute the SQL query
mao_all <- DBI::dbGetQuery(conn, query)

# Drop the temporary SQL table
DBI::dbExecute(conn, "DROP TABLE #TempPids")

rm(MIR_pids)
gc()

# Prepare MAO data
mao_all <- mao_all %>%
  rename(mao_startdate=startdate, mao_enddate=enddate) %>% 
  dplyr::select(Pid, mao_startdate, mao_enddate)  %>% 
  mutate(mao_enddate = as.Date(mao_enddate)) %>%
  distinct(Pid, mao_startdate, .keep_all = TRUE)

mao_cov <- dplyr::select(mao_all, Pid, mao_startdate, mao_enddate)

# Merge MAO data with MIR cohort
MIR_cohort <- merge_by_time(main_df = MIR_cohort, ref_df = mao_cov, main_date_col = "date_admission_lag", ref_date_col = "mao_startdate", 
                             main_by_col = "Pid", ref_by_col = "Pid", roll_value = Inf)

# Count the number of patients with MAO before 6 months of index date 
missing_count <- sum(difftime(MIR_cohort$date_admission, MIR_cohort$mao_enddate, units = "days")<6*30.5, na.rm=T)

# Print the number of patients with MAO before 6 months of index date
cat('Exclusion MAO within 6 months of index date, N:', missing_count, '\n')

# Exclude patients with MAO before 6 months of index date (MAO on the same day as index date are not excluded)
MIR_cohort <- MIR_cohort[!(difftime(MIR_cohort$date_admission, MIR_cohort$mao_enddate, units = "days")<6*30.5) | is.na(mao_enddate)]
MIR_cohort <- MIR_cohort[, c('mao_startdate', 'mao_enddate', 'join_time') := NULL]
MIR_cohort <- as_tibble(MIR_cohort)

rm(mao_all, mao_cov)

# Add IBD & MC
MIR_cohort <- left_join(MIR_cohort,dplyr::select(OUTCOME,PID,First_IBD_DATE, Second_IBD_DATE, 
                                                 IBD_Type,MC_D_DATE, MC_Type, NM_D_DATE),by=c("Pid"="PID"),copy=T)

# Exclude IBD prior to baseline
exclusion_count <- sum(!(as.numeric(MIR_cohort$date_admission) - as.numeric(MIR_cohort$First_IBD_DATE) < 0 | is.na(MIR_cohort$First_IBD_DATE)))
cat('Exclusion, prior IBD, N:', exclusion_count, '\n')
MIR_cohort <- filter(MIR_cohort,as.numeric(date_admission) - as.numeric(First_IBD_DATE) < 0 | is.na(First_IBD_DATE))

# Exclude MC prior to baseline
exclusion_count <- sum(!(as.numeric(MIR_cohort$date_admission) - as.numeric(MIR_cohort$MC_D_DATE) < 0 | is.na(MIR_cohort$MC_D_DATE)))
cat('Exclusion, prior MC, N:', exclusion_count, '\n')
MIR_cohort <- filter(MIR_cohort,as.numeric(date_admission) - as.numeric(MC_D_DATE) < 0 | is.na(MC_D_DATE))

# Add Death and Migration dates
MIR_cohort <- left_join(MIR_cohort, enddates, by=c("Pid"="pid"))

# Exclude Migration prior to baseline (if any)
exclusion_count <- sum(!(as.numeric(MIR_cohort$date_admission) - as.numeric(MIR_cohort$migr.date) < 0 | is.na(MIR_cohort$migr.date)))
cat('Exclusion, migration on or prior to date of admission, N:', exclusion_count, '\n')
MIR_cohort <- filter(MIR_cohort,as.numeric(date_admission) - as.numeric(migr.date) < 0 | is.na(migr.date))

# Exclude Death on index date
exclusion_count <- sum(!(as.numeric(MIR_cohort$date_admission) - as.numeric(MIR_cohort$x_death_date) < 0 | is.na(MIR_cohort$x_death_date)))
cat('Exclusion, death on or prior to date of admission, N:', exclusion_count, '\n')
MIR_cohort <- filter(MIR_cohort,as.numeric(date_admission) - as.numeric(x_death_date) < 0 | is.na(x_death_date))

####### Combine cohorts ####### 

cat('Remaining individuals in the user group, N:', nrow(SSRI_cohort), '\n')
cat('Remaining individuals in the non-user group, N:', nrow(MIR_cohort), '\n')

# Rename columns SSRI.begin and SSRI.end to med.begin and med.end
# Add a new column 'med_type' and fill it with 'SSRI' in one step
SSRI_cohort <- SSRI_cohort %>%
  mutate(med.begin = SSRI.begin,
         med.end = SSRI.end,
         med_type = 'SSRI') %>%
  select(-SSRI.begin, -SSRI.end)

# Rename columns MIR.begin and MIR.end to med.begin and med.end
# Add a new column 'med_type' and fill it with 'MIR' in one step
MIR_cohort <- MIR_cohort %>%
  mutate(med.begin = MIR.begin,
         med.end = MIR.end,
         med_type = 'MIR') %>%
  select(-MIR.begin, -MIR.end)

SSRI_combined <- rbind(SSRI_cohort, MIR_cohort)
cat('Number of individuals in the combined cohort, N:', nrow(SSRI_combined), '\n')

rm(outcome, enddates)
rm(SSRI_cohort, MIR_cohort)
gc()
Sys.sleep(60)

####### Add covariates ####### 

## Patients for cohort
patients <- read_feather("Data/patients.ft")
patients <- patients %>%
  filter(date_admission < '2018-01-01')
patients_SSRI_MIR <- inner_join(patients,dplyr::select(SSRI_combined,Pid),by=c("pid"="Pid"))
rm(patients)
gc()


## Colonoscopy
colonoscopy_cov <- patients_SSRI_MIR %>%
  filter(grepl("\\<9011|\\<9023|\\<4688|\\<4689|\\<4674|\\<4684|\\<UJF32|\\<UJF35", patients_SSRI_MIR$op)) %>%
  rename(Pid=pid,colonoscopy_date = date_admission) %>% 
  dplyr::select(Pid, colonoscopy_date)


### Education & country
educountry <- inner_join(COHORT,dplyr::select(SSRI_combined,Pid),by=c("pid"="Pid"),copy=TRUE)
educountry <- collect(dplyr::select(educountry,pid,x_highest_educ, country_group4))

educountry <- mutate(educountry, educ=ifelse(x_highest_educ=="Data not available",0,
                                             ifelse(x_highest_educ=="Pre-secondary education shorter than 9 years"|x_highest_educ=="Pre-secondary education 9 years",1,
                                                    ifelse(x_highest_educ=="Postgraduate education"|x_highest_educ=="Post-secondary education shorter than 3 years",2,3))))
educountry <- mutate(educountry, nordic=ifelse(country_group4=="Norden utom Sverige och Finland"|country_group4=="Finland"|country_group4=="Sverige",1,0))
educountry <- dplyr::select(educountry,pid,educ,nordic)        

SSRI_combined <- left_join(SSRI_combined, educountry, by=c("Pid"="pid"))

cat("Number of missing values in nordic (assign 1):", sum(is.na(SSRI_combined$nordic)), "\n")

# Convert NAs to 1 in the 'nordic' column (if any)
SSRI_combined$nordic[is.na(SSRI_combined$nordic)] <- 1

rm(educountry, COHORT)
gc()

### Study end dates
SSRI_combined <- SSRI_combined %>%
  mutate(
    comp_startdate_c = format(as.Date(comp_startdate, format = "%Y-%m-%d"), "%Y%m%d"),
    comp_enddate_c = format(as.Date(comp_enddate, format = "%Y-%m-%d"), "%Y%m%d")
  )

SSRI_combined$enddateMC <- pmin(SSRI_combined$MCdate, SSRI_combined$comp_startdate_c, SSRI_combined$migrationdate, SSRI_combined$deathdate, 20171231,na.rm=TRUE) #
SSRI_combined$eventMC   <- if_else(SSRI_combined$enddateMC==SSRI_combined$MCdate,1,0,missing=0)
SSRI_combined$enddateMC <- as.Date(paste(substr(SSRI_combined$enddateMC,1,4),'-',substr(SSRI_combined$enddateMC,5,6),'-',substr(SSRI_combined$enddateMC,7,8),sep=''))
SSRI_combined$timeMC    <- as.numeric(SSRI_combined$enddateMC - as.Date(SSRI_combined$date_admission))
SSRI_combined$impDMC    <- if_else(SSRI_combined$timeMC<=0,1,0)
SSRI_combined$timeMC    <- if_else(SSRI_combined$timeMC==0,1,SSRI_combined$timeMC)
SSRI_combined$time_yrMC <- SSRI_combined$timeMC/365.25
SSRI_combined <- mutate(SSRI_combined, ExcludeMC=ifelse(impDMC==1, 1, 0))

# Set NM date to NA if NM occurs after MC
SSRI_combined$NMdate[SSRI_combined$NM_DATE > SSRI_combined$MCdate] <- NA
# Set NM date to NA if NM occurs before baseline
SSRI_combined$NMdate[SSRI_combined$NM_DATE < SSRI_combined$date_admission] <- NA

# Set NM date to NA if NM occurs after MC
SSRI_combined$NMdate[SSRI_combined$NMdate > SSRI_combined$MCdate] <- NA
SSRI_combined$enddateNM <- pmin(SSRI_combined$NMdate, SSRI_combined$comp_startdate_c, SSRI_combined$migrationdate, SSRI_combined$deathdate, 20171231,na.rm=TRUE) #
SSRI_combined$eventNM   <- if_else(SSRI_combined$enddateNM==SSRI_combined$NMdate,1,0,missing=0)
SSRI_combined$enddateNM <- as.Date(paste(substr(SSRI_combined$enddateNM,1,4),'-',substr(SSRI_combined$enddateNM,5,6),'-',substr(SSRI_combined$enddateNM,7,8),sep=''))
SSRI_combined$timeNM    <- as.numeric(SSRI_combined$enddateNM - as.Date(SSRI_combined$date_admission))
SSRI_combined$impDNM    <- if_else(SSRI_combined$timeNM<=0,1,0)
SSRI_combined$timeNM    <- if_else(SSRI_combined$timeNM==0,1,SSRI_combined$timeNM)
SSRI_combined$time_yrNM <- SSRI_combined$timeNM/365.25
SSRI_combined <- mutate(SSRI_combined, ExcludeNM=ifelse(impDNM==1, 1, 0))

SSRI_data <- dplyr::select(SSRI_combined,Pid, date_admission, med.begin, med.end, med_type,
                          comp_startdate, comp_enddate,
                          Age, Sex, educ, nordic,
                          MCdate, First_IBD_DATE, Second_IBD_DATE, eventMC,timeMC, enddateMC, MC_Type,
                          NMdate, enddateNM, eventNM, timeNM, migrationdate, deathdate)

rm(SSRI_combined)
gc()
Sys.sleep(60)

### Convert data to long format
t_events<-seq(0,max(SSRI_data$timeMC),30.5)
times<-data.frame("tevent"=t_events,"ID_t"=seq(1:length(t_events)))

### MC & NM Outcome ### 

cat('Excluding past outcome events, N:', sum(SSRI_data[event_time]<=0), '\n')
SSRI_data <- SSRI_data[SSRI_data[event_time] > 0, ]

chunk_size <- 50000  # Adjust the chunk size as needed
data_long <- process_data_in_chunks(SSRI_data, chunk_size, t_events, event_time, event_name, "ID")

data_long <- data_long %>%
  group_by(Pid) %>%
  arrange(Pid, date_admission) %>%
  mutate(current_date = date_admission + (row_number() - 1) * 30.5,
         cal_time=(row_number()-1), cal_timesqr=(row_number()-1)**2,
         tstart=(row_number()-1), tstop=row_number())


## MAO-I

# Get the unique Pids from data_long data frame
data_long_pids <- unique(data_long$Pid)

# Drop the temporary table if it exists
tryCatch({
  DBI::dbExecute(conn, "BEGIN TRY DROP TABLE #TempPids END TRY BEGIN CATCH END CATCH")
}, error = function(e) {})

# Create a temporary table in the database to store the Pids
DBI::dbExecute(conn, "CREATE TABLE #TempPids (Pid INT)")
DBI::dbWriteTable(conn, "#TempPids", data.frame(Pid = data_long_pids), overwrite = TRUE)

# Construct the SQL query to join the temporary table with DRUG_MAIN
query <- "SELECT dm.Pid, dm.Atc, dm.Age, dm.Sex, dm.startdate, dm.enddate, dm.Antal
          FROM SCHEMA_X.DATASET_Y dm
          INNER JOIN #TempPids tl ON dm.Pid = tl.Pid
          WHERE dm.Age >= 64 
            AND dm.startdate >= '2005-01-01' AND dm.startdate < '2018-01-01'
            AND dm.startdate <= dm.enddate
            AND dm.Atc LIKE 'N04BD%'"

# Execute the SQL query
mao_all <- DBI::dbGetQuery(conn, query)

# Drop the temporary SQL table
DBI::dbExecute(conn, "DROP TABLE #TempPids")

rm(data_long_pids)
gc()

# Prepare MAO data
mao_all <- mao_all %>%
  rename(mao_startdate=startdate, mao_enddate=enddate) %>% 
  dplyr::select(Pid, mao_startdate, mao_enddate)  %>% 
  mutate(mao_enddate = as.Date(mao_enddate)) %>%
  distinct(Pid, mao_startdate, .keep_all = TRUE)

mao_cov <- dplyr::select(mao_all, Pid, mao_startdate, mao_enddate)

sixmonths <- 6*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = mao_cov, main_date_col = "current_date", ref_date_col = "mao_startdate", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = Inf)

# Calculate difference in days
data_long$mao_sixmonths <- 0
data_long[difftime(data_long$current_date, data_long$mao_enddate, units = "days") < 6 * 30.5, 'mao_sixmonths'] <- 1

rm(mao_all, mao_cov)
gc()

## Charlson Comorbidity Index (CCI)

# Extract CCI data

# Myocardial_infarction
icd10 <- "\\<I21|\\<I22|\\<I252"
cci_prefix <- 'MI'
MI <- extract_event_dates(icd10, patients_SSRI_MIR, 
                          sprintf('%s_date', cci_prefix), 
                          sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE) 

# Congestive_heart_failure
icd10 <- "\\<I110|\\<I130|\\<I132|\\<I255|\\<I420|\\<I426|\\<I427|\\<I428|\\<I429|\\<I43|\\<I50"
cci_prefix <- 'CHF'
CHF <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE) 

# Peripheral_vascular_disease
icd10 <- "\\<I70|\\<I71|\\<I731|\\<I738|\\<I739|\\<I771|\\<I790|\\<I792|\\<K55"
cci_prefix <- 'PVD'
PVD <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Cerebrovascular_disease
icd10 <- "\\<G45|\\<I60|\\<I61|\\<I62|\\<I63|\\<I64|\\<I67|\\<I69"
cci_prefix <- 'CVD'
CVD <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)   

# Chronic_obstructive_pulmonary_disease
icd10 <- "\\<J43|\\<J44"
cci_prefix <- 'COPD'
COPD <- extract_event_dates(icd10, patients_SSRI_MIR, 
                            sprintf('%s_date', cci_prefix), 
                            sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Chronic_other_pulmonary_disease
icd10 <- paste(c("\\<J45",41,42,46,47,60:70),collapse="|\\<J")
cci_prefix <- 'CPD'
CPD <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Rheumatic_disease
icd10 <- paste(c("\\<M05","06",123,"070","071","072","073","08",13,30,313:316,32:34,350:351,353,45:46),collapse="|\\<M")
cci_prefix <- 'RD'
RD <- extract_event_dates(icd10, patients_SSRI_MIR, 
                          sprintf('%s_date', cci_prefix), 
                          sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Dementia
icd10 <- "\\<F00|\\<F01|\\<F02|\\<F03|\\<F051|\\<G30|\\<G311|\\<G319"
cci_prefix <- 'DEM'
DEM <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Hemiplegia
icd10 <- "\\<G114|\\<G80|\\<G81|\\<G82|\\<G830|\\<G831|\\<G832|\\<G833|\\<G838"
cci_prefix <- 'HEM'
HEM <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Diabetes_without_chronic_complication
icd10 <- "\\<E100|\\<E101|\\<E110|\\<E111|\\<E120|\\<E121|\\<E130|\\<E131|\\<E140|\\<E141"
cci_prefix <- 'DIA_WO'
DIA_WO <- extract_event_dates(icd10, patients_SSRI_MIR, 
                              sprintf('%s_date', cci_prefix), 
                              sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Diabetes_with_chronic_complication
icd10 <- "\\<E102|\\<E103|\\<E104|\\<E105|\\<E107|\\<E112|\\<E113|\\<E114|\\<E115|\\<E116|\\<E117|\\<E122|\\<E123|\\<E124|\\<E125|\\<E126|\\<E127|\\<E132|\\<E133|\\<E134|\\<E135|\\<E136|\\<E137|\\<E142|\\<E143|\\<E144|\\<E145|\\<E146|\\<E147"
cci_prefix <- 'DIA_W'
DIA_W <- extract_event_dates(icd10, patients_SSRI_MIR, 
                             sprintf('%s_date', cci_prefix), 
                             sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Renal_disease
cci_prefix <- 'REN'
REN <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Mild_liver_disease
cci_prefix <- 'MLD'
MLD <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# liver special
icd10 <- "\\<R18"
cci_prefix <- 'LIVSP'
LIVSP <- extract_event_dates(icd10, patients_SSRI_MIR, 
                             sprintf('%s_date', cci_prefix), 
                             sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# moderate severe liver disease
icd10 <-  "\\<I850|\\<I859|\\<I982|\\<I983"
cci_prefix <- 'MSLIV'
MSLIV <- extract_event_dates(icd10, patients_SSRI_MIR, 
                             sprintf('%s_date', cci_prefix), 
                             sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Peptic_ulcer_disease
icd10 <-"\\<K25|\\<K26|\\<K27|\\<K28"
cci_prefix <- 'PUD'
PUD <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Malignancy
icd10  <- paste("\\<C00|\\<C0",paste(1:9,collapse = "|\\<C0",sep=""),paste("|\\<C",paste(10:76,collapse = "|\\<C"),sep=""),paste("|\\<C",paste(81:86,collapse = "|\\<C"),sep=""),paste("|\\<C",paste(88:97,collapse = "|\\<C"),sep=""),sep="")
cci_prefix <- 'MAL'
MAL <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Metastatic_cancer
icd10 <- "\\<C77|\\<C78|\\<C79|\\<C80"
cci_prefix <- 'MET'
MET <- extract_event_dates(icd10, patients_SSRI_MIR, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Aids
icd10 <- "\\<B20|\\<B21|\\<B22|\\<B23|\\<B24|\\<F024|\\<O987|\\<R75|\\<Z114|\\<Z219|\\<Z711"
cci_prefix <- 'AIDS'
AIDS <- extract_event_dates(icd10, patients_SSRI_MIR, 
                            sprintf('%s_date', cci_prefix), 
                            sprintf("ETT/SSRI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Merge CCI data

data_long <- as.data.table(data_long)
data_long <- merge_cci(MI, 'MI_date', data_long, current_date, "MI")
data_long <- merge_cci(CHF, 'CHF_date', data_long, current_date, 'CHF')
data_long <- merge_cci(PVD, 'PVD_date', data_long, current_date, "PVD")
data_long <- merge_cci(CVD, 'CVD_date', data_long, current_date, "CVD")
data_long <- merge_cci(COPD, 'COPD_date', data_long, current_date, "COPD")
data_long <- merge_cci(CPD, 'CPD_date', data_long, current_date, "CPD")
data_long <- merge_cci(RD, 'RD_date', data_long, current_date, "RD")
data_long <- merge_cci(DEM, 'DEM_date', data_long, current_date, "DEM")
data_long <- merge_cci(HEM, 'HEM_date', data_long, current_date, "HEM")
data_long <- merge_cci(DIA_WO, 'DIA_WO_date', data_long, current_date, "DIA_WO")
data_long <- merge_cci(DIA_W, 'DIA_W_date', data_long, current_date, "DIA_W")
data_long <- merge_cci(REN, 'REN_date', data_long, current_date, "REN")
data_long <- merge_cci(MLD, 'MLD_date', data_long, current_date, "MLD")
data_long <- merge_cci(LIVSP, 'LIVSP_date', data_long, current_date, "LIVSP")
data_long <- merge_cci(MSLIV, 'MSLIV_date', data_long, current_date, "MSLIV")
data_long <- merge_cci(PUD, 'PUD_date', data_long, current_date, "PUD")
data_long <- merge_cci(MAL, 'MAL_date', data_long, current_date, "MAL")
data_long <- merge_cci(MET, 'MET_date', data_long, current_date, "MET")
data_long <- merge_cci(AIDS, 'AIDS_date', data_long, current_date, "AIDS")

data_long <- data_long %>% mutate(MSLIV=if_else(MLD==1 & LIVSP==1,1,MSLIV))
data_long <- data_long %>% mutate(MLD=if_else(MSLIV==1,0,MLD))

# Generate weighted CCI
data_long$CCIW <- data_long$MI + data_long$CHF + data_long$PVD + 
  data_long$CVD + data_long$COPD + data_long$CPD + 
  data_long$RD + data_long$DEM + 2 * data_long$HEM + data_long$DIA_WO + 
  2 * data_long$DIA_W + 2 * data_long$RD + data_long$MLD + 3 * data_long$MSLIV + 
  data_long$PUD + 2 * data_long$MAL + 6 * data_long$MET + 6 * data_long$AIDS

# CCI baseline
data_long <- data_long %>%
  group_by(Pid) %>%
  mutate(CCIW_sb = ifelse(tstart == 0, CCIW[tstart == 0], NA)) %>%
  tidyr::fill(CCIW_sb)

rm(MI, CHF, PVD, CVD, COPD, CPD, RD, DEM, HEM, DIA_WO, DIA_W, REN, MLD, LIVSP, MSLIV, PUD, MAL, MET, AIDS)
gc()
Sys.sleep(60)

### GI comorbidities

## select gastrobleeding patients and unique dates
Vh.patients  <- patients_SSRI_MIR[grep("\\<K92",patients_SSRI_MIR$DIA),]

GIbleeding <- Vh.patients %>% group_by(pid) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(pid, date_admission) %>% rename(GI_bld_date=date_admission)
rm(Vh.patients)
gc()

## select IBS patients and unique dates
Vh10  <- patients_SSRI_MIR[grep("\\<K58",patients_SSRI_MIR$DIA),]
IBS <- Vh10 %>% group_by(pid) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(pid, date_admission) %>% rename(IBS_date=date_admission)
rm(Vh10)
gc()

## select DD patients and unique dates
Vh10  <- patients_SSRI_MIR[grep("\\<K572|\\<K573|\\<K574|\\<K575|\\<K578|\\<K579",patients_SSRI_MIR$DIA),]
DD <- Vh10 %>% group_by(pid) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(pid, date_admission) %>% rename(DD_date=date_admission)
rm(Vh10)
gc()

## select diarrhea patients and unique dates
Vh59  <- patients_SSRI_MIR[grep("\\<K591",patients_SSRI_MIR$DIA),]
Diarrhea <- Vh59 %>% group_by(pid) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(pid, date_admission) %>% rename(Diarrh_date=date_admission)
rm(Vh59)
gc()

## select celiac cases based on TOPO and SNOMED
Celiac_data <- as.data.frame(OUTCOME) %>%
  mutate(Celiac_D_DATE = ifelse(Celiac_D_DATE == 'NA', NA, Celiac_D_DATE))

Celiac <- dplyr::select(Celiac_data, PID, Celiac_D_DATE) %>% collect %>% filter(!is.na(Celiac_D_DATE))
Celiac <- Celiac %>%
  mutate(Celiac_date=as.Date(Celiac_D_DATE, format = "%Y%m%d")) %>%
  dplyr::select(-Celiac_D_DATE) %>% rename(pid=PID)

rm(Celiac_data)
gc()
Sys.sleep(60)

# Add GI features
data_long <- as.data.table(data_long)
data_long <- merge_cci(GIbleeding, 'GI_bld_date', data_long, current_date, "GI_bld")
data_long <- merge_cci(IBS, 'IBS_date', data_long, current_date, "IBS")
data_long <- merge_cci(DD, 'DD_date', data_long, current_date, "DD")
data_long <- merge_cci(Diarrhea, 'Diarrh_date', data_long, current_date, "Diarrh")
data_long <- merge_cci(Celiac, 'Celiac_date', data_long, current_date, "Celiac")

# Number of gastrointestinal comorbidities

data_long$N.GI <- data_long$GI_bld + data_long$IBS + data_long$DD + 
  data_long$Diarrh + data_long$Celiac

# GI comorbidities baseline
data_long <- data_long %>%
  group_by(Pid) %>%
  mutate(N.GI_sb = ifelse(tstart == 0, N.GI[tstart == 0], NA)) %>%
  tidyr::fill(N.GI_sb)

rm(GIbleeding, IBS, DD, Diarrhea, Celiac)
gc()
Sys.sleep(60)

# Indication of use

# HTN
htn <- patients_SSRI_MIR[grep("\\<I10|\\<I11|\\<I12|\\<I13|\\<I14|\\<I15",patients_SSRI_MIR$DIA),]  %>% group_by(pid) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

dm <- patients_SSRI_MIR[grep("\\<E10|\\<E11|\\<E12|\\<E13|\\<E14",patients_SSRI_MIR$DIA),]  %>% group_by(pid) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

ihd <- patients_SSRI_MIR[grep("\\<I20|\\<I21|\\<I22|\\<I23|\\<I24|\\<I25",patients_SSRI_MIR$DIA),]  %>% group_by(pid) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

chf <- patients_SSRI_MIR[grep("\\<I110|\\<I130|\\<I132|\\<I255|\\<I420|\\<I426|\\<I427|\\<I428|\\<I429|\\<I43|\\<I50",patients_SSRI_MIR$DIA),]  %>% group_by(pid) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

pvd <- patients_SSRI_MIR[grep("\\<I70|\\<I71|\\<I731|\\<I738|\\<I739|\\<I771|\\<I790|\\<I792|\\<K55",patients_SSRI_MIR$DIA),]  %>% group_by(pid) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup 

cvd <- patients_SSRI_MIR[grep("\\<G45|\\<I60|\\<I61|\\<I62|\\<I63|\\<I64|\\<I67|\\<I69",patients_SSRI_MIR$DIA),]  %>% group_by(pid) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup 

indication_df <- rbind(htn, dm, ihd, chf, pvd, cvd) %>% group_by(pid) %>%
  filter(row_number(date_admission_c)==1) %>% 
  dplyr::select(pid, date_admission) %>% rename(indication_date=date_admission) %>%
  ungroup

data_long <- as.data.table(data_long)
data_long <- merge_cci(indication_df, 'indication_date', data_long, current_date, "ind_SSRI")

# Baseline indication
data_long <- data_long %>%
  group_by(Pid) %>%
  mutate(ind_SSRI_sb = ifelse(tstart == 0, ind_SSRI[tstart == 0], NA)) %>%
  tidyr::fill(ind_SSRI_sb)

rm(htn, dm, ihd, chf, pvd, cvd,indication_df)
gc()
Sys.sleep(60)

# Convert IBD dates into datetime format
data_long$First_IBD_DATE <- as.Date(as.character(data_long$First_IBD_DATE), format = "%Y%m%d")
data_long$Second_IBD_DATE <- as.Date(as.character(data_long$Second_IBD_DATE), format = "%Y%m%d")


## Endoscopy counts (can have multiple endoscopies per day)

## select endoscopy patients and unique dates
Vh.patients  <- patients_SSRI_MIR[grep("\\<2861|\\<2880|\\<2881|\\<4480|\\<4483|\\<4486|\\<4487|\\<4488|\\<4489|\\<4490|\\<9021|\\<4686|\\<4687|\\<9003|\\<9004|\\<9021|\\<UJC|\\<UJD|\\<UJF02|\\<UJF05|\\<9011|\\<9012|\\<9023|\\<4685|\\<4688|\\<4689|\\<4674|\\<4684|\\<UJF32|\\<UJF35|\\<UJF42|\\<UJF45",
                                      patients_SSRI_MIR$op),]

# Filter out rows with NA values in the 'datum' column, select specific columns, and rename them
endo <- Vh.patients %>% filter(!is.na(date_admission)) %>%
  dplyr::select(pid, date_admission)  %>% rename(event_date=date_admission, Pid=pid) # %>% distinct(across(all_of(c('Pid', 'event_date'))))
rm(Vh.patients)

# Create endoscopy sequence data
endo_sequence <- create_rolling_sequence(endo)

endo_sequence <- endo_sequence %>% rename(N.endo = N_rolling)

# Create a new column 'month_year' based on the 'current_date' column
data_long$month_year <- format(data_long$current_date, "%m_%Y")

# Merge endoscopy sequence data
data_long <- data_long %>% 
  left_join(endo_sequence, by=c("Pid", "month_year")) %>%
  mutate(N.endo = ifelse(is.na(N.endo), 0, N.endo))

rm(endo, endo_sequence)
gc()
Sys.sleep(60)

## Prescription counts (can have multiple prescriptions per day)

# Get the unique Pids from data_long dataframe
data_long_pids <- unique(data_long$Pid)

# Drop the temporary table if it exists
tryCatch({
  DBI::dbExecute(conn, "BEGIN TRY DROP TABLE #TempPids END TRY BEGIN CATCH END CATCH")
}, error = function(e) {})


# Create a temporary table in the database to store the Pids
DBI::dbExecute(conn, "CREATE TABLE #TempPids (Pid INT)")
DBI::dbWriteTable(conn, "#TempPids", data.frame(Pid = data_long_pids), overwrite = TRUE)

# Construct the SQL query to join the temporary table with DRUG_MAIN
query <- "SELECT dm.Pid, dm.Atc, dm.Age, dm.Sex, dm.startdate, dm.enddate, dm.Antal
          FROM SCHEMA_X.DATASET_Y dm
          INNER JOIN #TempPids tl ON dm.Pid = tl.Pid
          WHERE dm.Age >= 64
            AND dm.startdate >= '2005-01-01' AND dm.startdate < '2018-01-01'
            AND dm.startdate <= dm.enddate
          ORDER BY dm.Pid, dm.startdate"

# Execute the SQL query and store the result as drugs_N
drugs_N <- DBI::dbGetQuery(conn, query)

# Drop the temporary SQL table
DBI::dbExecute(conn, "DROP TABLE #TempPids")

rm(data_long_pids)
gc()

drugs_N <- drugs_N %>%
  filter(!is.na(startdate)) %>% distinct(Pid, startdate, Atc, .keep_all = TRUE) %>%
  dplyr::select(Pid, startdate, Atc)  %>% rename(event_date=startdate) %>%
  mutate (event_date=as.Date(event_date))

# Create prescription sequence data
drugs_sequence <- calculate_unique_atc_counts(drugs_N)

# Merge prescription sequence
data_long <- data_long %>% 
  left_join(drugs_sequence, by=c("Pid", "month_year")) %>%
  mutate(N.drugs = ifelse(is.na(N.drugs), 0, N.drugs))

rm(drugs_sequence)
gc()
Sys.sleep(60)

## Visit counts (cannot have multiple visits per day)

visits_N <- patients_SSRI_MIR %>% filter(!is.na(date_admission)) %>%
  dplyr::select(pid, date_admission)  %>% rename(event_date=date_admission, Pid=pid) %>% distinct(across(all_of(c('Pid', 'event_date'))))

# Create visits sequence data
visits_sequence <- create_rolling_sequence(visits_N)

visits_sequence <- visits_sequence %>% rename(N.visits = N_rolling)

# Merge visit sequence data
data_long <- data_long %>% 
  left_join(visits_sequence, by=c("Pid", "month_year")) %>%
  mutate(N.visits = ifelse(is.na(N.visits), 0, N.visits))

rm(visits_N, visits_sequence)
gc()
Sys.sleep(60)

# Baseline GI covariates
data_long <- data_long %>%
  group_by(Pid) %>%
  mutate(N.GI_sb = ifelse(tstart == 0, N.GI[tstart == 0], NA),
         N.endo_sb = ifelse(tstart == 0, N.endo[tstart == 0], NA),
         N.drugs_sb = ifelse(tstart == 0, N.drugs[tstart == 0], NA),
         N.visits_sb = ifelse(tstart == 0, N.visits[tstart == 0], NA)) %>%
  tidyr::fill(N.GI_sb) %>% tidyr::fill(N.endo_sb) %>% tidyr::fill(N.drugs_sb) %>% tidyr::fill(N.visits_sb)


## Encounters within the past 24, 12 and 6 months - start date of the month as reference

# Create combined encounters
combined_encounters <- rbind(drugs_N %>% dplyr::select(Pid, event_date) %>% rename(encounter_date=event_date),
                             patients_SSRI_MIR %>% dplyr::select(pid, date_admission) %>% rename(Pid=pid, encounter_date=date_admission))

combined_encounters <-combined_encounters %>% 
  distinct(Pid, encounter_date, .keep_all = TRUE) %>%
  rename(encounter_date24 = encounter_date)

# Encounter within the past 24 months
twentyfourmonths <- 24*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date", ref_date_col = "encounter_date24", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = twentyfourmonths)

combined_encounters <-combined_encounters %>% 
  rename(encounter_date12 = encounter_date24)

# Encounter within the past 12 months
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date", ref_date_col = "encounter_date12", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = twelvemonths)


combined_encounters <-combined_encounters %>% 
  rename(encounter_date6 = encounter_date12)

# Encounter within the past 6 months
sixmonths <- 6*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date", ref_date_col = "encounter_date6", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = sixmonths)


## Encounters within the past 24, 12 and 6 months - end date of the month as reference

combined_encounters <- rbind(drugs_N %>% dplyr::select(Pid, event_date) %>% rename(encounter_e_date=event_date),
                             patients_SSRI_MIR %>% dplyr::select(pid, date_admission) %>% rename(Pid=pid, encounter_e_date=date_admission))

# Create end date for each trail in a given month
data_long <-data_long %>%
  mutate(
    current_date_end = pmin(current_date + days(30), as.Date("2017-12-31"))
  )

# Encounter within the past 24 months
combined_encounters <-combined_encounters %>% 
  distinct(Pid, encounter_e_date, .keep_all = TRUE) %>%
  rename(encounter_e_date24 = encounter_e_date)

twentyfourmonths <- 24*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "encounter_e_date24", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = twentyfourmonths)

combined_encounters <-combined_encounters %>% 
  rename(encounter_e_date12 = encounter_e_date24)

# Encounter within the past 12 months
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "encounter_e_date12", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = twelvemonths)

combined_encounters <-combined_encounters %>% 
  rename(encounter_e_date6 = encounter_e_date12)

# Encounter within the past 6 months
sixmonths <- 6*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "encounter_e_date6", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = sixmonths)


## Clinic visits within the past 24, 12 and 6 months of end date of a given month to assess HU eligibility in the baseline month

combined_encounters <- patients_SSRI_MIR %>% 
  dplyr::select(pid, date_admission) %>% 
  rename(Pid=pid, visit_e_date=date_admission)

# Visit within the past 24 months
combined_encounters <-combined_encounters %>% 
  distinct(Pid, visit_e_date, .keep_all = TRUE) %>%
  rename(visit_e_date24 = visit_e_date)

twentyfourmonths <- 24*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "visit_e_date24", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = twentyfourmonths)

# Visit within the past 12 months
combined_encounters <-combined_encounters %>% 
  rename(visit_e_date12 = visit_e_date24)

twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "visit_e_date12", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = twelvemonths)

# Visit within the past 6 months
combined_encounters <-combined_encounters %>% 
  rename(visit_e_date6 = visit_e_date12)

sixmonths <- 6*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "visit_e_date6", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = sixmonths)

rm(combined_encounters)
rm(patients_SSRI_MIR, drugs_N)
gc()

## Colonoscopy
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = colonoscopy_cov, 
                           main_date_col = "current_date", ref_date_col = "colonoscopy_date", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = twelvemonths)
rm(colonoscopy_cov)


## Treatment variable
SSRI_cov <- SSRI_all %>%
  rename(med_rollstart = SSRI_startdate, med_rollend = SSRI_endate)

MIR_cov <- MIR_all %>%
  rename(med_rollstart = MIR_startdate, med_rollend = MIR_endate)

# Combine the data frames SSRI_cov and MIR_cov
med_cov <- bind_rows(SSRI_cov, MIR_cov)
rm(SSRI_cov, MIR_cov, SSRI_all, MIR_all)
gc()
Sys.sleep(60)

data_long <- merge_by_time(main_df = data_long, ref_df = med_cov, main_date_col = "current_date_end", ref_date_col = "med_rollstart", 
                           main_by_col = "Pid", ref_by_col = "Pid", roll_value = Inf)

data_long <- data_long %>%
  mutate(treat = if_else(
    is.na(med_rollstart) | is.na(med_rollend),  # Check if either med_rollstart or med_rollend is NA
    0,                                          # If true, set treat to 0
    if_else((med_rollstart <= current_date_end) & (med_rollend >= current_date), 1, 0)  # If not NA, check date conditions
  ))

rm(med_cov)
gc()
Sys.sleep(60)

## Create a binary variable indicating whether PPI was used in the past 12 months
library(parallel)
data_long <- apply_twelvemonth_lag(data_long)
gc()
Sys.sleep(60)

## Eligibility

# Applied multiple eligibility criteria to assess the effects
data_long <- data_long %>%
  mutate(
    elig6 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (med_twelvemonths==0) &
        (mao_sixmonths==0) &
        (!is.na(encounter_date6)),
      1, 0
    ),
    elig12 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (med_twelvemonths==0) &
        (mao_sixmonths==0) &
        (!is.na(encounter_date12)),
      1, 0
    ),
    elig24 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (med_twelvemonths==0) &
        (mao_sixmonths==0) &
        (!is.na(encounter_date24)),
      1, 0
    )
  )

# Death, Migration
data_long$death_date <- as.Date(as.character(data_long$deathdate), format = "%Y%m%d")
data_long$migration_date <- as.Date(as.character(data_long$migrationdate), format = "%Y%m%d")

data_long <- data_long %>%
  mutate(death = ifelse(!is.na(death_date) & death_date <= current_date_end, 1, 0))
data_long <- data_long %>%
  mutate(censor = ifelse(!is.na(migration_date) & migration_date <= current_date_end, 1, 0))

# Final list of variables
final_variables <- c(
  "Pid", "date_admission", "colonoscopy_date", "mao_sixmonths", "mao_startdate", "mao_enddate",
  "encounter_date24", "encounter_date12", "encounter_date6", "med_twelvemonths",
  "encounter_e_date24", "encounter_e_date12", "encounter_e_date6", "med_type",
  "visit_e_date24", "visit_e_date12", "visit_e_date6", 
  "med_rollstart", "med_rollend", "comp_startdate", "comp_enddate",
  "current_date", "current_date_end", "elig24", "elig12", "elig6", "treat", 
  "death", "censor", "eventMC", "timeMC", "MC_Type",
  "NMdate", "enddateNM", "eventNM", "timeNM", "death_date", "migration_date", 
  "cal_time", "cal_timesqr","tstart", "tstop", "Age", "Sex", "ind_SSRI", "ind_SSRI_sb", "indication_date",
  "N.GI", "N.endo", "N.drugs", "N.visits", "N.GI_sb", "N.endo_sb", "N.drugs_sb", "N.visits_sb", 
  "educ", "nordic", "CCIW", "CCIW_sb", "MCdate", "First_IBD_DATE", "Second_IBD_DATE",
  "migrationdate", "deathdate", "enddateMC", "ID"
)

# select the variables from the data_long dataframe
data_long <- data_long[, final_variables, with = FALSE]

write_feather(data_long, sprintf("ETT/SSRI/Processed_files/SSRI_longdata_%s.ft", event_suffix))
