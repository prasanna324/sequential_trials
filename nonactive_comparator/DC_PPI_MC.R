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
if (!require("DBI")) install.packages("DBI")
library(DBI)

setwd('~/Medications')

set.seed(333)

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
random_sample <- function(treated_data, untreated_data, sample_size, n) {
  # Add treated column to both treated and untreated data
  treated_data$treated <- TRUE
  untreated_data$treated <- FALSE
  
  # Combine treated and untreated data
  data <- rbind(treated_data, untreated_data)
  
  # Calculate the number of samples before sampling
  data_count_before <- nrow(data)
  
  # Create an empty list to store sampled datasets
  sampled_datasets <- list()
  
  # Loop to create n datasets
  for (i in 1:n) {
    # Sample data
    sampled_data <- data %>%
      sample_n(sample_size)
    
    # Add a column specifying dataset name
    sampled_data$dataset <- paste0("DS_", i)
    
    # Store the sampled dataset in the list
    sampled_datasets[[i]] <- sampled_data
    
    # Remove sampled rows from the original data to avoid overlap
    data <- anti_join(data, sampled_data, by = "Lopnr")
  }
  
  # Combine sampled datasets into a single dataframe
  combined_df <- do.call(rbind, sampled_datasets)
  
  # Calculate the number of samples after sampling
  data_count_after <- nrow(combined_df)
  treated_count_after <- sum(combined_df$treated)
  untreated_count_after <- data_count_after - treated_count_after
  
  # Print summary statistics
  cat("Before selection:\n")
  cat("Total samples:", data_count_before, "\n")
  cat("Treated individuals:", nrow(treated_data), "\n")
  cat("Untreated individuals:", nrow(untreated_data), "\n\n")
  
  cat("After selection:\n")
  cat("Total samples:", data_count_after, "\n")
  cat("Treated individuals:", treated_count_after, "\n")
  cat("Untreated individuals:", untreated_count_after, "\n\n")
  
  return(combined_df)
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

# Applies the 12-month lag function within each individual (Lopnr)
apply_twelvemonth_lag <- function(data) {
  # Split the data frame by the grouping variable
  df_list <- split(data, data[['Lopnr']])
  
  # Define the function to be applied in parallel
  parallel_function <- function(sub_df) {
    sub_df[["PPI_twelvemonths"]] <- twelvemonth_lag(sub_df$treat)
    return(sub_df)
  }

  # Apply the function in parallel
  result_list <- mclapply(df_list, parallel_function)
  
  # Combine the results back into a data frame
  result_df <- bind_rows(result_list)
  
  return(result_df)
}

##  Function to create rolling counts for endo and visits
create_rolling_sequence <- function(data_frame) {
  # Create a new column 'month_year' based on the 'event_date' column
  data_frame$month_year <- format(data_frame$event_date, "%m_%Y")
  
  # Count the occurrences of each unique combination of 'Lopnr' and 'month_year'
  counts_df <- data_frame %>%
    group_by(Lopnr, month_year) %>%
    summarise(count = n())
  
  # Create a sequence of all unique 'Lopnr' and 'event_date' combinations
  sequence_df <- expand.grid(
    Lopnr = unique(data_frame$Lopnr),
    event_date = seq(as.Date("2005-01-01"), as.Date("2017-12-01"), by = "month")
  )
  
  # Create a new column 'month_year' based on the 'event_date' column
  sequence_df$month_year <- format(sequence_df$event_date, "%m_%Y")
  
  # Arrange the data frame based on 'Lopnr' and 'event_date'
  sequence_df <- sequence_df %>% arrange(Lopnr, event_date)
  
  # Left join 'sequence_df' with 'counts_df' on 'Lopnr' and 'month_year'
  sequence_df <- sequence_df %>% 
    left_join(counts_df, by=c("Lopnr", "month_year"))
  
  # Replace NA values in the count column with 0
  sequence_df$count[is.na(sequence_df$count)] <- 0
  
  # Create a new column 'N_rolling_lag' representing the rolling sum of 'count' for the previous 12 rows within the same 'Lopnr'
  sequence_df$N_rolling_lag <- ave(sequence_df$count, sequence_df$Lopnr, 
                                   FUN = function(x) zoo::rollapply(x, width = 12, align = "right", fill = 0, sum, na.rm = TRUE))
  
  # Group 'sequence_df' by 'Lopnr' and create a lagged version of 'N_rolling' to represent counts starting from the previous month
  sequence_df <- sequence_df %>%
    group_by(Lopnr) %>%
    mutate(N_rolling = lag(N_rolling_lag, n = 1, default = NA)) %>%
    dplyr::select(Lopnr, month_year, N_rolling)
  
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
    
    # Count unique ATC codes for each Lopnr
    counts <- subset_df %>%
      group_by(Lopnr) %>%
      summarize(unique_atc_count = n_distinct(Atc))
    
    # Create a new column for event_month_year
    counts$event_month_year <- format(as.Date(month), "%m_%Y")
    
    # Append the counts to the result list
    result_list[[length(result_list) + 1]] <- counts
  }
  
  # Combine the results into a single dataframe
  final_df <- do.call(rbind, result_list)
  
  # Ensure all combinations of Lopnr and event_month_year exist
  final_df <- merge(expand.grid(unique(final_df$Lopnr), unique(final_df$event_month_year)), final_df, by.x = c("Var1", "Var2"), by.y = c("Lopnr", "event_month_year"), all.x = TRUE)
  
  # Replace missing counts with 0
  final_df[is.na(final_df)] <- 0
  
  # Rename columns
  colnames(final_df) <- c("Lopnr", "month_year", "N.drugs")
  
  # Order the dataframe by Lopnr and event_month_year
  final_df <- final_df[order(final_df$Lopnr, as.Date(final_df$month_year, format = "%m_%Y")), ]
  
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
      group_by(lopnr) %>%
      filter(row_number(date_admission) == 1) %>%
      ungroup %>%
      rename(!!event_name := date_admission) %>%
      dplyr::select(lopnr, !!event_name)
    
    return(event_data)
  }
}

# Merges CCI or event data by rolling join on Lopnr and date, and creates an indicator variable
merge_cci <- function(input_df, date_col, result_df, current_date, new_col_name) {
  input_df <- input_df %>% rename(Lopnr = lopnr) %>% distinct(across(all_of(c('Lopnr', date_col))))
  input_df <- as.data.table(input_df)
  input_df[, join_time := input_df[[date_col]]]
  result_df[, join_time := current_date]
  setkey(input_df, Lopnr, join_time)
  setkey(result_df, Lopnr, join_time)
  result_df <- input_df[result_df, roll = Inf]
  result_df[, (new_col_name) := as.integer(!is.na(result_df[[date_col]]) & !is.na(result_df$current_date))]
  
  return(result_df)
}

#############################################################
#############################################################

conn <- DBI::dbConnect(
  odbc::odbc(), 
  .connection_string = "Driver=Driver_Name;server=Server_Name;database=DB_Name;schema=Schema_Name;trusted_connection=yes;TrustServerCertificate=yes")

# MC event data
ESPR <- read_feather("Data/ESPRESSO.ft")

# Convert dates in ESPR df to datetime format
ESPR <- ESPR %>%
  mutate(
    First_IBD_DATE = as.Date(First_IBD_DATE, format = "%Y%m%d"),
    Second_IBD_DATE = as.Date(Second_IBD_DATE, format = "%Y%m%d"),
    MC_D_DATE = as.Date(MC_D_DATE, format = "%Y%m%d"),
    MC_DATE = as.Date(MC_D_DATE, format = "%Y-%m-%d"),
    NM_DATE = as.Date(NM_D_DATE, format = "%Y%m%d")
  )

espr <- dplyr::select(ESPR, LopNr, MC_D_DATE, MC_DATE, MC_SNOMED, MC_TOPO, MC_Type, NM_D_DATE, NM_DATE, NM_Type) %>% 
  collect %>% filter(!is.na(MC_D_DATE) | !is.na(NM_D_DATE))

# Patient demographic data
COHORT <- read_feather("Data/COHORT.ft")

enddates <- dplyr::select(COHORT,lopnr,birth_year, birth_month, x_death_date, x_last_imm_date, x_last_emi_date) %>% collect
enddates <- mutate(enddates, migr.date=if_else(as.Date(x_last_emi_date) > as.Date(x_last_imm_date) | 
                                                 (!is.na(x_last_emi_date) & is.na(x_last_imm_date)),x_last_emi_date,NA_character_,missing=NA_character_))

enddates <- inner_join(enddates, espr,by=c("lopnr"="LopNr"))

enddates <- enddates %>%
  mutate(NMdate = as.numeric(gsub("-", "", NM_DATE)),
         MCdate = as.numeric(gsub("-", "", MC_DATE)),
         migrationdate = as.numeric(gsub("-", "", migr.date)),
         deathdate = as.numeric(gsub("-", "", x_death_date)),
         migr.date=as.Date(migr.date),
         x_death_date = as.Date(x_death_date))


enddates <- dplyr::select(enddates, lopnr, birth_year, MCdate, NMdate,migrationdate, deathdate, MC_DATE, NM_DATE, migr.date, x_death_date)

### Medications file

### PPI cohort

# Construct the SQL query to extract all PPI records from the drugs table
query <- "WITH RankedData AS (
              SELECT dm.Lopnr, dm.Atc, dm.Age, dm.Sex, dm.Edatum, dm.enddate, dm.Antal,
                     ROW_NUMBER() OVER (PARTITION BY dm.Lopnr, dm.Edatum ORDER BY dm.enddate DESC) AS RowNum
              FROM SCHEMA_X.DATASET_Y dm
              WHERE dm.Age >= 64 
                AND dm.Edatum >= '2005-01-01' 
                AND dm.Edatum < '2018-01-01'
                AND dm.Edatum <= dm.enddate
                AND dm.Atc LIKE 'A02BC%'
          )
          SELECT Lopnr, Atc, Age, Sex, Edatum, enddate, Antal
          FROM RankedData
          WHERE RowNum = 1;"

# Execute the SQL query
PPI_cohort <- DBI::dbGetQuery(conn, query)

print(dim(PPI_cohort))

PPI_cohort <- PPI_cohort %>%
  mutate(date_admission=Edatum, PPI.begin=Edatum,PPI.end=enddate, year = year(date_admission)) %>%
  dplyr::select(Lopnr, Atc, Edatum, enddate, Age, Sex, year, date_admission, PPI.begin, PPI.end) %>%
  arrange(Lopnr, PPI.begin, PPI.end)

# Keep all medication records for treatment assignment
PPI_all <- PPI_cohort  %>%
  rename(PPI_startdate = Edatum, PPI_endate = enddate) %>% 
  dplyr::select(Lopnr, PPI_startdate, PPI_endate)

### Non-PPI cohort

# Get the unique Lopnrs from APD_cohort data frame
PPI_lopnrs <- unique(PPI_cohort$Lopnr)

# Drop the temporary table if it exists
tryCatch({
  DBI::dbExecute(conn, "BEGIN TRY DROP TABLE #TempLopnrs END TRY BEGIN CATCH END CATCH")
}, error = function(e) {})


# Create a temporary table in the database to store the Lopnrs
DBI::dbExecute(conn, "CREATE TABLE #TempLopnrs (Lopnr INT)")
DBI::dbWriteTable(conn, "#TempLopnrs", data.frame(Lopnr = PPI_lopnrs), overwrite = TRUE)


# Construct the SQL query to join the temporary table with DRUG_MAIN
query <- "WITH RankedData AS (
              SELECT dm.Lopnr, dm.Atc, dm.Age, dm.Sex, dm.Edatum, dm.enddate, dm.Antal,
                     ROW_NUMBER() OVER (PARTITION BY dm.Lopnr ORDER BY dm.Edatum) AS RowNum
              FROM SCHEMA_X.DATASET_Y dm
              LEFT JOIN #TempLopnrs tl ON tl.Lopnr = dm.Lopnr
              WHERE dm.Age >= 65
                AND dm.Edatum >= '2006-01-01' 
                AND dm.Edatum < '2018-01-01'
                AND dm.Edatum <= dm.enddate
                AND tl.Lopnr IS NULL -- Exclude rows present in #TempLopnrs
          )
          SELECT Lopnr, Atc, Age, Sex, Edatum, enddate, Antal
          FROM RankedData
          WHERE RowNum = 1
          ORDER BY Lopnr, Edatum;"

# Execute the SQL query
non_PPI_cohort <- DBI::dbGetQuery(conn, query)

# Drop the temporary SQL table
DBI::dbExecute(conn, "DROP TABLE #TempLopnrs")

rm(PPI_lopnrs)
gc()

non_PPI_cohort <- non_PPI_cohort %>%
  mutate(date_admission=Edatum, PPI.begin=NA,PPI.end=NA, year = year(date_admission)) %>%
  dplyr::select(Lopnr, Atc, Edatum, enddate, Age, Sex, year, date_admission, PPI.begin, PPI.end) %>%
  arrange(Lopnr, PPI.begin, PPI.end) %>%
  group_by(Lopnr, PPI.begin) %>%
  slice_tail(n = 1) %>%
  ungroup()

###### PPI cohort - exclusions ###### 

# Preprocess PPI cohort to have one row per patient

unique_participants <- PPI_cohort %>%
  summarise(unique_count = n_distinct(Lopnr))
unique_participants <- unique_participants$unique_count
cat('Number of participants with PPI use:', unique_participants, '\n')

PPI_cohort <- filter(PPI_cohort, year > 2005 & Age > 64)
PPI_cohort <- subset(PPI_cohort,select = -year)
exclusion_count <-  unique_participants - length(unique(PPI_cohort$Lopnr))
cat('Exclusion, PPI usage only prior to 2006 or under the age of 65, N:', exclusion_count, '\n')

PPI_cohort <- PPI_cohort %>%
  arrange(Lopnr, Edatum) %>%  # Arrange by Lopnr and Edatum
  distinct(Lopnr, .keep_all = TRUE)  # Drop duplicates, keeping the first occurrence


## Prior MC
## Prior IBD
## Prior PPI

## Create a lag variable for admission date for merging
PPI_cohort$date_admission_lag <- PPI_cohort$date_admission - 1

## PPI before 6 months

PPI_cov <- dplyr::select(PPI_all, Lopnr, PPI_startdate, PPI_endate)

PPI_cohort <- merge_by_time(main_df = PPI_cohort, ref_df = PPI_cov, main_date_col = "date_admission_lag", ref_date_col = "PPI_startdate", 
                     main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = Inf)

# Count the number of patients with PPI before 6 months of index date
missing_count <- sum(difftime(PPI_cohort$date_admission, PPI_cohort$PPI_endate, units = "days")<6*30.5, na.rm=T)

# Print the number of patients with PPI before 6 months of index date
cat('Exclusion PPI within 6 months of index date, N:', missing_count, '\n')

# Exclude patients with PPI before 6 months of index date (PPI on the same day as index date are not excluded)
PPI_cohort <- PPI_cohort[!(difftime(PPI_cohort$date_admission, PPI_cohort$PPI_endate, units = "days")<6*30.5) | is.na(PPI_endate)]
PPI_cohort <- PPI_cohort[, c('PPI_startdate', 'PPI_endate', 'join_time') := NULL]
PPI_cohort <- as_tibble(PPI_cohort)

# Create an indicator for treatment
PPI_cohort <- mutate(PPI_cohort,PPI=if_else(is.na(Atc),0,1))
dim(PPI_cohort)

# Add IBD & MC
PPI_cohort <- left_join(PPI_cohort,dplyr::select(ESPR,LopNr,First_IBD_DATE, Second_IBD_DATE, 
                                   IBD_Type,MC_D_DATE, MC_Type, NM_D_DATE),by=c("Lopnr"="LopNr"),copy=T)

# Exclude IBD prior to baseline
exclusion_count <- sum(!(as.numeric(PPI_cohort$date_admission) - as.numeric(PPI_cohort$First_IBD_DATE) < 0 | is.na(PPI_cohort$First_IBD_DATE)))
cat('Exclusion, prior IBD, N:', exclusion_count, '\n')
PPI_cohort <- filter(PPI_cohort,as.numeric(date_admission) - as.numeric(First_IBD_DATE) < 0 | is.na(First_IBD_DATE))

# Exclude MC prior to baseline
exclusion_count <- sum(!(as.numeric(PPI_cohort$date_admission) - as.numeric(PPI_cohort$MC_D_DATE) < 0 | is.na(PPI_cohort$MC_D_DATE)))
cat('Exclusion, prior MC, N:', exclusion_count, '\n')
PPI_cohort <- filter(PPI_cohort,as.numeric(date_admission) - as.numeric(MC_D_DATE) < 0 | is.na(MC_D_DATE))

# Add Death and Migration dates
PPI_cohort <- left_join(PPI_cohort, enddates, by=c("Lopnr"="lopnr"))

# Exclude Migration prior to baseline (if any)
exclusion_count <- sum(!(as.numeric(PPI_cohort$date_admission) - as.numeric(PPI_cohort$migr.date) < 0 | is.na(PPI_cohort$migr.date)))
cat('Exclusion, migration on or prior to date of admission, N:', exclusion_count, '\n')
PPI_cohort <- filter(PPI_cohort,as.numeric(date_admission) - as.numeric(migr.date) < 0 | is.na(migr.date))

# Exclude Death on index date
exclusion_count <- sum(!(as.numeric(PPI_cohort$date_admission) - as.numeric(PPI_cohort$x_death_date) < 0 | is.na(PPI_cohort$x_death_date)))
cat('Exclusion, death on or prior to date of admission, N:', exclusion_count, '\n')
PPI_cohort <- filter(PPI_cohort,as.numeric(date_admission) - as.numeric(x_death_date) < 0 | is.na(x_death_date))



###### Non-PPI cohort - exclusions ###### 

# Preprocess PPI cohort to have one row per patient

unique_participants <- non_PPI_cohort %>%
  summarise(unique_count = n_distinct(Lopnr))
unique_participants <- unique_participants$unique_count
cat('Number of participants with no PPI use:', unique_participants, '\n')

non_PPI_cohort <- filter(non_PPI_cohort, year > 2005 & Age > 64)
non_PPI_cohort <- subset(non_PPI_cohort,select = -year)
exclusion_count <-  unique_participants - length(unique(non_PPI_cohort$Lopnr))
cat('Exclusion, drugs usage only prior to 2006 or under the age of 65, N:', exclusion_count, '\n')

# Preprocess non-PPI cohort to have one row per patient

dim(non_PPI_cohort)

non_PPI_cohort <- non_PPI_cohort %>%
  arrange(Lopnr, Edatum) %>%  # Arrange by Lopnr and Edatum
  distinct(Lopnr, .keep_all = TRUE)  # Drop duplicates, keeping the first occurrence

dim(non_PPI_cohort)

## Prior MC
## Prior IBD
## Prior PPI

## Create a lag variable for admission date for merging
non_PPI_cohort$date_admission_lag <- non_PPI_cohort$date_admission - 1

## PPI before 6 months

PPI_cov <- dplyr::select(PPI_all, Lopnr, PPI_startdate, PPI_endate)

non_PPI_cohort <- merge_by_time(main_df = non_PPI_cohort, ref_df = PPI_cov, main_date_col = "date_admission_lag", ref_date_col = "PPI_startdate", 
                            main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = Inf)

# Count the number of patients with PPI before 6 months of index date
missing_count <- sum(difftime(non_PPI_cohort$date_admission, non_PPI_cohort$PPI_endate, units = "days")<6*30.5, na.rm=T)

# Print the number of patients with PPI before 6 months of index date
cat('Exclusion PPI within 6 months of index date, N:', missing_count, '\n')

# Exclude patients with PPI before 6 months of index date (PPI on the same day as index date are not excluded)
non_PPI_cohort <- non_PPI_cohort[!(difftime(non_PPI_cohort$date_admission, non_PPI_cohort$PPI_endate, units = "days")<6*30.5) | is.na(PPI_endate)]
non_PPI_cohort <- non_PPI_cohort[, c('PPI_startdate', 'PPI_endate', 'join_time') := NULL]
non_PPI_cohort <- as_tibble(non_PPI_cohort)


# Create an indicator for treatment
non_PPI_cohort <- mutate(non_PPI_cohort,PPI=if_else(is.na(Atc),0,1))
dim(non_PPI_cohort)

# Add IBD & MC
non_PPI_cohort <- left_join(non_PPI_cohort,dplyr::select(ESPR,LopNr,First_IBD_DATE, Second_IBD_DATE, 
                                                 IBD_Type,MC_D_DATE, MC_Type, NM_D_DATE),by=c("Lopnr"="LopNr"),copy=T)

# Exclude IBD prior to baseline
exclusion_count <- sum(!(as.numeric(non_PPI_cohort$date_admission) - as.numeric(non_PPI_cohort$First_IBD_DATE) < 0 | is.na(non_PPI_cohort$First_IBD_DATE)))
cat('Exclusion, prior IBD, N:', exclusion_count, '\n')

non_PPI_cohort <- filter(non_PPI_cohort,as.numeric(date_admission) - as.numeric(First_IBD_DATE) < 0 | is.na(First_IBD_DATE))

# Exclude MC prior to baseline
exclusion_count <- sum(!(as.numeric(non_PPI_cohort$date_admission) - as.numeric(non_PPI_cohort$MC_D_DATE) < 0 | is.na(non_PPI_cohort$MC_D_DATE)))
cat('Exclusion, prior MC, N:', exclusion_count, '\n')
non_PPI_cohort <- filter(non_PPI_cohort,as.numeric(date_admission) - as.numeric(MC_D_DATE) < 0 | is.na(MC_D_DATE))

## Add Death and Migration dates
non_PPI_cohort <- left_join(non_PPI_cohort, enddates, by=c("Lopnr"="lopnr"))

# Exclude Migration prior to baseline (if any)
exclusion_count <- sum(!(as.numeric(non_PPI_cohort$date_admission) - as.numeric(non_PPI_cohort$migr.date) < 0 | is.na(non_PPI_cohort$migr.date)))
cat('Exclusion, migration on or prior to date of admission, N:', exclusion_count, '\n')
non_PPI_cohort <- filter(non_PPI_cohort,as.numeric(date_admission) - as.numeric(migr.date) < 0 | is.na(migr.date))

# Exclude Death on index date
exclusion_count <- sum(!(as.numeric(non_PPI_cohort$date_admission) - as.numeric(non_PPI_cohort$x_death_date) < 0 | is.na(non_PPI_cohort$x_death_date)))
cat('Exclusion, death on or prior to date of admission, N:', exclusion_count, '\n')
non_PPI_cohort <- filter(non_PPI_cohort,as.numeric(date_admission) - as.numeric(x_death_date) < 0 | is.na(x_death_date))

####### Combine cohorts ####### 

cat('Remaining individuals in the user group, N:', nrow(PPI_cohort), '\n')
cat('Remaining individuals in the non-user group, N:', nrow(non_PPI_cohort), '\n')

# PPI_combined <- random_sample(treated_data=PPI_cohort, untreated_data=non_PPI_cohort, 
#                               sample_size = 25000, n=12)

PPI_combined <- rbind(PPI_cohort, non_PPI_cohort)

cat('Number of individuals in the combined cohort:', nrow(PPI_combined), '\n')

rm(espr, enddates)
rm(PPI_cohort, non_PPI_cohort)
gc()

####### Add covariates ####### 

## Patient encounter data
patients <- read_feather("Data/patients.ft")
patients <- patients %>%
  filter(date_admission < '2018-01-01')
patients_PPI <- inner_join(patients,dplyr::select(PPI_combined,Lopnr),by=c("lopnr"="Lopnr"))
rm(patients)
gc()


## Colonoscopy
colonoscopy_cov <- patients_PPI %>%
  filter(grepl("\\<9011|\\<9023|\\<4688|\\<4689|\\<4674|\\<4684|\\<UJF32|\\<UJF35", patients_PPI$op)) %>%
  rename(Lopnr=lopnr,colonoscopy_date = date_admission) %>% 
  dplyr::select(Lopnr, colonoscopy_date)
  
## Endoscopy
endoscopy_cov  <- patients_PPI %>%
  filter(grepl("\\<2861|\\<2880|\\<2881|\\<4480|\\<4483|\\<4486|\\<4487|\\<4488|\\<4489|\\<4490|\\<9021|\\<4686|\\<4687|\\<9003|\\<9004|\\<9021|\\<UJC|\\<UJD|\\<UJF02|\\<UJF05|\\<9011|\\<9012|\\<9023|\\<4685|\\<4688|\\<4689|\\<4674|\\<4684|\\<UJF32|\\<UJF35|\\<UJF42|\\<UJF45", patients_PPI$op)) %>%
  rename(Lopnr=lopnr,endoscopy_date = date_admission) %>% 
  dplyr::select(Lopnr, endoscopy_date)

### Education & country
educountry <- inner_join(COHORT,dplyr::select(PPI_combined,Lopnr),by=c("lopnr"="Lopnr"),copy=TRUE)
educountry <- collect(dplyr::select(educountry,lopnr,x_highest_educ, country_group4))

educountry <- mutate(educountry, educ=ifelse(x_highest_educ=="Data not available",0,
                                             ifelse(x_highest_educ=="Pre-secondary education shorter than 9 years"|x_highest_educ=="Pre-secondary education 9 years",1,
                                                    ifelse(x_highest_educ=="Postgraduate education"|x_highest_educ=="Post-secondary education shorter than 3 years",2,3))))
educountry <- mutate(educountry, nordic=ifelse(country_group4=="Norden utom Sverige och Finland"|country_group4=="Finland"|country_group4=="Sverige",1,0))
educountry <- dplyr::select(educountry,lopnr,educ,nordic)

PPI_combined <- left_join(PPI_combined, educountry, by=c("Lopnr"="lopnr"))

cat("Number of missing values in nordic (assign 1):", sum(is.na(PPI_combined$nordic)), "\n")

# Convert NAs to 1 in the 'nordic' column (if any)
PPI_combined$nordic[is.na(PPI_combined$nordic)] <- 1

rm(educountry, COHORT)
gc()

### Study end dates
PPI_combined$enddateMC <- pmin(PPI_combined$MCdate, PPI_combined$migrationdate, PPI_combined$deathdate, 20171231,na.rm=TRUE) #
PPI_combined$eventMC   <- if_else(PPI_combined$enddateMC==PPI_combined$MCdate,1,0,missing=0)
PPI_combined$enddateMC <- as.Date(paste(substr(PPI_combined$enddateMC,1,4),'-',substr(PPI_combined$enddateMC,5,6),'-',substr(PPI_combined$enddateMC,7,8),sep=''))
PPI_combined$timeMC    <- as.numeric(PPI_combined$enddateMC - as.Date(PPI_combined$date_admission))
PPI_combined$impDMC    <- if_else(PPI_combined$timeMC<=0,1,0)
PPI_combined$timeMC    <- if_else(PPI_combined$timeMC==0,1,PPI_combined$timeMC)
PPI_combined$time_yrMC <- PPI_combined$timeMC/365.25
PPI_combined <- mutate(PPI_combined, ExcludeMC=ifelse(impDMC==1, 1, 0))

# Set NM date to NA if NM occurs after MC
PPI_combined$NMdate[PPI_combined$NM_DATE > PPI_combined$MCdate] <- NA
# Set NM date to NA if NM occurs before baseline
PPI_combined$NMdate[PPI_combined$NM_DATE < PPI_combined$date_admission] <- NA

# Set NM date to NA if NM occurs after MC
PPI_combined$NMdate[PPI_combined$NMdate > PPI_combined$MCdate] <- NA
PPI_combined$enddateNM <- pmin(PPI_combined$NMdate, PPI_combined$migrationdate, PPI_combined$deathdate, 20171231,na.rm=TRUE) #
PPI_combined$eventNM   <- if_else(PPI_combined$enddateNM==PPI_combined$NMdate,1,0,missing=0)
PPI_combined$enddateNM <- as.Date(paste(substr(PPI_combined$enddateNM,1,4),'-',substr(PPI_combined$enddateNM,5,6),'-',substr(PPI_combined$enddateNM,7,8),sep=''))
PPI_combined$timeNM    <- as.numeric(PPI_combined$enddateNM - as.Date(PPI_combined$date_admission))
PPI_combined$impDNM    <- if_else(PPI_combined$timeNM<=0,1,0)
PPI_combined$timeNM    <- if_else(PPI_combined$timeNM==0,1,PPI_combined$timeNM)
PPI_combined$time_yrNM <- PPI_combined$timeNM/365.25
PPI_combined <- mutate(PPI_combined, ExcludeNM=ifelse(impDNM==1, 1, 0))

PPI_data <- dplyr::select(PPI_combined,Lopnr, date_admission, PPI, PPI.begin, PPI.end, 
                          Age, Sex, educ, nordic,
                          MCdate, First_IBD_DATE, Second_IBD_DATE, eventMC,timeMC, enddateMC, MC_Type,
                          NMdate, enddateNM, eventNM, timeNM, migrationdate, deathdate)

rm(PPI_combined)
gc()


### Convert data to long format
t_events<-seq(0,max(PPI_data$timeMC),30.5)
times<-data.frame("tevent"=t_events,"ID_t"=seq(1:length(t_events)))

## Exclude events prior to baseline (if any)
# cat('Excluding past outcome events, N:', sum(PPI_data[event_time]<=0), '\n')
# PPI_data <- PPI_data[event_time > 0, ]

chunk_size <- 50000  # Adjust the chunk size based on memory availability
data_long <- process_data_in_chunks(PPI_data, chunk_size, t_events, event_time, event_name, "ID")
rm(PPI_data)

data_long <- data_long %>%
  group_by(Lopnr) %>%
  arrange(Lopnr, date_admission) %>%
  mutate(current_date = date_admission + (row_number() - 1) * 30.5,
         cal_time=(row_number()-1), cal_timesqr=(row_number()-1)**2,
         tstart=(row_number()-1), tstop=row_number())


## Charlson Comorbidity Index (CCI)

# Extract CCI data

# Myocardial_infarction
icd10 <- "\\<I21|\\<I22|\\<I252"
cci_prefix <- 'MI'
MI <- extract_event_dates(icd10, patients_PPI, 
                                    sprintf('%s_date', cci_prefix), 
                                    sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE) 

# Congestive_heart_failure
icd10 <- "\\<I110|\\<I130|\\<I132|\\<I255|\\<I420|\\<I426|\\<I427|\\<I428|\\<I429|\\<I43|\\<I50"
cci_prefix <- 'CHF'
CHF <- extract_event_dates(icd10, patients_PPI, 
                          sprintf('%s_date', cci_prefix), 
                          sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE) 

# Peripheral_vascular_disease
icd10 <- "\\<I70|\\<I71|\\<I731|\\<I738|\\<I739|\\<I771|\\<I790|\\<I792|\\<K55"
cci_prefix <- 'PVD'
PVD <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Cerebrovascular_disease
icd10 <- "\\<G45|\\<I60|\\<I61|\\<I62|\\<I63|\\<I64|\\<I67|\\<I69"
cci_prefix <- 'CVD'
CVD <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)   

# Chronic_obstructive_pulmonary_disease
icd10 <- "\\<J43|\\<J44"
cci_prefix <- 'COPD'
COPD <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Chronic_other_pulmonary_disease
icd10 <- paste(c("\\<J45",41,42,46,47,60:70),collapse="|\\<J")
cci_prefix <- 'CPD'
CPD <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Rheumatic_disease
icd10 <- paste(c("\\<M05","06",123,"070","071","072","073","08",13,30,313:316,32:34,350:351,353,45:46),collapse="|\\<M")
cci_prefix <- 'RD'
RD <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Dementia
icd10 <- "\\<F00|\\<F01|\\<F02|\\<F03|\\<F051|\\<G30|\\<G311|\\<G319"
cci_prefix <- 'DEM'
DEM <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Hemiplegia
icd10 <- "\\<G114|\\<G80|\\<G81|\\<G82|\\<G830|\\<G831|\\<G832|\\<G833|\\<G838"
cci_prefix <- 'HEM'
HEM <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Diabetes_without_chronic_complication
icd10 <- "\\<E100|\\<E101|\\<E110|\\<E111|\\<E120|\\<E121|\\<E130|\\<E131|\\<E140|\\<E141"
cci_prefix <- 'DIA_WO'
DIA_WO <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Diabetes_with_chronic_complication
icd10 <- "\\<E102|\\<E103|\\<E104|\\<E105|\\<E107|\\<E112|\\<E113|\\<E114|\\<E115|\\<E116|\\<E117|\\<E122|\\<E123|\\<E124|\\<E125|\\<E126|\\<E127|\\<E132|\\<E133|\\<E134|\\<E135|\\<E136|\\<E137|\\<E142|\\<E143|\\<E144|\\<E145|\\<E146|\\<E147"
cci_prefix <- 'DIA_W'
DIA_W <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Renal_disease
cci_prefix <- 'REN'
REN <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Mild_liver_disease
cci_prefix <- 'MLD'
MLD <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# liver special
icd10 <- "\\<R18"
cci_prefix <- 'LIVSP'
LIVSP <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# moderate severe liver disease
icd10 <-  "\\<I850|\\<I859|\\<I982|\\<I983"
cci_prefix <- 'MSLIV'
MSLIV <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Peptic_ulcer_disease
icd10 <-"\\<K25|\\<K26|\\<K27|\\<K28"
cci_prefix <- 'PUD'
PUD <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Malignancy
icd10  <- paste("\\<C00|\\<C0",paste(1:9,collapse = "|\\<C0",sep=""),paste("|\\<C",paste(10:76,collapse = "|\\<C"),sep=""),paste("|\\<C",paste(81:86,collapse = "|\\<C"),sep=""),paste("|\\<C",paste(88:97,collapse = "|\\<C"),sep=""),sep="")
cci_prefix <- 'MAL'
MAL <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Metastatic_cancer
icd10 <- "\\<C77|\\<C78|\\<C79|\\<C80"
cci_prefix <- 'MET'
MET <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  

# Aids
icd10 <- "\\<B20|\\<B21|\\<B22|\\<B23|\\<B24|\\<F024|\\<O987|\\<R75|\\<Z114|\\<Z219|\\<Z711"
cci_prefix <- 'AIDS'
AIDS <- extract_event_dates(icd10, patients_PPI, 
                           sprintf('%s_date', cci_prefix), 
                           sprintf("ETT/PPI/Processed_files/ptnts_%s.ft", cci_prefix), read_data = FALSE)  


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
  group_by(Lopnr) %>%
  mutate(CCIW_sb = ifelse(tstart == 0, CCIW[tstart == 0], NA)) %>%
  tidyr::fill(CCIW_sb)

rm(MI, CHF, PVD, CVD, COPD, CPD, RD, DEM, HEM, DIA_WO, DIA_W, REN, MLD, LIVSP, MSLIV, PUD, MAL, MET, AIDS)
gc()

### GI comorbidities

## select gastrobleeding patients and unique dates
Vh.patients  <- patients_PPI[grep("\\<K92",patients_PPI$DIA),]

GIbleeding <- Vh.patients %>% group_by(lopnr) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(lopnr, date_admission) %>% rename(GI_bld_date=date_admission)

GIbleeding_all <- Vh.patients %>%
  dplyr::select(lopnr, date_admission) %>%
  distinct(lopnr, date_admission, .keep_all = TRUE) %>%
  rename(GI_encounter=date_admission)

rm(Vh.patients)
gc()

## select IBS patients and unique dates
Vh10  <- patients_PPI[grep("\\<K58",patients_PPI$DIA),]

IBS <- Vh10 %>% group_by(lopnr) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(lopnr, date_admission) %>% rename(IBS_date=date_admission)

IBS_all <- Vh10 %>% group_by(lopnr) %>%
  dplyr::select(lopnr, date_admission) %>% 
  distinct(lopnr, date_admission, .keep_all = TRUE) %>% 
  rename(GI_encounter=date_admission)

rm(Vh10)
gc()

## select DD patients and unique dates
Vh10  <- patients_PPI[grep("\\<K572|\\<K573|\\<K574|\\<K575|\\<K578|\\<K579",patients_PPI$DIA),]

DD <- Vh10 %>% group_by(lopnr) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(lopnr, date_admission) %>% rename(DD_date=date_admission)

DD_all <- Vh10 %>% group_by(lopnr) %>% 
  dplyr::select(lopnr, date_admission) %>% 
  distinct(lopnr, date_admission, .keep_all = TRUE) %>% 
  rename(GI_encounter=date_admission)

rm(Vh10)
gc()

## select diarrhea patients and unique dates
Vh59  <- patients_PPI[grep("\\<K591",patients_PPI$DIA),]
Diarrhea <- Vh59 %>% group_by(lopnr) %>% 
  filter(row_number(date_admission_c)==1) %>% ungroup() %>% 
  dplyr::select(lopnr, date_admission) %>% rename(Diarrh_date=date_admission)

Diarrhea_all <- Vh59 %>% group_by(lopnr) %>% 
  dplyr::select(lopnr, date_admission) %>% 
  distinct(lopnr, date_admission, .keep_all = TRUE) %>% 
  rename(GI_encounter=date_admission)

rm(Vh59)
gc()

## select celiac cases based on TOPO and SNOMED
Celiac_ESPR <- as.data.frame(ESPR) %>%
  mutate(Celiac_D_DATE = ifelse(Celiac_D_DATE == 'NA', NA, Celiac_D_DATE))

Celiac <- dplyr::select(Celiac_ESPR, LopNr, Celiac_D_DATE) %>% collect %>% filter(!is.na(Celiac_D_DATE)) %>%
  mutate(Celiac_date=as.Date(Celiac_D_DATE, format = "%Y%m%d")) %>%
  dplyr::select(-Celiac_D_DATE) %>% rename(lopnr=LopNr)

Celiac_all <- dplyr::select(Celiac_ESPR, LopNr, Celiac_D_DATE) %>% collect %>% filter(!is.na(Celiac_D_DATE)) %>%
  mutate(Celiac_date=as.Date(Celiac_D_DATE, format = "%Y%m%d")) %>%
  dplyr::select(-Celiac_D_DATE) %>% rename(lopnr=LopNr) %>%
  distinct(lopnr, Celiac_date, .keep_all = TRUE) %>%
  rename(GI_encounter=Celiac_date)

rm(ESPR, Celiac_ESPR)
gc()

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
  group_by(Lopnr) %>%
  mutate(N.GI_b = ifelse(tstart == 0, N.GI[tstart == 0], NA)) %>%
  tidyr::fill(N.GI_b)

rm(GIbleeding, IBS, DD, Diarrhea, Celiac)
gc()

## GI encounters
gi_encounters <- rbind(GIbleeding_all, IBS_all, DD_all, Diarrhea_all, Celiac_all) %>%
  dplyr::select(lopnr, GI_encounter) %>% distinct(lopnr, GI_encounter, .keep_all = TRUE) %>%
  rename(Lopnr=lopnr)

rm(GIbleeding_all, IBS_all, DD_all, Diarrhea_all, Celiac_all)
gc()

gi_encounters <- gi_encounters %>% 
  rename(GI_encounter12 = GI_encounter)

# GI encounter within the past 12 months
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = gi_encounters, 
                           main_date_col = "current_date", ref_date_col = "GI_encounter12", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twelvemonths)


gi_encounters <- gi_encounters %>% 
  rename(GI_encounter24 = GI_encounter12)

# GI encounter within the past 24 months
twentyfourmonths <- 24*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = gi_encounters, 
                           main_date_col = "current_date", ref_date_col = "GI_encounter24", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twentyfourmonths)

##############################################################

# Count number of unique ids before any exclusions
unique_ids_before <- length(unique(data_long$Lopnr))

# Filter data_long to include only those patients in gi_encounters who have at least one non-NA GI_encounter date
# Note: Specific to PPI cohort
data_long <- data_long %>%
  filter(Lopnr %in% gi_encounters$Lopnr[!is.na(gi_encounters$GI_encounter24)])

# Count number of unique ids after exclusion
unique_ids_after <- length(unique(data_long$Lopnr))

# Print the number of unique ids before and after exclusion
cat("Number of unique ids before exclusion:", unique_ids_before, "\n")
cat("Number of unique ids after exclusion:", unique_ids_after, "\n")

rm(gi_encounters)
gc()

##############################################################

## Indication of medication use for PPI

# GERD
gerd <- patients_PPI[grep("\\<K21",patients_PPI$DIA),]  %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

# Ulcers
ulcer <- patients_PPI[grep("\\<K25|\\<K26|\\<K27|\\<K28",patients_PPI$DIA),]  %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

# Gastritis and duodenitis
K29 <- patients_PPI[grep("\\<K29",patients_PPI$DIA),]  %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

# Barrett esophagus
barrett <- patients_PPI[grep("\\<K227",patients_PPI$DIA),]  %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

# Hematemesis
hema <- patients_PPI[grep("\\<K920",patients_PPI$DIA),]  %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

# Melena
melena <- patients_PPI[grep("\\<K921",patients_PPI$DIA),]  %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

# Gastrointestinal haemorrhage, unspecified
unsp <- patients_PPI[grep("\\<K922",patients_PPI$DIA),]  %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% ungroup

indication_df <- rbind(gerd,ulcer,K29,barrett,hema,melena,unsp) %>% group_by(lopnr) %>%
  filter(row_number(date_admission_c)==1) %>% 
  dplyr::select(lopnr, date_admission) %>% rename(indication_date=date_admission) %>%
  ungroup

data_long <- as.data.table(data_long)
data_long <- merge_cci(indication_df, 'indication_date', data_long, current_date, "ind_PPI")

# Baseline indication
data_long <- data_long %>%
  group_by(Lopnr) %>%
  mutate(ind_PPI_sb = ifelse(tstart == 0, ind_PPI[tstart == 0], NA)) %>%
  tidyr::fill(ind_PPI_sb)

rm(gerd,ulcer,K29,barrett,hema,melena,unsp,indication_df)
gc()

# Convert IBD dates into datetime format
data_long$First_IBD_DATE <- as.Date(as.character(data_long$First_IBD_DATE), format = "%Y%m%d")
data_long$Second_IBD_DATE <- as.Date(as.character(data_long$Second_IBD_DATE), format = "%Y%m%d")


## Endoscopy counts (can have multiple endoscopies per day)

## select endoscopy patients and unique dates
Vh.patients  <- patients_PPI[grep("\\<2861|\\<2880|\\<2881|\\<4480|\\<4483|\\<4486|\\<4487|\\<4488|\\<4489|\\<4490|\\<9021|\\<4686|\\<4687|\\<9003|\\<9004|\\<9021|\\<UJC|\\<UJD|\\<UJF02|\\<UJF05|\\<9011|\\<9012|\\<9023|\\<4685|\\<4688|\\<4689|\\<4674|\\<4684|\\<UJF32|\\<UJF35|\\<UJF42|\\<UJF45",
                                  patients_PPI$op),]

# Filter out rows with NA values in the 'datum' column, select specific columns, and rename them
endo <- Vh.patients %>% filter(!is.na(date_admission)) %>%
  dplyr::select(lopnr, date_admission)  %>% rename(event_date=date_admission, Lopnr=lopnr) # %>% distinct(across(all_of(c('Lopnr', 'event_date'))))
rm(Vh.patients)

# Create endoscopy sequence data
endo_sequence <- create_rolling_sequence(endo)

endo_sequence <- endo_sequence %>% rename(N.endo = N_rolling)

# Create a new column 'month_year' based on the 'current_date' column
data_long$month_year <- format(data_long$current_date, "%m_%Y")

# Merge endoscopy sequence data
data_long <- data_long %>% 
  left_join(endo_sequence, by=c("Lopnr", "month_year")) %>%
  mutate(N.endo = ifelse(is.na(N.endo), 0, N.endo))

rm(endo, endo_sequence)
gc()

## Prescription counts (can have multiple prescriptions per day)

# Get the unique Lopnrs from data_long dataframe
data_long_lopnrs <- unique(data_long$Lopnr)

# Drop the temporary table if it exists
tryCatch({
  DBI::dbExecute(conn, "BEGIN TRY DROP TABLE #TempLopnrs END TRY BEGIN CATCH END CATCH")
}, error = function(e) {})


# Create a temporary table in the database to store the Lopnrs
DBI::dbExecute(conn, "CREATE TABLE #TempLopnrs (Lopnr INT)")
DBI::dbWriteTable(conn, "#TempLopnrs", data.frame(Lopnr = data_long_lopnrs), overwrite = TRUE)

# Construct the SQL query to join the temporary table with DRUG_MAIN
query <- "SELECT dm.Lopnr, dm.Atc, dm.Age, dm.Sex, dm.Edatum, dm.enddate, dm.Antal
          FROM SCHEMA_X.DATASET_Y dm
          INNER JOIN #TempLopnrs tl ON dm.Lopnr = tl.Lopnr
          WHERE dm.Age >= 64
            AND dm.Edatum >= '2005-01-01' AND dm.Edatum < '2018-01-01'
            AND dm.Edatum <= dm.enddate
          ORDER BY dm.Lopnr, dm.Edatum"

# Execute the SQL query and store the result as drugs_N
drugs_N <- DBI::dbGetQuery(conn, query)

# Drop the temporary SQL table
DBI::dbExecute(conn, "DROP TABLE #TempLopnrs")

rm(data_long_lopnrs)
gc()

drugs_N <- drugs_N %>%
  filter(!is.na(Edatum)) %>% distinct(Lopnr, Edatum, Atc, .keep_all = TRUE) %>%
  dplyr::select(Lopnr, Edatum, Atc)  %>% rename(event_date=Edatum) %>%
  mutate (event_date=as.Date(event_date))

# Create prescription sequence data
drugs_sequence <- calculate_unique_atc_counts(drugs_N)

# Merge prescription sequence
data_long <- data_long %>% 
  left_join(drugs_sequence, by=c("Lopnr", "month_year")) %>%
  mutate(N.drugs = ifelse(is.na(N.drugs), 0, N.drugs))

rm(drugs_sequence)
gc()

## Visit counts (cannot have multiple visits per day)

visits_N <- patients_PPI %>% filter(!is.na(date_admission)) %>%
  dplyr::select(lopnr, date_admission)  %>% rename(event_date=date_admission, Lopnr=lopnr) %>% distinct(across(all_of(c('Lopnr', 'event_date'))))

# Create visits sequence data
visits_sequence <- create_rolling_sequence(visits_N)

visits_sequence <- visits_sequence %>% rename(N.visits = N_rolling)

# Merge visit sequence data
data_long <- data_long %>% 
  left_join(visits_sequence, by=c("Lopnr", "month_year")) %>%
  mutate(N.visits = ifelse(is.na(N.visits), 0, N.visits))

rm(visits_N, visits_sequence)
gc()

## Baseline GI covariates
data_long <- data_long %>%
  group_by(Lopnr) %>%
  mutate(N.GI_sb = ifelse(tstart == 0, N.GI[tstart == 0], NA),
         N.endo_sb = ifelse(tstart == 0, N.endo[tstart == 0], NA),
         N.drugs_sb = ifelse(tstart == 0, N.drugs[tstart == 0], NA),
         N.visits_sb = ifelse(tstart == 0, N.visits[tstart == 0], NA)) %>%
  tidyr::fill(N.GI_sb) %>% tidyr::fill(N.endo_sb) %>% tidyr::fill(N.drugs_sb) %>% tidyr::fill(N.visits_sb)

## Encounters within the past 24, 12 and 6 months - start date of the month as reference

# Create combined encounters
combined_encounters <- rbind(drugs_N %>% dplyr::select(Lopnr, event_date) %>% rename(encounter_date=event_date),
                             patients_PPI %>% dplyr::select(lopnr, date_admission) %>% rename(Lopnr=lopnr, encounter_date=date_admission))

combined_encounters <-combined_encounters %>% 
  distinct(Lopnr, encounter_date, .keep_all = TRUE) %>%
  rename(encounter_date24 = encounter_date)

# Encounter within the past 24 months
twentyfourmonths <- 24*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date", ref_date_col = "encounter_date24", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twentyfourmonths)

combined_encounters <-combined_encounters %>% 
  rename(encounter_date12 = encounter_date24)

# Encounter within the past 12 months
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date", ref_date_col = "encounter_date12", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twelvemonths)


combined_encounters <-combined_encounters %>% 
  rename(encounter_date6 = encounter_date12)

# Encounter within the past 6 months
sixmonths <- 6*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date", ref_date_col = "encounter_date6", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = sixmonths)

## Encounters within the past 24, 12 and 6 months - end date of the month as reference

combined_encounters <- rbind(drugs_N %>% dplyr::select(Lopnr, event_date) %>% rename(encounter_e_date=event_date),
                             patients_PPI %>% dplyr::select(lopnr, date_admission) %>% rename(Lopnr=lopnr, encounter_e_date=date_admission))

# Create end date for each trail in a given month
data_long <-data_long %>%
  mutate(
    current_date_end = pmin(current_date + days(30), as.Date("2017-12-31"))
  )

# Encounter within the past 24 months
combined_encounters <-combined_encounters %>% 
  distinct(Lopnr, encounter_e_date, .keep_all = TRUE) %>%
  rename(encounter_e_date24 = encounter_e_date)

twentyfourmonths <- 24*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "encounter_e_date24", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twentyfourmonths)

combined_encounters <-combined_encounters %>% 
  rename(encounter_e_date12 = encounter_e_date24)

# Encounter within the past 12 months
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "encounter_e_date12", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twelvemonths)

combined_encounters <-combined_encounters %>% 
  rename(encounter_e_date6 = encounter_e_date12)

# Encounter within the past 6 months
sixmonths <- 6*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "encounter_e_date6", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = sixmonths)

rm(combined_encounters)

## Clinic visits within the past 24, 12 and 6 months of end date of a given month to assess HU eligibility in the baseline month

combined_encounters <- patients_PPI %>% 
  dplyr::select(lopnr, date_admission) %>% 
  rename(Lopnr=lopnr, visit_e_date=date_admission)

# Visit within the past 24 months
combined_encounters <-combined_encounters %>% 
  distinct(Lopnr, visit_e_date, .keep_all = TRUE) %>%
  rename(visit_e_date24 = visit_e_date)

twentyfourmonths <- 24*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "visit_e_date24", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twentyfourmonths)

combined_encounters <-combined_encounters %>% 
  rename(visit_e_date12 = visit_e_date24)

# Visit within the past 12 months
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "visit_e_date12", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twelvemonths)

combined_encounters <-combined_encounters %>% 
  rename(visit_e_date6 = visit_e_date12)

# Visit within the past 6 months
sixmonths <- 6*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = combined_encounters, 
                           main_date_col = "current_date_end", ref_date_col = "visit_e_date6", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = sixmonths)

rm(combined_encounters)
rm(patients_PPI, drugs_N)
gc()


## Colonoscopy
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = colonoscopy_cov, 
                           main_date_col = "current_date_end", ref_date_col = "colonoscopy_date", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twelvemonths)

## Endoscopy
twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = endoscopy_cov, 
                           main_date_col = "current_date_end", ref_date_col = "endoscopy_date", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twelvemonths)


## Either colonoscopy or endoscopy
colonoscopy_cov  <- colonoscopy_cov %>%
  rename(scope_date = colonoscopy_date) %>% 
  dplyr::select(Lopnr, scope_date)

endoscopy_cov  <- endoscopy_cov %>%
  rename(scope_date = endoscopy_date) %>% 
  dplyr::select(Lopnr, scope_date)

scope_cov <- rbind(colonoscopy_cov, endoscopy_cov) %>% distinct(across(all_of(c('Lopnr', 'scope_date'))))

twelvemonths <- 12*30.5
data_long <- merge_by_time(main_df = data_long, ref_df = scope_cov, 
                           main_date_col = "current_date_end", ref_date_col = "scope_date", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = twelvemonths)

rm(scope_cov, colonoscopy_cov, endoscopy_cov)
gc()

## Treatment variable - PPI
PPI_cov <- PPI_all %>%
  rename(med_rollstart = PPI_startdate, med_rollend = PPI_endate)

data_long <- merge_by_time(main_df = data_long, ref_df = PPI_cov, main_date_col = "current_date_end", ref_date_col = "med_rollstart", 
                           main_by_col = "Lopnr", ref_by_col = "Lopnr", roll_value = Inf)

data_long <- data_long %>%
  mutate(treat = if_else(
    is.na(med_rollstart) | is.na(med_rollend),  # Check if either med_rollstart or med_rollend is NA
    0,                                          # If true, set treat to 0
    if_else((med_rollstart <= current_date_end) & (med_rollend >= current_date), 1, 0)  # If not NA, check date conditions
  ))

## Create a binary variable indicating whether PPI was used in the past 12 months
library(parallel)
data_long <- apply_twelvemonth_lag(data_long)
gc()

## Eligibility

# Applied multiple eligibility criteria to assess the effects
data_long <- data_long %>%
  mutate(
    elig6 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (PPI_twelvemonths==0) &
        (!is.na(encounter_date6)),
      1, 0
    ),
    elig12 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (PPI_twelvemonths==0) &
        (!(N.visits == 0 & N.drugs == 0)) &
        (!is.na(encounter_date12)),
      1, 0
    ),
    elig24 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (PPI_twelvemonths==0) &
        (!(N.visits == 0 & N.drugs == 0)) &
        (!is.na(encounter_date24)),
      1, 0
    ),
    elig_gi12 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (PPI_twelvemonths==0) &
        (!(N.visits == 0 & N.drugs == 0)) &
        (!is.na(encounter_date12)) &
        (!is.na(GI_encounter12)),
      1, 0
    ),
    elig_gi24 = ifelse(
      (current_date < get(event_enddate)) &
        (current_date < First_IBD_DATE | is.na(First_IBD_DATE)) &
        (PPI_twelvemonths==0) &
        (!(N.visits == 0 & N.drugs == 0)) &
        (!is.na(encounter_date12)) &
        (!is.na(GI_encounter24)),
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
  "Lopnr", "date_admission", "colonoscopy_date", "endoscopy_date", "scope_date", 
  "encounter_date24", "encounter_date12", "encounter_date6",
  "PPI_twelvemonths", "encounter_e_date24", "encounter_e_date12", "encounter_e_date6", 
  "visit_e_date24", "visit_e_date12", "visit_e_date6", 
  "GI_encounter24", "GI_encounter12", "elig_gi24", "elig_gi12",
  "med_rollstart", "med_rollend", "current_date", "current_date_end", "elig24", "elig12", "elig6", "treat", 
  "death", "censor", "eventMC", "timeMC", "enddateMC", "MC_Type",
  "NMdate", "enddateNM", "eventNM", "timeNM", "death_date", "migration_date", 
  "cal_time", "cal_timesqr","tstart", "tstop", "Age", "Sex", "ind_PPI", "ind_PPI_sb", "indication_date",
  "N.GI", "N.endo", "N.drugs", "N.visits", "N.GI_sb", "N.endo_sb", "N.drugs_sb", "N.visits_sb", 
  "educ", "nordic", "CCIW", "CCIW_sb", "MCdate", "First_IBD_DATE", "Second_IBD_DATE",
  "migrationdate", "deathdate", "ID"
)

# select the variables from the data_long dataframe
data_long <- data_long[, final_variables, with = FALSE]

write_feather(data_long, sprintf("ETT/PPI/Processed_files/PPI_longdata_GI_%s.ft", event_suffix))

