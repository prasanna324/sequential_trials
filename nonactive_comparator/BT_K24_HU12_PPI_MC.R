# Load required packages

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

med_elig <- readRDS(sprintf("ETT/PPI/Processed_files/PPI_BT_GI_%d_%s.rds", K, event_name))

print(nrow(med_elig))

## Bootstrapping ###

# Create input list of ids (eligible persons)
med_elig_ids <- data.table(id = unique(med_elig$id))

# Create the main treat_all df from med_elig
treat_all <- med_elig[time == 0]
row_count <- nrow(treat_all)

treat_all[, treat_b := NULL]
gc()

# Generate all time points for each individual
treat_all <- treat_all[, .SD[rep(1:.N, each = K)], by = .(id)]
treat_all[, time := rep(seq(0, K-1), row_count)]
treat_all[, timesqr := time^2]
gc()
Sys.sleep(3)

treat_all_arrow <- as_arrow_table(treat_all)
rm(treat_all)
gc()
Sys.sleep(3)

# Create a function to obtain risks, RD, and RR from each bootstrap sample
std.boot <- function(data, indices) {
  
  tryCatch({
    # Select individuals into each bootstrapped sample
    ids <- data$id
    boot.ids <- data.table(id = ids[indices])
    boot.ids[, bid := .I]
    
    # Subset person-time data to individuals selected into the bootstrapped sample
    d <- med_elig[boot.ids, on = "id", nomatch = 0, allow.cartesian = TRUE]
    
    # Fit the model
    fit.pool3 <- tryCatch({
      speedglm(formula = eventMC==1 ~ treat_b + time + timesqr + period + periodsqr + Year_cat +
                 I(treat_b*time) +  I(treat_b*timesqr) +
                 as.factor(Age_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) + 
                 as.factor(N.endo_cat_b) + as.factor(N.visits_cat_b) + 
                 as.factor(N.drugs_cat_b) + as.factor(CCIW_cat_b) + as.factor(ind_cat_b),
               family = binomial(link = 'logit'),
               data = d, weights = d$sw_99, sparse=FALSE)
    }, error = function(e) {
      cat("Error in model fitting:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(fit.pool3)) {
      return(rep(NA, (K + 1) * 4))
    }
    
    # Remove the d object
    rm(d)
    gc()
    Sys.sleep(3)
    
    # Subset the precomputed treat_all for the selected individuals
    treat_all_d <- as.data.table(treat_all_arrow)[id %in% boot.ids$id]
    treat_all_d[, bid := boot.ids$bid[match(treat_all_d$id, boot.ids$id)]]
    
    # Ensure the data is sorted by id, trial_num, and time before calculating cumprod
    treat_all_d <- treat_all_d[order(id, trial_num, time)]
    
    # Generate estimates for treat0
    treat_all_d[, treat_b := 0]
    treat_all_d[, p.event0 := predict(fit.pool3, treat_all_d, type = "response")]
    treat_all_d[, surv0 := cumprod(1 - p.event0), by = .(id, trial_num)]
    treat_all_d[, risk0 := 1 - surv0]
    
    treat_all_d[, c("p.event0", "surv0") := NULL]
    gc()
    
    # Modify the dataset for treat1
    treat_all_d[, treat_b := 1]
    
    # Predict probabilities and compute for treat1
    treat_all_d[, p.event1 := predict(fit.pool3, treat_all_d, type = "response")]
    treat_all_d[, surv1 := cumprod(1 - p.event1), by = .(id, trial_num)]
    treat_all_d[, risk1 := 1 - surv1]
    
    treat_all_d[, c("p.event1", "surv1") := NULL]
    gc()
    
    # Extract and compute for treat0
    risk0 <- treat_all_d[, .(risk0 = mean(risk0)), by = time]
    
    # Extract and compute for treat1
    risk1 <- treat_all_d[, .(risk1 = mean(risk1)), by = time]
    rm(treat_all_d)
    gc()
    Sys.sleep(3)
    
    # Prepare data
    graph.pred <- merge(risk0, risk1, by = "time")
    graph.pred <- as.data.frame(graph.pred)
    
    # Edit data frame to reflect that risks are estimated at the END of each interval
    graph.pred$time_0 <- graph.pred$time + 1
    zero <- data.frame(cbind(0, 0, 0, 0))
    zero <- setNames(zero, names(graph.pred))
    graph <- rbind(zero, graph.pred)
    
    graph$rd <- graph$risk1 - graph$risk0
    graph$rr <- graph$risk1 / graph$risk0
    
    print('Iteration done')
    
    return(c(graph$risk0, graph$risk1, graph$rd, graph$rr))
  }, error = function(e) {
    cat("Error in bootstrap iteration:", conditionMessage(e), "\n")
    return(rep(NA, (K + 1) * 4))
  })
}

# Define the function for performing bootstrap iterations and saving combined results
perform_bootstrap <- function(data, statistic, total_iterations, iterations_per_save, event_name) {
  # Construct the file path
  file_path_template <- "ETT/PPI/Processed_files/PPI_RRList_%d_%s.rds"
  file_path <- sprintf(file_path_template, K, event_name)
  
  # Initialize combined_results
  if (file.exists(file_path)) {
    combined_results <- readRDS(file_path)
    saved_iterations <- combined_results$R  # Number of iterations already saved
    
    if (saved_iterations >= total_iterations) {
      print(paste("All", saved_iterations, "iterations have already been saved. No additional iterations required."))
      return(invisible())  # Exit the function if no additional iterations are needed
    }
    
    total_iterations <- total_iterations - saved_iterations  # Adjust remaining iterations
    print(paste("Retrieved", saved_iterations, "saved iterations. Continuing with", total_iterations, "additional iterations."))
  } else {
    combined_results <- list(
      t = NULL,
      t0 = NULL,
      R = 0,
      call = NULL,
      seed = NULL
    )
    print(paste("No saved iterations found. Starting with", total_iterations, "iterations."))
  }
  
  for (i in seq(1, total_iterations, by = iterations_per_save)) {
    # Perform bootstrapping
    risk.results <- boot(
      data = data,
      statistic = statistic,
      R = iterations_per_save
    )
    
    if (is.null(combined_results$t)) {
      combined_results$t <- risk.results$t
      combined_results$t0 <- risk.results$t0
      combined_results$call <- risk.results$call
      combined_results$seed <- risk.results$seed
    } else {
      combined_results$t <- rbind(combined_results$t, risk.results$t)
      combined_results$t0 <- c(combined_results$t0, risk.results$t0)
    }
    
    # Update the number of bootstrap replicates
    combined_results$R <- nrow(combined_results$t)
    
    # Save combined results externally, overwriting the existing file
    saveRDS(combined_results, file_path)
    
    # Print the current iteration number
    print(paste("Saved iteration:", combined_results$R))
  }
}


# Call the function with the appropriate parameters
set.seed(332)
perform_bootstrap(
  data = med_elig_ids,
  statistic = std.boot,
  total_iterations = 500,
  iterations_per_save = 25,
  event_name = event_name
)


# Read the final combined results from the external path
file_path_template <- "ETT/PPI/Processed_files/PPI_RRList_%d_%s.rds"
file_path <- sprintf(file_path_template, K, event_name)
final_combined_results <- readRDS(file_path)


calculate_CI_metrics <- function(risk_results, K) {
  # Initialize an empty data frame
  result_df <- data.frame(HR = numeric(0), Lower_CI = numeric(0), Upper_CI = numeric(0), metric = character(0))
  
  # Loop through indexes from 1 to (K+1)*4
  for (i in 1:((K+1)*4)) {
    # Use tryCatch to handle errors
    temp_res <- tryCatch(
      {
        # Calculate boot.ci for each index
        boot.ci(risk_results,
                conf = 0.95,
                type = "perc",
                index = i)
      },
      error = function(e) {
        # If an error occurs (e.g., all values are NA's), assign NA to the result
        list(t0 = NA, percent = c(NA, NA, NA, NA, NA))
      }
    )
    # Determine the metric based on the index
    if (i <= (K+1)) {
      metric <- "NonUser"
      k_val <- i
    } else if (i <= 2*(K+1)) {
      metric <- "User"
      k_val <- i - (K+1)
    } else if (i <= 3*(K+1)) {
      metric <- "RD"
      k_val <- i - 2*(K+1)
    } else {
      metric <- "RR"
      k_val <- i - 3*(K+1)
    }
    # Extract HR, Lower_CI, Upper_CI, metric, and K and store them in the data frame
    result_df <- rbind(result_df, c(temp_res$t0, temp_res$percent[4], temp_res$percent[5], metric, k_val))
  }
  
  # Rename the rows of the data frame
  rownames(result_df) <- 1:nrow(result_df)
  
  # Rename the columns of the data frame
  colnames(result_df) <- c("HR", "Lower_CI", "Upper_CI", "metric", "K")
  
  return(result_df)
}


result_df <- calculate_CI_metrics(final_combined_results, K)
print(result_df)
write_feather(result_df, sprintf("ETT/PPI/Processed_files/PPI_bootstrap_DS_%d_%s.ft", K, event_name))

