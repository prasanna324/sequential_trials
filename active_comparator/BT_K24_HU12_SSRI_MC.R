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

med_elig <- read_feather(sprintf("ETT/SSRI/Processed_files/SSRI_med_elig_%d_%s.ft", K, event_name))

## Bootstrapping ###

# Create input list of ids (eligible persons)
med_elig_ids <- data.frame(id = unique(med_elig$id))


# Create a function to obtain risks, RD, and RR from each bootstrap sample
std.boot <- function(data, indices) {
  
  # Select individuals into each bootstrapped sample
  ids <- data$id
  boot.ids <- data.frame(id = ids[indices])
  boot.ids$bid <- 1:nrow(boot.ids)
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- left_join(boot.ids, med_elig, by = "id")
  d$bid_new <- interaction(d$bid, d$trial_num)
  
  # Fit pooled logistic model to estimate discrete hazards
  fit.pool3 <- speedglm(formula = eventMC==1 ~ treat_b + time + timesqr +
                          I(treat_b*time) +  I(treat_b*timesqr) +
                          as.factor(Age_cat) + as.factor(Year_cat) + as.factor(Sex) + as.factor(educ) + as.factor(nordic) + 
                          as.factor(N.endo_cat_sb) + as.factor(N.GI_cat_sb) + as.factor(N.visits_cat_sb) + 
                          as.factor(N.drugs_cat_sb) + as.factor(CCIW_cat_sb) + as.factor(ind_cat_sb),
                        family = binomial(link = 'logit'),
                        data = d, weights = d$sw_99)
  
  # Create dataset with all time points for each individual under each treatment level
  treat0 <- expandRows(d[which(d$time==0),], count=K, count.is.col=F)
  treat0$time <- rep(seq(0, K-1), nrow(d[which(d$time==0),]))
  treat0$timesqr <- treat0$time^2
  
  # Create "treat_b" variable under no baseline vaccination
  treat0$treat_b <- 0
  
  # Create "treat_b" variable under baseline CROWN vaccination
  treat1 <- treat0
  treat1$treat_b <- 1
  
  # Extract predicted values from pooled logistic regression model for each person-time row
  # Predicted values correspond to discrete-time hazards
  treat0$p.event0 <- predict(fit.pool3, treat0, type = "response")
  treat1$p.event1 <- predict(fit.pool3, treat1, type = "response")
  # The above creates a person-time dataset where we have predicted discrete-time hazards
  # For each person-time row in the dataset
  
  # Obtain predicted survival probabilities from discrete-time hazards
  treat0.surv <- treat0 %>% group_by(bid_new) %>% mutate(surv0 = cumprod(1 - p.event0)) %>% ungroup()
  treat1.surv <- treat1 %>% group_by(bid_new) %>% mutate(surv1 = cumprod(1 - p.event1)) %>% ungroup()
  
  # Estimate risks from survival probabilities
  # Risk = 1 - S(t)
  treat0.surv$risk0 <- 1 - treat0.surv$surv0
  treat1.surv$risk1 <- 1 - treat1.surv$surv1
  
  # Get the mean in each treatment group at each time point from 0 to 23 (24 time points in total)
  risk0 <- aggregate(treat0.surv[c("treat_b", "time", "risk0")], by=list(treat0.surv$time), FUN=mean)[c("treat_b", "time", "risk0")]
  risk1 <- aggregate(treat1.surv[c("treat_b", "time", "risk1")], by=list(treat1.surv$time), FUN=mean)[c("treat_b", "time", "risk1")]
  
  # Prepare data
  graph.pred <- merge(risk0, risk1, by=c("time"))
  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$time + 1
  zero <- data.frame(cbind(0,0,0,1,0,0))
  zero <- setNames(zero,names(graph.pred))
  graph <- rbind(zero, graph.pred)
  
  graph$rd <- graph$risk1-graph$risk0
  graph$rr <- graph$risk1/graph$risk0
  
  print('Iteration done')
  
  return(c(graph$risk0,
           graph$risk1,
           graph$rd,
           graph$rr))
}


# Run 2 bootstrap samples (2 used for illustration)
set.seed(322)
gc()

risk.results <- boot(
  data = med_elig_ids,
  statistic = std.boot,
  R = 500
)

print("Boostraps ran successfully")

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

result_df <- calculate_CI_metrics(risk.results, K)

write_feather(result_df, sprintf("ETT/SSRI/Processed_files/SSRI_bootstrap_%d_%s.ft", K, event_name))

print(result_df)
