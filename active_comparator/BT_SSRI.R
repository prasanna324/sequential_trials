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

process_event_data <- function(event_name, K_values) {
  # File path based on event_name
  file_path <- paste0('ETT/SSRI/Processed_files/SSRI_bootstrap_24_', event_name, '.ft')
  
  # Read the file
  result_df <- read_feather(file_path)
  
  # Create the df for CI plot
  graph_data <- result_df %>%
    filter(metric %in% c("NonUser", "User")) %>%
    group_by(K) %>%
    summarise(
      risk0 = HR[metric == "NonUser"],
      risk1 = HR[metric == "User"],
      .groups = "drop"
    ) %>%
    mutate(
      time_0 = as.numeric(K) - 1,  # Create time column as K - 1
      risk0 = as.numeric(risk0),
      risk1 = as.numeric(risk1)
    ) %>%
    arrange(time_0)
  
  # Loop through K_values to print the HR, Lower_CI, and Upper_CI for each event_name
  for (K_val in K_values) {
    cat(sprintf('\n%s results for month %s\n', event_name, K_val - 1))
    print(result_df[result_df$K == K_val, c("HR", "Lower_CI", "Upper_CI")])
  }
  
  # Return the processed graph data for further use
  return(graph_data)
}

# Example usage:
# Define K values you are interested in
K_values <- c(13, 25)

# Process and retrieve data for each event name
graph_MC <- process_event_data("MC", K_values)
graph_CC <- process_event_data("CC", K_values)
graph_LC <- process_event_data("LC", K_values)
graph_NM <- process_event_data("NM", K_values)

#################################################################################
#################################################################################

### Construct marginal parametric cumulative incidence (risk) curves ###

create_plot <- function(event_name, K, graph_MC, HU_time) {
  if (event_name == "MC" || event_name == "LC" || event_name == "CC") {
    CI_title <- switch(event_name,
                       "MC" = "Microscopic Colitis",
                       "LC" = "Lymphocytic Colitis",
                       "CC" = "Collagenous Colitis")
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
  plot.plr <- ggplot(graph_MC, 
                     aes(x = time_0, y = risk)) + 
    geom_line(aes(y = risk1, color = "SSRI"), linewidth = 1.5) + 
    geom_line(aes(y = risk0, color = "Mirtazapine"), linewidth = 1.5) +
    xlab("Months") + 
    scale_x_continuous(limits = c(0, K), breaks = seq(0, K, by = CI_x_by)) + 
    ylab("Cumulative Incidence (%)") + 
    scale_y_continuous(limits = c(0, CI_y_max), 
                       breaks = seq(0, CI_y_max, by = CI_y_by),
                       labels = sprintf("%.2f%%", seq(0, CI_y_max, by = CI_y_by) * 100)) + 
    labs(title = CI_title) + 
    theme_minimal() + 
    theme(axis.text = element_text(size = 14), 
          legend.position = c(0.2, 0.8),
          axis.line = element_line(colour = "black"),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    font("xlab", size = 14) +
    font("ylab", size = 14) +
    font("legend.text", size = 10) +
    scale_color_manual(values = c("#E7B800", "#2E9FDF"), 
                       breaks = c('Mirtazapine', 'SSRI')) + 
    theme(plot.margin = margin(0.8, 0.8, 0.8, 0.8, "cm"))
  
  # Save the plot
  ggsave <- function(..., bg = 'white') ggplot2::ggsave(..., bg = bg)
  ggsave(paste0("ETT/SSRI/BT_Output/CI_K", K, "_HU", HU_time, "_", event_name, '_BT', ".png"), plot = plot.plr, width = 8, height = 8)
  
  # Close plot device
  while (!is.null(dev.list())) dev.off()
}

create_plot("MC", 24, graph_MC, 12)
create_plot("CC", 24, graph_CC, 12)
create_plot("LC", 24, graph_LC, 12)
create_plot("NM", 24, graph_NM, 12)


