# Sequential Trial Emulation: Comparative Safety of Medications

This repository contains R code for emulating sequential trials using observational healthcare data. The objective is to estimate gastrointestinal outcomes associated with different medication exposures in older adults in Sweden by applying a trial emulation framework to observational data.

The repository includes one example of each of two study design approaches: non-active comparator and active comparator, selected based on the characteristics of the medications being evaluated. Additional analyses for other exposures follow the same overall structure.

- **Non-Active Comparator**: Proton Pump Inhibitor (PPI) users vs. non-users
- **Active Comparator**: Selective Serotonin Reuptake Inhibitor (SSRI) users vs. Mirtazapine users

The analysis is structured into four stages per cohort: data preparation, sequential emulation, bootstrap estimation, and post-processing. Each stage is implemented as a separate R script.

---

## Scripts Overview

### Non-Active Comparator (PPI vs. Non-Users)

- `DC_PPI_MC.R`: Creates the cohort dataset from raw data. Applies eligibility criteria, generates exposure and outcome variables, and formats longitudinal data.
- `SE_K24_HU12_PPI_MC.R`: Implements monthly sequential trial emulation. Assigns trial entry dates, eligibility, censoring rules (healthcare utilization and follow-up), and baseline covariates.
- `BT_K24_HU12_PPI_MC.R`: Runs bootstrap estimation using a pooled logistic regression model to estimate discrete hazards. Calculates marginal risks, risk differences (RD), and risk ratios (RR).
- `BT_PPI.R`: Loads bootstrapped results, summarizes them by time, and produces cumulative incidence plots.

### Active Comparator (SSRI vs. Mirtazapine)

- `DC_SSRI_MC.R`: Assembles active comparator cohort by identifying SSRI and Mirtazapine users. Applies exclusion criteria and constructs follow-up time and treatment indicators.
- `SE_K24_HU12_SSRI_MC.R`: Applies the same sequential trial emulation logic to the active comparator cohort.
- `BT_K24_HU12_SSRI_MC.R`: Performs bootstrap estimation for marginal risks and computes RD and RR using pooled logistic regression models.
- `BT_SSRI.R`: Processes results and generates cumulative incidence plots for interpretation.

---

## Methods Summary

The analysis uses a sequential trial emulation approach with monthly staggered entry to mimic a randomized controlled trial using observational data. Eligible individuals at each time point are assigned to a new trial, with follow-up for 24 months or until censoring. Pooled logistic regression is used to estimate treatment effects, accounting for baseline covariates and time.

Bootstrapping (500 iterations) is used to compute confidence intervals for marginal risk estimates, risk differences, and risk ratios.

---

## Study Reference

This work is part of the study:

**Medications and Risk of Microscopic Colitis: A Nationwide Cohort Study of Older Adults in Sweden**  
Hamed Khalili MD MPH, Emma McGee PhD, Prasanna K Challa MS, Bjorn Roelstraete PhD, Kristina Johnell PhD, Sebastian Schneeweiss MD ScD, Jonas W. Wastesson PhD, Jonas F Ludvigsson MD PhD
