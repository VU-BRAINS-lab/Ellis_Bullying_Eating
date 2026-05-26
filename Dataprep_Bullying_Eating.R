# -----------------------------------------------------------------------------
# Data preparation for ED/FBSV/RSES/PTQ manuscript
# -----------------------------------------------------------------------------
# Purpose:
#   1. Read cleaned survey data.
#   2. Keep manuscript variables.
#   3. Restrict sample to ages 17 through 23.
#   4. Report missingness separately for:
#        a. primary analysis variables, which define the analytic sample, and
#        b. demographic/item variables, which are retained but do not define the
#           analytic sample.
#   5. Create final complete-case analytic dataset based only on the primary
#      analysis variables.
#   6. Save final dataset for descriptives and mediation scripts.
# -----------------------------------------------------------------------------

# Load packages
library(dplyr)

# -----------------------------------------------------------------------------
# User-defined paths
# -----------------------------------------------------------------------------

input_file <- "~/Documents/Survey Cleaning/Survey_cleaned_2026-05-23.rds"
output_dir <- "~/Desktop/ED Paper"
output_file <- file.path(output_dir, "ED_fbsv_rses_ptq.rds")

# -----------------------------------------------------------------------------
# Read the cleaned data
# -----------------------------------------------------------------------------

survey <- readRDS(input_file)

cat("Initial survey N:", nrow(survey), "\n")
cat("Initial number of variables:", ncol(survey), "\n\n")

# Optional: print all variable names in the source data
# names(survey)

# -----------------------------------------------------------------------------
# Select variables needed for the manuscript
# -----------------------------------------------------------------------------

vars_to_keep <- c(
  "record_id", "age", "sex_factor", "race_ethnicity", "education_factor", "household_income_factor",
  "edeqs_total", paste0("edeqs_", 1:12),
  "fbsv_total", paste0("fbsv_", 1:10),
  "ptq_total", paste0("ptq_", 1:15),
  "rses_sc_score", paste0("rses_sc_", 1:10)
)

# These variables cannot be missing for the primary analyses.
analysis_vars <- c(
  "fbsv_total",
  "rses_sc_score",
  "ptq_total",
  "edeqs_total",
  "age"
)

# These variables are used for sample description only. 
# Missingness is reported, but participants are not excluded because missing.
demographic_vars <- c(
  "sex_factor",
  "race_ethnicity",
  "education_factor",
  "household_income_factor"
)

# Item-level variables are retained for scale checks/descriptives if needed.
# Missingness is reported, but participants are not excluded because individual
# items are missing, as long as the scale totals needed for analyses are present.
item_vars <- setdiff(vars_to_keep, c("record_id", analysis_vars, demographic_vars))

# Check that all requested variables exist before subsetting.
missing_vars <- setdiff(vars_to_keep, names(survey))

if (length(missing_vars) > 0) {
  stop(
    "These variables are missing from survey: ",
    paste(missing_vars, collapse = ", ")
  )
}

EDData <- survey[, vars_to_keep]

cat("N after variable selection:", nrow(EDData), "\n")
cat("Number of variables selected:", ncol(EDData), "\n\n")

# -----------------------------------------------------------------------------
# Helper function for missingness summaries
# -----------------------------------------------------------------------------

missing_summary <- function(data, vars) {
  data.frame(
    variable = vars,
    n_missing = colSums(is.na(data[, vars, drop = FALSE])),
    pct_missing = round(colMeans(is.na(data[, vars, drop = FALSE])) * 100, 2),
    row.names = NULL
  )
}

# -----------------------------------------------------------------------------
# Missingness before sample restrictions
# -----------------------------------------------------------------------------

cat("Missingness in primary analysis variables before sample restrictions:\n")
print(missing_summary(EDData, analysis_vars))
cat("\n")

cat("Missingness in demographic variables before sample restrictions, for reporting only:\n")
print(missing_summary(EDData, demographic_vars))
cat("\n")

cat("Missingness in item-level variables before sample restrictions, for reporting/scale checks only:\n")
print(missing_summary(EDData, item_vars))
cat("\n")

# -----------------------------------------------------------------------------
# Restrict age range
# -----------------------------------------------------------------------------
# Keep participants ages 17 through 23 inclusive.
# The condition age < 24 includes anyone age 23.x but excludes age 24+.
# Because age is a primary analysis variable, participants with missing age are
# removed here.

n_before_age_restriction <- nrow(EDData)

EDData <- EDData[!is.na(EDData$age) & EDData$age >= 17 & EDData$age < 24, ]

n_after_age_restriction <- nrow(EDData)
n_removed_age_restriction <- n_before_age_restriction - n_after_age_restriction

cat("N before age restriction:", n_before_age_restriction, "\n")
cat("N after age restriction:", n_after_age_restriction, "\n")
cat("Removed because age was missing or outside 17-23:", n_removed_age_restriction, "\n\n")

# Check age-bin counts after age restriction.
age_count_table_after_age_restriction <- EDData %>%
  mutate(age_bin = floor(age)) %>%
  group_by(age_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(age_bin)

cat("Age-bin counts after age restriction:\n")
print(as.data.frame(age_count_table_after_age_restriction))
cat("\n")

# -----------------------------------------------------------------------------
# Complete-case filtering for primary analysis variables only
# -----------------------------------------------------------------------------

missing_analysis_vars <- setdiff(analysis_vars, names(EDData))

if (length(missing_analysis_vars) > 0) {
  stop(
    "These analysis variables are missing from EDData: ",
    paste(missing_analysis_vars, collapse = ", ")
  )
}

n_before_cc <- nrow(EDData)

EDData <- EDData[complete.cases(EDData[, analysis_vars]), ]

n_after_cc <- nrow(EDData)
n_removed_cc <- n_before_cc - n_after_cc

cat("N before complete-case filtering on primary analysis variables:", n_before_cc, "\n")
cat("N after complete-case filtering on primary analysis variables:", n_after_cc, "\n")
cat("Removed due to missing primary analysis variables only:", n_removed_cc, "\n\n")

# -----------------------------------------------------------------------------
# Final checks
# -----------------------------------------------------------------------------

age_count_table_final <- EDData %>%
  mutate(age_bin = floor(age)) %>%
  group_by(age_bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  arrange(age_bin)

cat("Final age-bin counts:\n")
print(as.data.frame(age_count_table_final))
cat("\n")

cat("Final missingness in primary analysis variables:\n")
print(missing_summary(EDData, analysis_vars))
cat("\n")

cat("Final missingness in demographic variables, for reporting only:\n")
print(missing_summary(EDData, demographic_vars))
cat("\n")

cat("Final missingness in item-level variables, for reporting/scale checks only:\n")
print(missing_summary(EDData, item_vars))
cat("\n")

cat("Final analytic N:", nrow(EDData), "\n")
cat("Final number of variables:", ncol(EDData), "\n\n")

# -----------------------------------------------------------------------------
# Save final dataset
# -----------------------------------------------------------------------------

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(EDData, output_file)

cat("Saved final dataset to:", output_file, "\n")
