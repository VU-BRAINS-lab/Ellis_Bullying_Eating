# Load packages
library(psych)
library(dplyr)

# Read the cleaned data
survey <- readRDS("~/Documents/Survey Cleaning/Survey_cleaned_2025-09-23.rds")

# See what variables you have and that the sample size is what you expect
names(survey)
nrow(survey)

# Trim to just variables needed
EDData <- survey[c(grep("record_id|^age$|^sex$|^education$|race_ethnicity|household_income|edeqs_total|edeqs_1|edeqs_2|edeqs_3|edeqs_4|edeqs_5|edeqs_6|edeqs_7|edeqs_8|edeqs_9|edeqs_10|edeqs_11|edeqs_12|
                        |fbsv_1|fbsv_2|fbsv_3|fbsv_4|fbsv_5|fbsv_6|fbsv_7|fbsv_8|fbsv_9|fbsv_10|fbsv_total|ptq_1|ptq_2|ptq_3|ptq_4|ptq_5|ptq_6|ptq_7|ptq_8|ptq_9|ptq_10
                        |ptq_11|ptq_12|ptq_13|ptq_14|ptq_15|ptq_total|rses_sc_1|rses_sc_2|rses_sc_3|rses_sc_4|rses_sc_5|rses_sc_6|rses_sc_7|rses_sc_8|rses_sc_9|rses_sc_10|rses_sc_score",names(survey)))]

# Check that it worked; the number of subjects should not have changed
names(EDData)
nrow(EDData)

# Trim to those 17-23 years of age and remove those with NAs for age
EDData <- EDData[!is.na(EDData$age) & EDData$age >= 17 & EDData$age < 24, ]
nrow(EDData)

## Check that you have the expected numbers in each age bin
# Create integer age bins by flooring the age values
data_clean <- EDData %>%
  filter(!is.na(age)) %>%
  mutate(age_bin = floor(age))
# Create a summary table of counts for each age
age_count_table <- data_clean %>%
  group_by(age_bin) %>%
  summarise(count = n()) %>%
  arrange(age_bin)
# Display the formatted table
print(as.data.frame(age_count_table))

# List of variables of interest
vars_of_interest <- c(
  "age",
  "sex",
  "fbsv_total",
  "ptq_total",
  "edeqs_total",
  "rses_sc_score"
)

# Keep only rows with no missing values in these columns
EDData <- EDData[complete.cases(EDData[, vars_of_interest]), ]
nrow(EDData)

# Based on visual inspection, log transform the FBSV (adding 1 to all scores is needed so we don't try and take the log of 0)
EDData <- EDData %>%
  mutate(fbsv_log = log(fbsv_total + 1))

# Only bullying has a relationship with age - regress out age from bullying and use the new variable in further analyses
mod_fbsv_res <- lm(fbsv_log ~ age, data = EDData)
EDData$fbsv_res <- residuals(mod_fbsv_res)

# Final sample size
nrow(EDData)

## Count the frequency in your final age bins to make sure you still have enough in each bin
# Create integer age bins by flooring the age values
data_final <- EDData %>%
  filter(!is.na(age)) %>%
  mutate(age_bin = floor(age))
# Create a summary table of counts for each age
age_count_table <- data_final %>%
  group_by(age_bin) %>%
  summarise(count = n()) %>%
  arrange(age_bin)
# Display the formatted table
print(as.data.frame(age_count_table))

# Look at your data to make sure everything looks like you would expect
View(EDData)

# Save final dataset
saveRDS(EDData,"~/Desktop/ED Paper/ED_fbsv_rses.rds")
