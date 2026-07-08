# -----------------------------------------------------------------------------
# Descriptives for ED/FBSV/RSES/PTQ manuscript
# -----------------------------------------------------------------------------
# Purpose:
#   Descriptive statistics, correlations, diagnostics, and
#   reliability estimates.
#
# Expected input:
#   ~/Desktop/ED Paper/ED_fbsv_rses_ptq.rds
#
# Expected output:
#   ~/Desktop/ED Paper/Descriptives Outputs/
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------

required_packages <- c("dplyr", "Hmisc", "moments", "EnvStats", "car", "psych")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install the following required package(s) before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(Hmisc)
library(moments)
library(EnvStats)
library(car)
library(psych)

# -----------------------------------------------------------------------------
# File paths
# -----------------------------------------------------------------------------

input_file <- "~/Desktop/ED Paper/ED_fbsv_rses_ptq.rds"
output_dir <- "~/Desktop/ED Paper/Descriptives Outputs"

if (!file.exists(path.expand(input_file))) {
  stop("Input file not found: ", input_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Read in data
# -----------------------------------------------------------------------------

EDData <- readRDS(input_file)

cat("Analytic N:", nrow(EDData), "\n\n")

# -----------------------------------------------------------------------------
# Variable definitions
# -----------------------------------------------------------------------------

analysis_vars <- c(
  "age",
  "edeqs_total",
  "fbsv_total",
  "ptq_total",
  "rses_sc_score"
)

model_vars <- c(
  "edeqs_total",
  "fbsv_total",
  "ptq_total",
  "rses_sc_score"
)

demographic_vars <- c(
  "sex_factor",
  "race_ethnicity",
  "household_income_factor",
  "education_factor"
)

all_needed_vars <- c(analysis_vars, demographic_vars)
missing_needed_vars <- setdiff(all_needed_vars, names(EDData))

if (length(missing_needed_vars) > 0) {
  stop(
    "These expected variables are missing from EDData: ",
    paste(missing_needed_vars, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

round2 <- function(x) {
  round(x, 2)
}

format_num <- function(x, digits = 2, omit_leading_zero = FALSE) {
  ifelse(
    is.na(x),
    NA_character_,
    {
      out <- formatC(x, format = "f", digits = digits)
      if (omit_leading_zero) {
        out <- sub("^-0\\.", "-.", out)
        out <- sub("^0\\.", ".", out)
      }
      out
    }
  )
}

p_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < .001 ~ "***",
    p < .01 ~ "**",
    p < .05 ~ "*",
    TRUE ~ ""
  )
}

format_category <- function(x) {
  x_chr <- as.character(x)
  x_chr <- trimws(x_chr)
  x_chr
}

format_corr <- function(r, p) {
  if (is.na(r)) {
    return(NA_character_)
  }
  paste0(format_num(r, digits = 2, omit_leading_zero = TRUE), p_stars(p))
}

make_frequency_table <- function(data, var, label) {
  var_name <- rlang::as_name(rlang::ensym(var))
  x <- data[[var_name]]
  x_label <- format_category(x)
  missing_flag <- is.na(x_label) | x_label == ""
  nonmissing_n <- sum(!missing_flag)
  missing_n <- sum(missing_flag)

  if (is.factor(x)) {
    category_order <- format_category(levels(x))
    category_order <- category_order[!is.na(category_order) & category_order != ""]
  } else {
    category_order <- unique(x_label[!missing_flag])
  }

  out <- data.frame(.label = x_label, stringsAsFactors = FALSE) %>%
    filter(!is.na(.data$.label), .data$.label != "") %>%
    mutate(.label = factor(.data$.label, levels = unique(category_order))) %>%
    count(.data$.label, name = "n", .drop = TRUE) %>%
    mutate(
      .label = as.character(.data$.label),
      Percent = 100 * .data$n / nonmissing_n,
      denominator = nonmissing_n,
      n_missing = missing_n
    ) %>%
    transmute(
      !!label := .data$.label,
      N = .data$n,
      Percent = round2(.data$Percent),
      denominator = .data$denominator,
      n_missing = .data$n_missing
    )

  out
}

make_demographic_rows <- function(data, var_name, group_label) {
  x <- data[[var_name]]
  x_label <- format_category(x)
  missing_flag <- is.na(x_label) | x_label == ""
  nonmissing_n <- sum(!missing_flag)
  missing_n <- sum(missing_flag)

  if (is.factor(x)) {
    category_order <- format_category(levels(x))
    category_order <- category_order[!is.na(category_order) & category_order != ""]
  } else {
    category_order <- unique(x_label[!missing_flag])
  }

  header <- data.frame(
    Characteristic = group_label,
    Mean = "",
    SD = "",
    N = "",
    Percent = "",
    stringsAsFactors = FALSE
  )

  rows <- data.frame(.label = x_label, stringsAsFactors = FALSE) %>%
    filter(!is.na(.data$.label), .data$.label != "") %>%
    mutate(.label = factor(.data$.label, levels = unique(category_order))) %>%
    count(.data$.label, name = "n", .drop = TRUE) %>%
    mutate(
      Characteristic = paste0("  ", as.character(.data$.label)),
      Mean = "",
      SD = "",
      Percent = format_num(100 * .data$n / nonmissing_n, digits = 2),
      N = as.character(.data$n)
    ) %>%
    select(Characteristic, Mean, SD, N, Percent)

  if (missing_n > 0) {
    missing_row <- data.frame(
      Characteristic = "  Missing",
      Mean = "",
      SD = "",
      N = as.character(missing_n),
      Percent = format_num(100 * missing_n / nrow(data), digits = 2),
      stringsAsFactors = FALSE
    )
    rows <- bind_rows(rows, missing_row)
  }

  bind_rows(header, rows)
}

# -----------------------------------------------------------------------------
# Continuous descriptives
# -----------------------------------------------------------------------------

continuous_descriptives <- data.frame(
  variable = analysis_vars,
  n = vapply(EDData[analysis_vars], function(x) sum(!is.na(x)), integer(1)),
  mean = vapply(EDData[analysis_vars], mean, numeric(1), na.rm = TRUE),
  sd = vapply(EDData[analysis_vars], sd, numeric(1), na.rm = TRUE),
  min = vapply(EDData[analysis_vars], min, numeric(1), na.rm = TRUE),
  max = vapply(EDData[analysis_vars], max, numeric(1), na.rm = TRUE),
  skewness = vapply(EDData[analysis_vars], moments::skewness, numeric(1), na.rm = TRUE),
  kurtosis = vapply(EDData[analysis_vars], moments::kurtosis, numeric(1), na.rm = TRUE),
  row.names = NULL
) %>%
  mutate(
    across(c(mean, sd, min, max, skewness, kurtosis), round2)
  )

write.csv(
  continuous_descriptives,
  file.path(output_dir, "continuous_descriptives.csv"),
  row.names = FALSE
)

cat("Continuous descriptives:\n")
print(continuous_descriptives)
cat("\n")

# -----------------------------------------------------------------------------
# Clinical cutoff for eating disorder symptoms
# -----------------------------------------------------------------------------

# Scores greater than or equal to 15 are classified as meeting the cutoff.
clinical_edeqs <- sum(EDData$edeqs_total >= 15, na.rm = TRUE)
clinical_denom <- sum(!is.na(EDData$edeqs_total))
clinical_perc <- 100 * clinical_edeqs / clinical_denom

clinical_cutoff <- data.frame(
  cutoff_variable = "edeqs_total",
  cutoff_rule = ">= 15",
  n_at_or_above_cutoff = clinical_edeqs,
  denominator = clinical_denom,
  percent_at_or_above_cutoff = round2(clinical_perc),
  row.names = NULL
)

write.csv(
  clinical_cutoff,
  file.path(output_dir, "clinical_cutoff.csv"),
  row.names = FALSE
)

cat("Clinical cutoff summary:\n")
print(clinical_cutoff)
cat("\n")

# -----------------------------------------------------------------------------
# Manuscript Table: Demographic characteristics
# -----------------------------------------------------------------------------

age_mean <- mean(EDData$age, na.rm = TRUE)
age_sd <- sd(EDData$age, na.rm = TRUE)

Table_demographics <- bind_rows(
  data.frame(
    Characteristic = "Age (years)",
    Mean = format_num(age_mean, digits = 2),
    SD = format_num(age_sd, digits = 2),
    N = "",
    Percent = "",
    stringsAsFactors = FALSE
  ),
  data.frame(
    Characteristic = "",
    Mean = "",
    SD = "",
    N = "",
    Percent = "",
    stringsAsFactors = FALSE
  ),
  make_demographic_rows(EDData, "sex_factor", "Sex"),
  make_demographic_rows(EDData, "race_ethnicity", "Race/Ethnicity"),
  make_demographic_rows(EDData, "household_income_factor", "Household Annual Income"),
  make_demographic_rows(EDData, "education_factor", "Education")
)

write.csv(
  Table_demographics,
  file.path(output_dir, "Table_demographics.csv"),
  row.names = FALSE
)

cat("Demographics Table output:\n")
print(Table_demographics)
cat("\n")

# -----------------------------------------------------------------------------
# Correlation matrix and manuscript Table
# -----------------------------------------------------------------------------

corr_vars <- c("edeqs_total", "fbsv_total", "ptq_total", "rses_sc_score")
corr_labels <- c(
  "Eating Disorder Symptoms",
  "Bullying",
  "Rumination",
  "Self-esteem"
)

EDData_subset <- EDData[, corr_vars]
EDData_corr <- Hmisc::rcorr(as.matrix(EDData_subset), type = "pearson")

correlation_matrix <- as.data.frame(EDData_corr$r)
correlation_p_values <- as.data.frame(EDData_corr$P)
correlation_n <- as.data.frame(EDData_corr$n)

write.csv(
  correlation_matrix,
  file.path(output_dir, "correlation_matrix_r.csv"),
  row.names = TRUE
)

write.csv(
  correlation_p_values,
  file.path(output_dir, "correlation_matrix_p.csv"),
  row.names = TRUE
)

write.csv(
  correlation_n,
  file.path(output_dir, "correlation_matrix_n.csv"),
  row.names = TRUE
)

# Optional individual cor.test objects, saved as a compact CSV.
cor_test_pairs <- expand.grid(
  var1 = corr_vars,
  var2 = corr_vars,
  stringsAsFactors = FALSE
) %>%
  filter(match(var1, corr_vars) > match(var2, corr_vars))

cor_test_results <- cor_test_pairs %>%
  rowwise() %>%
  mutate(
    r = cor.test(EDData[[var1]], EDData[[var2]], method = "pearson")$estimate[[1]],
    p = cor.test(EDData[[var1]], EDData[[var2]], method = "pearson")$p.value
  ) %>%
  ungroup() %>%
  mutate(
    r = round2(r),
    p = signif(p, 3)
  )

write.csv(
  cor_test_results,
  file.path(output_dir, "cor_test_pairwise_results.csv"),
  row.names = FALSE
)

Table_correlations <- data.frame(
  Variable = paste0(seq_along(corr_labels), ". ", corr_labels),
  M = format_num(vapply(EDData[corr_vars], mean, numeric(1), na.rm = TRUE), digits = 2),
  SD = format_num(vapply(EDData[corr_vars], sd, numeric(1), na.rm = TRUE), digits = 2),
  `1` = "-",
  `2` = "-",
  `3` = "-",
  `4` = "-",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

for (i in seq_along(corr_vars)) {
  for (j in seq_along(corr_vars)) {
    if (j < i) {
      Table_correlations[i, as.character(j)] <- format_corr(
        EDData_corr$r[i, j],
        EDData_corr$P[i, j]
      )
    }
  }
}

write.csv(
  Table_correlations,
  file.path(output_dir, "Table_correlations.csv"),
  row.names = FALSE
)

cat("Correlation Table output:\n")
print(Table_correlations)
cat("\n")

# -----------------------------------------------------------------------------
# Diagnostic plots
# -----------------------------------------------------------------------------

plot_dir <- file.path(output_dir, "Diagnostic Plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

plot_vars <- c("edeqs_total", "fbsv_total", "ptq_total", "rses_sc_score", "age")
plot_labels <- c(
  edeqs_total = "Eating Disorder Symptoms",
  fbsv_total = "Bullying",
  ptq_total = "Rumination",
  rses_sc_score = "Self-esteem",
  age = "Age"
)

for (v in plot_vars) {
  x <- EDData[[v]]
  x_nonmissing <- x[!is.na(x)]

  png(file.path(plot_dir, paste0(v, "_histogram.png")), width = 900, height = 700)
  hist(
    x_nonmissing,
    main = paste("Histogram of", plot_labels[[v]]),
    xlab = plot_labels[[v]],
    breaks = max(10, round(sqrt(length(x_nonmissing))))
  )
  dev.off()

  png(file.path(plot_dir, paste0(v, "_boxplot.png")), width = 900, height = 700)
  boxplot(
    x_nonmissing,
    main = paste("Boxplot of", plot_labels[[v]]),
    ylab = plot_labels[[v]]
  )
  dev.off()

  if (length(unique(x_nonmissing)) > 1) {
    png(file.path(plot_dir, paste0(v, "_density.png")), width = 900, height = 700)
    plot(
      density(x_nonmissing, na.rm = TRUE),
      main = paste("Density Plot of", plot_labels[[v]]),
      xlab = plot_labels[[v]]
    )
    dev.off()
  }
}

# -----------------------------------------------------------------------------
# Rosner's tests for outliers
# Note: Rosner's test assumes approximate normality. These outputs should be
# interpreted as diagnostics, especially for skewed variables such as bullying.
# -----------------------------------------------------------------------------

rosner_vars <- c("edeqs_total", "fbsv_total", "ptq_total", "rses_sc_score")
rosner_results <- list()

for (v in rosner_vars) {
  x <- EDData[[v]]
  x_nonmissing <- x[!is.na(x)]

  if (length(x_nonmissing) > 10 && length(unique(x_nonmissing)) > 1) {
    result <- EnvStats::rosnerTest(x_nonmissing, k = 5)
    rosner_results[[v]] <- result

    result_table <- as.data.frame(result$all.stats)
    write.csv(
      result_table,
      file.path(output_dir, paste0("rosner_", v, ".csv")),
      row.names = FALSE
    )

    cat("Rosner's test for", v, ":\n")
    print(result)
    cat("\n")
  } else {
    cat("Rosner's test skipped for", v, "because there are too few observations or no variance.\n\n")
  }
}

# -----------------------------------------------------------------------------
# Variance inflation factor (VIF)
# Matches the primary mediation outcome model with age as a covariate.
# -----------------------------------------------------------------------------

vif_model <- lm(
  edeqs_total ~ fbsv_total + ptq_total + rses_sc_score + age,
  data = EDData
)

vif_values <- car::vif(vif_model)
vif_table <- data.frame(
  predictor = names(vif_values),
  vif = round2(as.numeric(vif_values)),
  row.names = NULL
)

write.csv(
  vif_table,
  file.path(output_dir, "vif_primary_outcome_model.csv"),
  row.names = FALSE
)

cat("VIF for primary outcome model:\n")
print(vif_table)
cat("\n")

# -----------------------------------------------------------------------------
# Cronbach's alpha
# Alpha is computed using complete item-level data for each scale. Therefore,
# alpha N may be smaller than the analytic N if item-level data are missing.
# -----------------------------------------------------------------------------

get_alpha <- function(data, item_names, scale_name) {
  missing_items <- setdiff(item_names, names(data))
  if (length(missing_items) > 0) {
    warning(
      "Skipping alpha for ", scale_name,
      "; missing item(s): ", paste(missing_items, collapse = ", ")
    )
    return(data.frame(
      scale = scale_name,
      raw_alpha = NA_real_,
      standardized_alpha = NA_real_,
      n_obs = NA_integer_,
      n_items = length(item_names),
      stringsAsFactors = FALSE
    ))
  }

  data_items <- data[, item_names]
  complete_item_n <- sum(complete.cases(data_items))

  if (complete_item_n < 2) {
    warning(
      "Skipping alpha for ", scale_name,
      "; fewer than 2 complete item-level observations."
    )
    return(data.frame(
      scale = scale_name,
      raw_alpha = NA_real_,
      standardized_alpha = NA_real_,
      n_obs = complete_item_n,
      n_items = length(item_names),
      stringsAsFactors = FALSE
    ))
  }

  result <- psych::alpha(data_items, check.keys = FALSE, use = "complete.obs")

  data.frame(
    scale = scale_name,
    raw_alpha = round2(result$total$raw_alpha),
    standardized_alpha = round2(result$total$std.alpha),
    n_obs = complete_item_n,
    n_items = length(item_names),
    stringsAsFactors = FALSE
  )
}

edeq_items <- paste0("edeqs_", 1:12)
fbsv_items <- paste0("fbsv_", 1:10)
rses_items <- paste0("rses_sc_", 1:10)
ptq_items <- paste0("ptq_", 1:15)

alpha_table <- bind_rows(
  get_alpha(EDData, edeq_items, "Eating Disorder Symptoms"),
  get_alpha(EDData, fbsv_items, "Bullying"),
  get_alpha(EDData, ptq_items, "Rumination"),
  get_alpha(EDData, rses_items, "Self-esteem")
)

write.csv(
  alpha_table,
  file.path(output_dir, "reliability_alpha.csv"),
  row.names = FALSE
)

cat("Cronbach's alpha values:\n")
print(alpha_table)
cat("\n")

# -----------------------------------------------------------------------------
# Bullying skew diagnostics
# -----------------------------------------------------------------------------

bullying_skew_diagnostics <- data.frame(
  statistic = c(
    "N nonmissing",
    "Minimum",
    "Maximum",
    "Mean",
    "SD",
    "Median",
    "Percent equal to zero",
    "Percent equal to minimum"
  ),
  value = c(
    sum(!is.na(EDData$fbsv_total)),
    min(EDData$fbsv_total, na.rm = TRUE),
    max(EDData$fbsv_total, na.rm = TRUE),
    mean(EDData$fbsv_total, na.rm = TRUE),
    sd(EDData$fbsv_total, na.rm = TRUE),
    median(EDData$fbsv_total, na.rm = TRUE),
    100 * mean(EDData$fbsv_total == 0, na.rm = TRUE),
    100 * mean(EDData$fbsv_total == min(EDData$fbsv_total, na.rm = TRUE), na.rm = TRUE)
  )
) %>%
  mutate(value = round2(value))

write.csv(
  bullying_skew_diagnostics,
  file.path(output_dir, "bullying_skew_diagnostics.csv"),
  row.names = FALSE
)

cat("Bullying skew diagnostics:\n")
print(bullying_skew_diagnostics)
cat("\n")

cat("All descriptives outputs saved to:", output_dir, "\n")
