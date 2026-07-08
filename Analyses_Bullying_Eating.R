# =========================================================
# Bullying mediation analyses
# =========================================================

# ----------------------------
# Load libraries
# ----------------------------
library(mediation)
library(dplyr)
library(lavaan)
library(mice)

# ----------------------------
# Read in data
# ----------------------------
EDData <- readRDS("~/Desktop/ED Paper/ED_fbsv_rses_ptq.rds")

# ----------------------------
# Settings
# ----------------------------
set.seed(1234)

# Number of bootstrap resamples
B <- 5000

vars_primary <- c(
  "fbsv_total",
  "rses_sc_score",
  "ptq_total",
  "edeqs_total",
  "age"
)

# Keep complete cases for the variables used in these models.
# If EDData is already complete from the data-prep script, this should not change N.
EDData_model <- EDData %>%
  dplyr::select(dplyr::all_of(vars_primary)) %>%
  tidyr::drop_na()

# Standardized data for primary and sensitivity analyses.
EDData_std <- as.data.frame(scale(EDData_model))

# =========================================================
# Helper functions
# =========================================================
# The helper functions below are retained for the single-mediator
# sensitivity analyses, which continue to use the mediation package.

sig_code <- function(p_value) {
  if (length(p_value) == 0 || is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01)  return("**")
  if (p_value < 0.05)  return("*")
  if (p_value < 0.1)   return(".")
  return("")
}

format_pval <- function(p) {
  if (length(p) == 0 || is.na(p)) return(NA)
  if (p < 0.001) return("<.001")
  return(round(p, 3))
}

# Extract one coefficient from an lm model.
extract_path <- function(model, var_name) {
  coef_table <- summary(model)$coefficients

  if (!var_name %in% rownames(coef_table)) {
    stop(paste("Variable", var_name, "not found in model coefficients."))
  }

  coef_info <- coef_table[var_name, ]

  estimate <- coef_info["Estimate"]
  se <- coef_info["Std. Error"]
  t_val <- coef_info["t value"]
  p_val <- coef_info["Pr(>|t|)"]
  ci <- estimate + c(-1.96, 1.96) * se

  data.frame(
    B = round(estimate, 3),
    SE = round(se, 3),
    Lower = round(ci[1], 3),
    Upper = round(ci[2], 3),
    t = round(t_val, 3),
    p = format_pval(p_val),
    Signif = sig_code(p_val),
    stringsAsFactors = FALSE
  )
}

# Table for mediation output from the mediation package.
make_mediation_table <- function(model_m, model_y, med_out, treat, mediator) {
  path_a <- extract_path(model_m, treat)
  path_b <- extract_path(model_y, mediator)
  path_direct <- extract_path(model_y, treat)

  s <- summary(med_out)

  est_ind <- s$d0
  se_ind <- sd(s$d0.sims, na.rm = TRUE)
  t_ind <- est_ind / se_ind
  p_ind <- s$d0.p
  ci_ind <- s$d0.ci

  est_ade <- s$z0
  se_ade <- sd(s$z0.sims, na.rm = TRUE)
  t_ade <- est_ade / se_ade
  p_ade <- s$z0.p
  ci_ade <- s$z0.ci

  est_tot <- s$tau.coef
  se_tot <- sd(s$tau.sims, na.rm = TRUE)
  t_tot <- est_tot / se_tot
  p_tot <- s$tau.p
  ci_tot <- quantile(s$tau.sims, c(.025, .975), na.rm = TRUE)

  data.frame(
    Effect = c(
      "Indirect (ACME)",
      "Component a",
      "Component b",
      "Direct effect (ADE)",
      "Total effect"
    ),
    B = round(c(
      est_ind,
      path_a$B,
      path_b$B,
      est_ade,
      est_tot
    ), 4),
    SE = round(c(
      se_ind,
      path_a$SE,
      path_b$SE,
      se_ade,
      se_tot
    ), 3),
    Lower = round(c(
      ci_ind[1],
      path_a$Lower,
      path_b$Lower,
      ci_ade[1],
      ci_tot[1]
    ), 3),
    Upper = round(c(
      ci_ind[2],
      path_a$Upper,
      path_b$Upper,
      ci_ade[2],
      ci_tot[2]
    ), 3),
    t = round(c(
      t_ind,
      path_a$t,
      path_b$t,
      t_ade,
      t_tot
    ), 3),
    p = c(
      format_pval(p_ind),
      path_a$p,
      path_b$p,
      format_pval(p_ade),
      format_pval(p_tot)
    ),
    Signif = c(
      sig_code(p_ind),
      path_a$Signif,
      path_b$Signif,
      sig_code(p_ade),
      sig_code(p_tot)
    ),
    stringsAsFactors = FALSE
  )
}

# Explicitly extract the proportion mediated from a mediation object.
# Percent_Mediated is the value to report in manuscript text.
make_proportion_mediated_table <- function(med_out, effect_label) {
  s <- summary(med_out)

  acme <- as.numeric(s$d0)
  total_effect <- as.numeric(s$tau.coef)
  prop_mediated <- as.numeric(s$n0)
  prop_ci <- as.numeric(s$n0.ci)
  prop_p <- as.numeric(s$n0.p)

  data.frame(
    Effect = effect_label,
    ACME = round(acme, 4),
    Total_Effect = round(total_effect, 4),
    Prop_Mediated = round(prop_mediated, 4),
    Percent_Mediated = round(100 * prop_mediated, 2),
    Percent_Mediated_from_ACME_over_Total = round(100 * acme / total_effect, 2),
    Lower_Percent = round(100 * prop_ci[1], 2),
    Upper_Percent = round(100 * prop_ci[2], 2),
    p = format_pval(prop_p),
    Signif = sig_code(prop_p),
    stringsAsFactors = FALSE
  )
}

make_cooks_summary <- function(model, model_name) {
  cooks_d <- cooks.distance(model)
  cutoff_4n <- 4 / length(cooks_d)

  data.frame(
    Model = model_name,
    N = length(cooks_d),
    Cutoff_4_over_N = cutoff_4n,
    N_above_4_over_N = sum(cooks_d > cutoff_4n, na.rm = TRUE),
    N_above_0_5 = sum(cooks_d > 0.5, na.rm = TRUE),
    N_above_1 = sum(cooks_d > 1, na.rm = TRUE),
    Max_Cooks_D = max(cooks_d, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

make_cooks_manuscript_summary <- function(cooks_summary, analysis_label) {
  max_cooks <- max(cooks_summary$Max_Cooks_D, na.rm = TRUE)
  n_above_0_5 <- sum(cooks_summary$N_above_0_5, na.rm = TRUE)
  n_above_1 <- sum(cooks_summary$N_above_1, na.rm = TRUE)

  data.frame(
    Analysis = analysis_label,
    Max_Cooks_D = round(max_cooks, 3),
    N_above_0_5 = n_above_0_5,
    N_above_1 = n_above_1,
    Manuscript_Sentence = sprintf(
      "Influence diagnostics did not indicate that the results were driven by highly influential observations: the Cook's distance values were equal to or lower than %.3f across the mediator and outcome models, with no observations exceeding conventional thresholds of 0.50 or 1.00.",
      max_cooks
    ),
    stringsAsFactors = FALSE
  )
}

# =========================================================
# Age as a predictor
# =========================================================

mod_age_edeqs <- lm(edeqs_total ~ age, data = EDData_std)
mod_age_fbsv <- lm(fbsv_total ~ age, data = EDData_std)
mod_age_ptq <- lm(ptq_total ~ age, data = EDData_std)
mod_age_rses <- lm(rses_sc_score ~ age, data = EDData_std)

summary(mod_age_edeqs)
summary(mod_age_fbsv)
summary(mod_age_ptq)
summary(mod_age_rses)

# =========================================================
# Primary parallel mediation, standardized scale
# =========================================================

parallel_model <- '

  # Mediator models
  rses_sc_score ~ a1*fbsv_total + age
  ptq_total ~ a2*fbsv_total + age

  # Outcome model
  edeqs_total ~ cprime*fbsv_total +
                b1*rses_sc_score +
                b2*ptq_total +
                age

# Specific indirect effects
ind_rses := a1*b1
ind_ptq  := a2*b2

# Total indirect effect
ind_total := ind_rses + ind_ptq

# Total effect
total := cprime + ind_total

# Contrast between indirect effects
contrast := ind_rses - ind_ptq

# Proportion mediated
prop_rses := ind_rses / total
prop_ptq := ind_ptq / total
prop_total := ind_total / total

'

fit_parallel <- sem(
  parallel_model,
  data = EDData_std,
  se = "bootstrap",
  bootstrap = B
)

summary(
  fit_parallel,
  standardized = TRUE,
  ci = TRUE
)

parallel_results <- parameterEstimates(
  fit_parallel,
  standardized = TRUE,
  ci = TRUE
)

parallel_results

# Extract estimates
parallel_prop_mediated_table <- subset(
  parallel_results,
  label %in% c(
    "prop_rses",
    "prop_ptq",
    "prop_total"
  )
) %>%
  mutate(
    Effect = c(
      "Primary parallel mediation: self-esteem",
      "Primary parallel mediation: repetitive negative thinking",
      "Primary parallel mediation: total indirect effect"
    ),
    Prop_Mediated = round(est, 4),
    Percent_Mediated = round(est * 100, 3),
    Lower_Percent = round(ci.lower * 100, 3),
    Upper_Percent = round(ci.upper * 100, 3),
    p = ifelse(
      pvalue < .001,
      "<.001",
      sprintf("%.3f", pvalue)
    ),
    Signif = case_when(
      pvalue < .001 ~ "***",
      pvalue < .01 ~ "**",
      pvalue < .05 ~ "*",
      pvalue < .10 ~ ".",
      TRUE ~ ""
    )
  )

parallel_prop_mediated_table

# Primary effects

primary_effects <- subset(
  parallel_results,
  label %in% c(
    "cprime",
    "ind_rses",
    "ind_ptq",
    "ind_total",
    "total",
    "contrast",
    "prop_rses",
    "prop_ptq",
    "prop_total"
  )
)

primary_effects

# =========================================================
# Single-mediator sensitivity analyses, standardized scale
# =========================================================

# Single-mediator RSES model
single_mediator_rses_mediator_model <- lm(
  rses_sc_score ~ fbsv_total + age,
  data = EDData_std
)

single_mediator_rses_outcome_model <- lm(
  edeqs_total ~ fbsv_total + rses_sc_score + age,
  data = EDData_std
)

single_mediator_rses_mediation <- mediate(
  single_mediator_rses_mediator_model,
  single_mediator_rses_outcome_model,
  treat = "fbsv_total",
  mediator = "rses_sc_score",
  boot = TRUE,
  sims = B
)

single_mediator_rses_table <- make_mediation_table(
  model_m = single_mediator_rses_mediator_model,
  model_y = single_mediator_rses_outcome_model,
  med_out = single_mediator_rses_mediation,
  treat = "fbsv_total",
  mediator = "rses_sc_score"
)

# Single-mediator PTQ model
single_mediator_ptq_mediator_model <- lm(
  ptq_total ~ fbsv_total + age,
  data = EDData_std
)

single_mediator_ptq_outcome_model <- lm(
  edeqs_total ~ fbsv_total + ptq_total + age,
  data = EDData_std
)

single_mediator_ptq_mediation <- mediate(
  single_mediator_ptq_mediator_model,
  single_mediator_ptq_outcome_model,
  treat = "fbsv_total",
  mediator = "ptq_total",
  boot = TRUE,
  sims = B
)

single_mediator_ptq_table <- make_mediation_table(
  model_m = single_mediator_ptq_mediator_model,
  model_y = single_mediator_ptq_outcome_model,
  med_out = single_mediator_ptq_mediation,
  treat = "fbsv_total",
  mediator = "ptq_total"
)

single_mediator_rses_table
single_mediator_ptq_table

single_mediator_prop_mediated_table <- bind_rows(
  make_proportion_mediated_table(
    single_mediator_rses_mediation,
    "Single-mediator sensitivity analysis: self-esteem"
  ),
  make_proportion_mediated_table(
    single_mediator_ptq_mediation,
    "Single-mediator sensitivity analysis: repetitive negative thinking"
  )
)

single_mediator_prop_mediated_table

# =========================================================
# FBSV log-transformed parallel mediation, standardized scale
# =========================================================
EDData_log_model <- EDData_model %>%
  mutate(fbsv_log = log(fbsv_total + 1)) %>%
  dplyr::select(
    fbsv_log,
    rses_sc_score,
    ptq_total,
    edeqs_total,
    age
  )

EDData_log_std <- as.data.frame(scale(EDData_log_model))

log_parallel_model <- '

  # Mediator models
  rses_sc_score ~ a1*fbsv_log + age
  ptq_total ~ a2*fbsv_log + age

  # Outcome model
  edeqs_total ~ cprime*fbsv_log +
                b1*rses_sc_score +
                b2*ptq_total +
                age

  # Specific indirect effects
  ind_rses := a1*b1
  ind_ptq  := a2*b2

  # Total indirect effect
  ind_total := ind_rses + ind_ptq

  # Total effect
  total := cprime + ind_total

  # Contrast between indirect effects
  contrast := ind_rses - ind_ptq

'

fit_log_parallel <- sem(
  log_parallel_model,
  data = EDData_log_std,
  se = "bootstrap",
  bootstrap = B
)

summary(
  fit_log_parallel,
  standardized = TRUE,
  ci = TRUE
)

log_parallel_results <- parameterEstimates(
  fit_log_parallel,
  standardized = TRUE,
  ci = TRUE
)

log_parallel_results

# ========================================================
# Multiple imputation sensitivity analysis
# =========================================================

# Read age-eligible dataset with missingness retained

EDData_MI <- readRDS("~/Desktop/ED Paper/ED_fbsv_rses_ptq_MI.rds")

# Variables used in the primary mediation model

mi_vars <- c(
  "fbsv_total",
  "rses_sc_score",
  "ptq_total",
  "edeqs_total",
  "age"
)

# Standardize variables so MI estimates are comparable to the
# standardized primary mediation model

EDData_MI_std <- EDData_MI

EDData_MI_std[, mi_vars] <- scale(
  EDData_MI_std[, mi_vars]
)

# Examine missing-data pattern

md.pattern(EDData_MI_std[, mi_vars])

# Multiple imputation using predictive mean matching

mi_imputations <- mice(
  EDData_MI_std[, mi_vars],
  m = 20,
  method = "pmm",
  seed = 1234
)

mi_imputations

# =========================================================
# Fit parallel mediation model in each imputed dataset
# =========================================================

mi_effects <- list()

for(i in 1:mi_imputations$m){
  
  dat_i <- complete(mi_imputations, i)
  
  fit_i <- sem(
    parallel_model,
    data = dat_i
  )
  
  pe <- parameterEstimates(fit_i)
  
  effects_i <- subset(
    pe,
    label %in% c(
      "a1",
      "a2",
      "b1",
      "b2",
      "cprime",
      "ind_rses",
      "ind_ptq",
      "ind_total",
      "total"
    )
  )[, c("label", "est", "se")]
  
  effects_i$imputation <- i
  
  mi_effects[[i]] <- effects_i
}

mi_effects_df <- bind_rows(mi_effects)

# =========================================================
# Rubin pooling function
# =========================================================

pool_rubin <- function(est, se){
  
  m <- length(est)
  
  qbar <- mean(est)
  
  Ubar <- mean(se^2)
  
  B <- var(est)
  
  Tvar <- Ubar + (1 + 1/m) * B
  
  pooled_se <- sqrt(Tvar)
  
  z <- qbar / pooled_se
  
  p <- 2 * (1 - pnorm(abs(z)))
  
  lower <- qbar - 1.96 * pooled_se
  upper <- qbar + 1.96 * pooled_se
  
  data.frame(
    Estimate = qbar,
    SE = pooled_se,
    Lower = lower,
    Upper = upper,
    z = z,
    p = p
  )
}

# =========================================================
# Pool estimates across imputations
# =========================================================

mi_pooled_effects <- mi_effects_df %>%
  group_by(label) %>%
  summarise(
    pooled = list(
      pool_rubin(est, se)
    ),
    .groups = "drop"
  ) %>%
  tidyr::unnest(pooled)

# =========================================================
# Create table
# =========================================================

mi_pooled_effects <- mi_pooled_effects %>%
  mutate(
    Parameter = recode(
      label,
      a1 = "a1: FBSV → Self-esteem",
      a2 = "a2: FBSV → Repetitive negative thinking",
      b1 = "b1: Self-esteem → Eating-disorder symptoms",
      b2 = "b2: Repetitive negative thinking → Eating-disorder symptoms",
      cprime = "Direct effect (c′)",
      ind_rses = "Indirect via self-esteem",
      ind_ptq = "Indirect via repetitive negative thinking",
      ind_total = "Total indirect effect",
      total = "Total effect"
    ),
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    Lower = round(Lower, 3),
    Upper = round(Upper, 3),
    z = round(z, 3),
    p = ifelse(p < .001, "<.001", round(p, 3))
  )

mi_pooled_effects

# =========================================================
# Influence diagnostics: Cook's distance
# =========================================================
# Regression models
primary_mediator_rses_model <- lm(
  rses_sc_score ~ fbsv_total + age,
  data = EDData_std
)

primary_mediator_ptq_model <- lm(
  ptq_total ~ fbsv_total + age,
  data = EDData_std
)

primary_outcome_model <- lm(
  edeqs_total ~ fbsv_total + rses_sc_score + ptq_total + age,
  data = EDData_std
)

# Regression models corresponding to the log-transformed sensitivity analysis

log_parallel_mediator_rses_model <- lm(
  rses_sc_score ~ fbsv_log + age,
  data = EDData_log_std
)

log_parallel_mediator_ptq_model <- lm(
  ptq_total ~ fbsv_log + age,
  data = EDData_log_std
)

log_parallel_outcome_model <- lm(
  edeqs_total ~ fbsv_log + rses_sc_score + ptq_total + age,
  data = EDData_log_std
)

# Concise Cook's distance summaries. These retain the maximum Cook's distance
# and counts above conventional thresholds without printing top-observation lists
# or diagnostic plots.
primary_cooks_summary <- bind_rows(
  make_cooks_summary(primary_mediator_rses_model, "Primary mediator model: self-esteem"),
  make_cooks_summary(primary_mediator_ptq_model, "Primary mediator model: repetitive negative thinking"),
  make_cooks_summary(primary_outcome_model, "Primary outcome model: parallel mediators")
)

single_mediator_cooks_summary <- bind_rows(
  make_cooks_summary(single_mediator_rses_mediator_model, "Single-mediator self-esteem: mediator model"),
  make_cooks_summary(single_mediator_rses_outcome_model, "Single-mediator self-esteem: outcome model"),
  make_cooks_summary(single_mediator_ptq_mediator_model, "Single-mediator repetitive negative thinking: mediator model"),
  make_cooks_summary(single_mediator_ptq_outcome_model, "Single-mediator repetitive negative thinking: outcome model")
)

log_parallel_cooks_summary <- bind_rows(
  make_cooks_summary(log_parallel_mediator_rses_model, "Log-FBSV mediator model: self-esteem"),
  make_cooks_summary(log_parallel_mediator_ptq_model, "Log-FBSV mediator model: repetitive negative thinking"),
  make_cooks_summary(log_parallel_outcome_model, "Log-FBSV outcome model: parallel mediators")
)

cooks_summary_all <- bind_rows(
  primary_cooks_summary,
  single_mediator_cooks_summary,
  log_parallel_cooks_summary
) %>%
  mutate(
    Cutoff_4_over_N = round(Cutoff_4_over_N, 4),
    Max_Cooks_D = round(Max_Cooks_D, 3)
  )

primary_cooks_manuscript_summary <- make_cooks_manuscript_summary(
  primary_cooks_summary,
  "Primary parallel mediation"
)

single_mediator_cooks_manuscript_summary <- make_cooks_manuscript_summary(
  single_mediator_cooks_summary,
  "Single-mediator sensitivity analyses"
)

log_parallel_cooks_manuscript_summary <- make_cooks_manuscript_summary(
  log_parallel_cooks_summary,
  "Log-FBSV sensitivity analysis"
)

cooks_manuscript_summaries <- bind_rows(
  primary_cooks_manuscript_summary,
  single_mediator_cooks_manuscript_summary,
  log_parallel_cooks_manuscript_summary
)

cooks_summary_all
cooks_manuscript_summaries

cat(
  "\nCook's distance sentence for manuscript text, primary model:\n",
  primary_cooks_manuscript_summary$Manuscript_Sentence,
  "\n",
  sep = ""
)
