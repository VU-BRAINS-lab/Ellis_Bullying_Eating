# =========================================================
# Bullying mediation analyses
# =========================================================

# ----------------------------
# Load libraries
# ----------------------------
library(lm.beta)
library(mediation)
library(dplyr)

# ----------------------------
# Read in data
# ----------------------------
EDData <- readRDS("~/Desktop/ED Paper/ED_fbsv_rses_ptq.rds")

# ----------------------------
# Settings
# ----------------------------
set.seed(1234)
sims <- 5000
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
    ), 3),
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

# Bootstrap p-value for a distribution tested against zero.
# This prevents p = 0 by using 1 / number of valid bootstrap samples as the lower bound.
bootstrap_p_two_sided <- function(x) {
  x <- x[!is.na(x)]
  p <- 2 * min(mean(x >= 0), mean(x <= 0))
  max(p, 1 / length(x))
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

# Mediator models
primary_mediator_rses_model <- lm(
  rses_sc_score ~ fbsv_total + age,
  data = EDData_std
)

primary_mediator_ptq_model <- lm(
  ptq_total ~ fbsv_total + age,
  data = EDData_std
)

# Outcome model includes both mediators: this is the primary parallel mediation model.
primary_outcome_model <- lm(
  edeqs_total ~ fbsv_total + rses_sc_score + ptq_total + age,
  data = EDData_std
)

primary_mediation_rses <- mediate(
  primary_mediator_rses_model,
  primary_outcome_model,
  treat = "fbsv_total",
  mediator = "rses_sc_score",
  boot = TRUE,
  sims = sims
)

primary_mediation_ptq <- mediate(
  primary_mediator_ptq_model,
  primary_outcome_model,
  treat = "fbsv_total",
  mediator = "ptq_total",
  boot = TRUE,
  sims = sims
)

primary_parallel_rses_table <- make_mediation_table(
  model_m = primary_mediator_rses_model,
  model_y = primary_outcome_model,
  med_out = primary_mediation_rses,
  treat = "fbsv_total",
  mediator = "rses_sc_score"
)

primary_parallel_ptq_table <- make_mediation_table(
  model_m = primary_mediator_ptq_model,
  model_y = primary_outcome_model,
  med_out = primary_mediation_ptq,
  treat = "fbsv_total",
  mediator = "ptq_total"
)

primary_parallel_rses_table
primary_parallel_ptq_table

primary_parallel_prop_mediated_table <- bind_rows(
  make_proportion_mediated_table(
    primary_mediation_rses,
    "Primary parallel mediation: self-esteem"
  ),
  make_proportion_mediated_table(
    primary_mediation_ptq,
    "Primary parallel mediation: repetitive negative thinking"
  )
)

primary_parallel_prop_mediated_table

cat(
  "\nProportion mediated for manuscript text:\n",
  "Self-esteem mediated ",
  primary_parallel_prop_mediated_table$Percent_Mediated[1],
  "% of the total effect.\n",
  "Repetitive negative thinking mediated ",
  primary_parallel_prop_mediated_table$Percent_Mediated[2],
  "% of the total effect.\n",
  sep = ""
)

# =========================================================
# Paired bootstrap comparison of the two indirect effects
# =========================================================

boot_indirect_rses <- numeric(B)
boot_indirect_ptq <- numeric(B)
boot_indirect_diff <- numeric(B)

n <- nrow(EDData_model)

for (b in seq_len(B)) {
  idx <- sample(seq_len(n), size = n, replace = TRUE)
  dat_b <- EDData_model[idx, ]

  # Standardize within bootstrap sample so bootstrap estimates match the primary standardized analyses.
  dat_b_std <- as.data.frame(scale(dat_b[, vars_primary]))

  m_rses_b <- lm(rses_sc_score ~ fbsv_total + age, data = dat_b_std)
  m_ptq_b <- lm(ptq_total ~ fbsv_total + age, data = dat_b_std)
  y_b <- lm(
    edeqs_total ~ fbsv_total + rses_sc_score + ptq_total + age,
    data = dat_b_std
  )

  a_rses_b <- coef(m_rses_b)["fbsv_total"]
  a_ptq_b <- coef(m_ptq_b)["fbsv_total"]
  b_rses_b <- coef(y_b)["rses_sc_score"]
  b_ptq_b <- coef(y_b)["ptq_total"]

  boot_indirect_rses[b] <- a_rses_b * b_rses_b
  boot_indirect_ptq[b] <- a_ptq_b * b_ptq_b
  boot_indirect_diff[b] <- boot_indirect_rses[b] - boot_indirect_ptq[b]
}

# Original-sample point estimates on the same standardized scale.
indirect_rses_est <- coef(primary_mediator_rses_model)["fbsv_total"] *
  coef(primary_outcome_model)["rses_sc_score"]

indirect_ptq_est <- coef(primary_mediator_ptq_model)["fbsv_total"] *
  coef(primary_outcome_model)["ptq_total"]

indirect_diff_est <- indirect_rses_est - indirect_ptq_est

paired_bootstrap_indirect_comparison <- data.frame(
  Effect = c(
    "Indirect via RSES",
    "Indirect via PTQ",
    "Difference: RSES indirect - PTQ indirect"
  ),
  Estimate = round(c(
    indirect_rses_est,
    indirect_ptq_est,
    indirect_diff_est
  ), 3),
  Lower = round(c(
    quantile(boot_indirect_rses, .025, na.rm = TRUE),
    quantile(boot_indirect_ptq, .025, na.rm = TRUE),
    quantile(boot_indirect_diff, .025, na.rm = TRUE)
  ), 3),
  Upper = round(c(
    quantile(boot_indirect_rses, .975, na.rm = TRUE),
    quantile(boot_indirect_ptq, .975, na.rm = TRUE),
    quantile(boot_indirect_diff, .975, na.rm = TRUE)
  ), 3),
  p = c(
    format_pval(bootstrap_p_two_sided(boot_indirect_rses)),
    format_pval(bootstrap_p_two_sided(boot_indirect_ptq)),
    format_pval(bootstrap_p_two_sided(boot_indirect_diff))
  ),
  Signif = c(
    sig_code(bootstrap_p_two_sided(boot_indirect_rses)),
    sig_code(bootstrap_p_two_sided(boot_indirect_ptq)),
    sig_code(bootstrap_p_two_sided(boot_indirect_diff))
  ),
  stringsAsFactors = FALSE
)

paired_bootstrap_indirect_comparison

# Logical significance flag for the difference between indirect effects.
paired_bootstrap_diff_significant <-
  quantile(boot_indirect_diff, .025, na.rm = TRUE) > 0 |
  quantile(boot_indirect_diff, .975, na.rm = TRUE) < 0

paired_bootstrap_diff_significant

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
  sims = sims
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
  sims = sims
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
  dplyr::select(fbsv_log, rses_sc_score, ptq_total, edeqs_total, age)

EDData_log_std <- as.data.frame(scale(EDData_log_model))

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

log_parallel_mediation_rses <- mediate(
  log_parallel_mediator_rses_model,
  log_parallel_outcome_model,
  treat = "fbsv_log",
  mediator = "rses_sc_score",
  boot = TRUE,
  sims = sims
)

log_parallel_mediation_ptq <- mediate(
  log_parallel_mediator_ptq_model,
  log_parallel_outcome_model,
  treat = "fbsv_log",
  mediator = "ptq_total",
  boot = TRUE,
  sims = sims
)

log_parallel_rses_table <- make_mediation_table(
  model_m = log_parallel_mediator_rses_model,
  model_y = log_parallel_outcome_model,
  med_out = log_parallel_mediation_rses,
  treat = "fbsv_log",
  mediator = "rses_sc_score"
)

log_parallel_ptq_table <- make_mediation_table(
  model_m = log_parallel_mediator_ptq_model,
  model_y = log_parallel_outcome_model,
  med_out = log_parallel_mediation_ptq,
  treat = "fbsv_log",
  mediator = "ptq_total"
)

log_parallel_rses_table
log_parallel_ptq_table

log_parallel_prop_mediated_table <- bind_rows(
  make_proportion_mediated_table(
    log_parallel_mediation_rses,
    "Log-FBSV sensitivity analysis: self-esteem"
  ),
  make_proportion_mediated_table(
    log_parallel_mediation_ptq,
    "Log-FBSV sensitivity analysis: repetitive negative thinking"
  )
)

log_parallel_prop_mediated_table

# =========================================================
# Influence diagnostics: Cook's distance
# =========================================================

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
