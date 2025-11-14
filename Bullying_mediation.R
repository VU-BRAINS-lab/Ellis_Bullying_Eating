# Load libraries
library(mediation)

# Read in data
EDData <- readRDS("~/Desktop/ED Paper/ED_fbsv_rses.rds")

# Age as a predictor
mod_age_edeqs <- lm(edeqs_total ~ age, data = EDData)
mod_age_fbsv <- lm(fbsv_log ~ age, data = EDData)
mod_age_ptq <- lm(ptq_total ~ age, data = EDData)
mod_age_rses <- lm(rses_sc_score ~ age, data = EDData)

# Print model results
summary(mod_age_edeqs)
summary(mod_age_fbsv)
summary(mod_age_ptq)
summary(mod_age_rses)

## Mediations

# ----------------------------
# Helper functions
# ----------------------------

set.seed(1234)

# Significance code
signif <- function(p_value) {
  if (length(p_value) == 0 || is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01)  return("**")
  if (p_value < 0.05)  return("*")
  if (p_value < 0.1)   return(".")
  return("")
}

# Format p-value for table
format_pval <- function(p) {
  if (is.na(p)) return(NA)
  if (p < 0.001) return("<.001")
  return(round(p, 3))
}

# Standardized beta
std_beta <- function(coef, sd_x, sd_y) {
  coef * (sd_x / sd_y)
}

# Extract path info from lm model
extract_path <- function(model, var_name) {
  coef_info <- summary(model)$coefficients[var_name, ]
  estimate <- coef_info["Estimate"]
  se <- coef_info["Std. Error"]
  t_val <- estimate / se
  p_val <- coef_info["Pr(>|t|)"]
  ci <- estimate + c(-1.96, 1.96) * se
  beta <- std_beta(estimate, sd(model$model[[var_name]]), sd(model$model[[1]]))
  
  data.frame(
    Estimate = round(estimate, 3),
    SE = round(se, 3),
    Lower = round(ci[1], 3),
    Upper = round(ci[2], 3),
    B = round(beta, 3),
    t = round(t_val, 3),
    p = format_pval(p_val),
    Signif = signif(p_val)
  )
}

# ----------------------------
# Full table function
# ----------------------------
make_full_table <- function(model.m, model.y, med.out, med.out.std, treat, mediator) {
  
  # Path coefficients a and b
  path_a <- extract_path(model.m, treat)
  path_b <- extract_path(model.y, mediator)
  path_y <- extract_path(model.y, treat)  # c' path for direct effect
  
  s <- summary(med.out)
  s.std <- summary(med.out.std)
  
  # -------------------
  # Indirect (ACME)
  # -------------------
  est_ind <- if (!is.null(s$d0)) s$d0 else NA
  se_ind  <- if (!is.null(s$d0.sims)) sd(s$d0.sims) else NA
  t_ind   <- if (!is.na(est_ind) && !is.na(se_ind) && se_ind != 0) est_ind / se_ind else NA
  p_ind   <- if (!is.na(t_ind)) 2 * (1 - pnorm(abs(t_ind))) else NA
  ci_ind  <- if (!is.null(s$d0.ci)) s$d0.ci else c(NA, NA)
  beta_ind <- if (!is.null(s.std$d0)) s.std$d0 else NA
  
  # -------------------
  # ADE (direct effect)
  # -------------------
  est_ade <- if (!is.null(s$z0)) s$z0 else NA
  se_ade  <- if (!is.null(s$z0.sims)) sd(s$z0.sims) else NA
  t_ade   <- if (!is.na(est_ade) && !is.na(se_ade) && se_ade != 0) est_ade / se_ade else NA
  p_ade   <- if (!is.na(t_ade)) 2 * (1 - pnorm(abs(t_ade))) else NA
  ci_ade  <- if (!is.null(s$z0.ci)) s$z0.ci else c(NA, NA)
  beta_ade <- if (!is.null(s.std$z0)) s.std$z0 else NA
  
  # -------------------
  # Total effect
  # -------------------
  est_tot <- if (!is.null(s$tau.coef)) s$tau.coef else NA
  se_tot <- if (!is.null(s$tau.sims)) sd(s$tau.sims) else NA
  t_tot <- if (!is.na(est_tot) && !is.na(se_tot) && se_tot != 0) est_tot / se_tot else NA
  p_tot <- if (!is.na(t_tot)) 2 * (1 - pnorm(abs(t_tot))) else NA
  ci_tot <- if (!is.null(s$tau.sims)) quantile(s$tau.sims, c(0.025,0.975)) else c(NA, NA)
  beta_tot <- if (!is.null(s.std$tau.coef)) s.std$tau.coef else NA
  
  # -------------------
  # Combine into table
  # -------------------
  data.frame(
    Effect = c("Indirect (ACME)", "Component a", "Component b", "ADE", "Total"),
    Estimate = round(c(est_ind, path_a$Estimate, path_b$Estimate, est_ade, est_tot), 3),
    SE       = round(c(se_ind, path_a$SE, path_b$SE, se_ade, se_tot), 3),
    Lower    = round(c(ci_ind[1], path_a$Lower, path_b$Lower, ci_ade[1], ci_tot[1]), 3),
    Upper    = round(c(ci_ind[2], path_a$Upper, path_b$Upper, ci_ade[2], ci_tot[2]), 3),
    B        = round(c(beta_ind, path_a$B, path_b$B, beta_ade, beta_tot), 3),
    t        = round(c(t_ind, path_a$t, path_b$t, t_ade, t_tot), 3),
    p        = c(format_pval(p_ind), path_a$p, path_b$p, format_pval(p_ade), format_pval(p_tot)),
    Signif   = c(signif(p_ind), path_a$Signif, path_b$Signif, signif(p_ade), signif(p_tot))
  )
}

# ----------------------------
# Fit models
# ----------------------------

# Mediator models
model.m1 <- lm(rses_sc_score ~ fbsv_res, data = EDData)
model.m2 <- lm(ptq_total ~ fbsv_res, data = EDData)

# Outcome model including both mediators
model.y <- lm(edeqs_total ~ fbsv_res + rses_sc_score + ptq_total, data = EDData)

# Standardized versions for betas
EDData_std <- as.data.frame(scale(EDData[, c("fbsv_res", "rses_sc_score", "ptq_total", "edeqs_total")]))
model.m1.std <- lm(rses_sc_score ~ fbsv_res, data = EDData_std)
model.m2.std <- lm(ptq_total ~ fbsv_res, data = EDData_std)
model.y.std  <- lm(edeqs_total ~ fbsv_res + rses_sc_score + ptq_total, data = EDData_std)

# ----------------------------
# Mediation analysis
# ----------------------------
med.out1 <- mediate(model.m1, model.y, treat = "fbsv_res", mediator = "rses_sc_score", boot = TRUE, sims = 5000)
med.out2 <- mediate(model.m2, model.y, treat = "fbsv_res", mediator = "ptq_total", boot = TRUE, sims = 5000)

med.out1.std <- mediate(model.m1.std, model.y.std, treat = "fbsv_res", mediator = "rses_sc_score", boot = TRUE, sims = 5000)
med.out2.std <- mediate(model.m2.std, model.y.std, treat = "fbsv_res", mediator = "ptq_total", boot = TRUE, sims = 5000)

# ----------------------------
# Generate clean tables
# ----------------------------
tab1 <- make_full_table(model.m1, model.y, med.out1, med.out1.std, "fbsv_res", "rses_sc_score")
tab2 <- make_full_table(model.m2, model.y, med.out2, med.out2.std, "fbsv_res", "ptq_total")

tab1
tab2

# Self-esteem
prop_self <- med.out1$d0 / med.out1$tau.coef
prop_self

# Rumination
prop_rum <- med.out2$d0 / med.out2$tau.coef
prop_rum
