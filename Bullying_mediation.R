# Load libraries
library(dplyr)
library(lm.beta)
library(mediation)

# Read in data
EDData <- readRDS("~/Desktop/ED Paper/ED_fbsv_rses.rds")

# Age as a predictor
mod_age_edeqs <- lm(edeqs_total ~ age, data = EDData)
mod_age_fbsv <- lm(fbsv_total ~ age, data = EDData)
mod_age_ptq <- lm(ptq_total ~ age, data = EDData)
mod_age_rses <- lm(rses_sc_score ~ age, data = EDData)

# Print model results
summary(mod_age_edeqs)
summary(mod_age_fbsv)
summary(mod_age_ptq)
summary(mod_age_rses)

# Standardized beta
lm.beta(mod_age_edeqs)
lm.beta(mod_age_fbsv)
lm.beta(mod_age_ptq)
lm.beta(mod_age_rses)

## Mediations

# ----------------------------
# Helper functions
# ----------------------------

set.seed(1234)

# Significance code
sig_code <- function(p_value) {
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

# ----------------------------
# Extract path info
# ----------------------------

extract_path <- function(model, var_name) {
  
  coef_info <- summary(model)$coefficients[var_name, ]
  
  estimate <- coef_info["Estimate"]
  se <- coef_info["Std. Error"]
  t_val <- estimate / se
  p_val <- coef_info["Pr(>|t|)"]
  ci <- estimate + c(-1.96, 1.96) * se
  
  data.frame(
    Beta = round(estimate, 3),
    SE = round(se, 3),
    Lower = round(ci[1], 3),
    Upper = round(ci[2], 3),
    t = round(t_val, 3),
    p = format_pval(p_val),
    Signif = sig_code(p_val)
  )
}

# ----------------------------
# Full table function
# ----------------------------

make_full_table <- function(model.m, model.y, med.out, treat, mediator) {
  
  # Path coefficients
  path_a <- extract_path(model.m, treat)
  path_b <- extract_path(model.y, mediator)
  path_y <- extract_path(model.y, treat)
  
  s <- summary(med.out)
  
  # -------------------
  # Indirect effect
  # -------------------
  
  est_ind <- s$d0
  se_ind  <- sd(s$d0.sims)
  t_ind   <- est_ind / se_ind
  p_ind   <- s$d0.p
  ci_ind  <- s$d0.ci
  
  # -------------------
  # Direct effect
  # -------------------
  
  est_ade <- s$z0
  se_ade  <- sd(s$z0.sims)
  t_ade   <- est_ade / se_ade
  p_ade   <- s$z0.p
  ci_ade  <- s$z0.ci
  
  # -------------------
  # Total effect
  # -------------------
  
  est_tot <- s$tau.coef
  se_tot  <- sd(s$tau.sims)
  t_tot   <- est_tot / se_tot
  p_tot   <- s$tau.p
  ci_tot  <- quantile(s$tau.sims, c(.025, .975))
  
  # -------------------
  # Combine table
  # -------------------
  
  data.frame(
    Effect = c(
      "Indirect (ACME)",
      "Component a",
      "Component b",
      "ADE",
      "Total"
    ),
    
    Beta = round(
      c(
        est_ind,
        path_a$Beta,
        path_b$Beta,
        est_ade,
        est_tot
      ),
      3
    ),
    
    SE = round(
      c(
        se_ind,
        path_a$SE,
        path_b$SE,
        se_ade,
        se_tot
      ),
      3
    ),
    
    Lower = round(
      c(
        ci_ind[1],
        path_a$Lower,
        path_b$Lower,
        ci_ade[1],
        ci_tot[1]
      ),
      3
    ),
    
    Upper = round(
      c(
        ci_ind[2],
        path_a$Upper,
        path_b$Upper,
        ci_ade[2],
        ci_tot[2]
      ),
      3
    ),
    
    t = round(
      c(
        t_ind,
        path_a$t,
        path_b$t,
        t_ade,
        t_tot
      ),
      3
    ),
    
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
    )
  )
}

# ----------------------------
# Standardize variables
# ----------------------------

EDData_std <- as.data.frame(
  scale(
    EDData[, c(
      "fbsv_total",
      "rses_sc_score",
      "ptq_total",
      "edeqs_total",
      "age"
    )]
  )
)

# ----------------------------
# Fit standardized models
# ----------------------------

# Mediator models
model.m1.std <- lm(
  rses_sc_score ~ fbsv_total + age,
  data = EDData_std
)

model.m2.std <- lm(
  ptq_total ~ fbsv_total + age,
  data = EDData_std
)

# Outcome model
model.y.std <- lm(
  edeqs_total ~ fbsv_total + rses_sc_score + ptq_total + age,
  data = EDData_std
)

# ----------------------------
# Mediation analyses
# ----------------------------

med.out1.std <- mediate(
  model.m1.std,
  model.y.std,
  treat = "fbsv_total",
  mediator = "rses_sc_score",
  boot = TRUE,
  sims = 5000
)

med.out2.std <- mediate(
  model.m2.std,
  model.y.std,
  treat = "fbsv_total",
  mediator = "ptq_total",
  boot = TRUE,
  sims = 5000
)

# ----------------------------
# Generate tables
# ----------------------------

tab1 <- make_full_table(
  model.m1.std,
  model.y.std,
  med.out1.std,
  "fbsv_total",
  "rses_sc_score"
)

tab2 <- make_full_table(
  model.m2.std,
  model.y.std,
  med.out2.std,
  "fbsv_total",
  "ptq_total"
)

tab1
tab2

##Paired Bootstrap Comparison
set.seed(1234)

# Number of bootstrap samples
B <- 5000

# Store bootstrap estimates
boot_ind1 <- numeric(B)
boot_ind2 <- numeric(B)
boot_diff <- numeric(B)

n <- nrow(EDData)

for (b in 1:B) {
  
  # ---------------------------------
  # Bootstrap resample
  # ---------------------------------
  idx <- sample(1:n, size = n, replace = TRUE)
  dat_b <- EDData[idx, ]
  
  # ---------------------------------
  # Mediator models
  # ---------------------------------
  m1_b <- lm(rses_sc_score ~ fbsv_total + age, data = dat_b)
  m2_b <- lm(ptq_total ~ fbsv_total + age, data = dat_b)
  
  # ---------------------------------
  # Outcome model
  # ---------------------------------
  y_b <- lm(
    edeqs_total ~ fbsv_total + rses_sc_score + ptq_total + age,
    data = dat_b
  )
  
  # ---------------------------------
  # Extract path coefficients
  # ---------------------------------
  
  # a paths
  a1 <- coef(m1_b)["fbsv_total"]
  a2 <- coef(m2_b)["fbsv_total"]
  
  # b paths
  b1 <- coef(y_b)["rses_sc_score"]
  b2 <- coef(y_b)["ptq_total"]
  
  # ---------------------------------
  # Indirect effects
  # ---------------------------------
  ind1 <- a1 * b1
  ind2 <- a2 * b2
  
  # Store
  boot_ind1[b] <- ind1
  boot_ind2[b] <- ind2
  boot_diff[b] <- ind1 - ind2
}

# ---------------------------------
# Point estimates
# ---------------------------------

ind1_est <- mean(boot_ind1)
ind2_est <- mean(boot_ind2)
diff_est <- mean(boot_diff)

# ---------------------------------
# Confidence intervals
# ---------------------------------

ind1_ci <- quantile(boot_ind1, c(.025, .975))
ind2_ci <- quantile(boot_ind2, c(.025, .975))
diff_ci <- quantile(boot_diff, c(.025, .975))

# ---------------------------------
# Two-tailed p-value
# ---------------------------------

p_diff <- 2 * min(
  mean(boot_diff >= 0),
  mean(boot_diff <= 0)
)

# ---------------------------------
# Results
# ---------------------------------

diff_est
diff_ci
p_diff



# ----------------------------
# Sensitivity Analyses
# ----------------------------

set.seed(1234)
sims <- 5000

# ----------------------------
# Standardize variables
# ----------------------------

vars <- c(
  "fbsv_total",
  "rses_sc_score",
  "ptq_total",
  "edeqs_total",
  "age"
)

EDData_std <- as.data.frame(scale(EDData[, vars]))

# ----------------------------
# Helper functions
# ----------------------------

# Significance stars
signif_code <- function(p) {
  if (is.na(p)) return("")
  if (p < .001) return("***")
  if (p < .01)  return("**")
  if (p < .05)  return("*")
  if (p < .10)  return(".")
  return("")
}

# Format p-values
format_p <- function(p) {
  if (is.na(p)) return(NA)
  if (p < .001) return("<.001")
  return(round(p, 3))
}

# ----------------------------
# Extract mediation results
# ----------------------------

extract_mediation <- function(
    med_model,
    out_model,
    med_obj,
    treat,
    mediator
) {
  
  s <- summary(med_obj)
  
  # ------------------------
  # a path
  # ------------------------
  
  a <- coef(summary(med_model))[treat, ]
  
  a_est <- a["Estimate"]
  a_se  <- a["Std. Error"]
  a_p   <- a["Pr(>|t|)"]
  a_ci  <- a_est + c(-1.96, 1.96) * a_se
  
  # ------------------------
  # b path
  # ------------------------
  
  b <- coef(summary(out_model))[mediator, ]
  
  b_est <- b["Estimate"]
  b_se  <- b["Std. Error"]
  b_p   <- b["Pr(>|t|)"]
  b_ci  <- b_est + c(-1.96, 1.96) * b_se
  
  # ------------------------
  # Indirect effect (ACME)
  # ------------------------
  
  ind_est <- s$d0
  ind_se  <- sd(s$d0.sims)
  ind_ci  <- s$d0.ci
  ind_p   <- s$d0.p
  
  # ------------------------
  # Direct effect (ADE)
  # ------------------------
  
  ade_est <- s$z0
  ade_se  <- sd(s$z0.sims)
  ade_ci  <- s$z0.ci
  ade_p   <- s$z0.p
  
  # ------------------------
  # Total effect
  # ------------------------
  
  tot_est <- s$tau.coef
  tot_se  <- sd(s$tau.sims)
  tot_ci  <- quantile(s$tau.sims, c(.025, .975))
  tot_p   <- s$tau.p
  
  # ------------------------
  # Combine table
  # ------------------------
  
  data.frame(
    Effect = c(
      "Indirect (ACME)",
      "a path",
      "b path",
      "ADE",
      "Total"
    ),
    
    Beta = round(
      c(ind_est, a_est, b_est, ade_est, tot_est),
      3
    ),
    
    SE = round(
      c(ind_se, a_se, b_se, ade_se, tot_se),
      3
    ),
    
    Lower = round(
      c(ind_ci[1], a_ci[1], b_ci[1], ade_ci[1], tot_ci[1]),
      3
    ),
    
    Upper = round(
      c(ind_ci[2], a_ci[2], b_ci[2], ade_ci[2], tot_ci[2]),
      3
    ),
    
    p = c(
      format_p(ind_p),
      format_p(a_p),
      format_p(b_p),
      format_p(ade_p),
      format_p(tot_p)
    ),
    
    Signif = c(
      signif_code(ind_p),
      signif_code(a_p),
      signif_code(b_p),
      signif_code(ade_p),
      signif_code(tot_p)
    )
  )
}

# =========================================================
# MODEL 1:
# fbsv_total → rses_sc_score → edeqs_total
# =========================================================

# Mediator model
med_model1_std <- lm(
  rses_sc_score ~ fbsv_total + age,
  data = EDData_std
)

# Outcome model
out_model1_std <- lm(
  edeqs_total ~ fbsv_total + rses_sc_score + age,
  data = EDData_std
)

# Mediation analysis
med1_std <- mediate(
  med_model1_std,
  out_model1_std,
  treat = "fbsv_total",
  mediator = "rses_sc_score",
  boot = TRUE,
  sims = sims
)

# Results table
tab1 <- extract_mediation(
  med_model1_std,
  out_model1_std,
  med1_std,
  treat = "fbsv_total",
  mediator = "rses_sc_score"
)

tab1

# =========================================================
# MODEL 2:
# fbsv_total → ptq_total → edeqs_total
# =========================================================

# Mediator model
med_model2_std <- lm(
  ptq_total ~ fbsv_total + age,
  data = EDData_std
)

# Outcome model
out_model2_std <- lm(
  edeqs_total ~ fbsv_total + ptq_total + age,
  data = EDData_std
)

# Mediation analysis
med2_std <- mediate(
  med_model2_std,
  out_model2_std,
  treat = "fbsv_total",
  mediator = "ptq_total",
  boot = TRUE,
  sims = sims
)

# Results table
tab2 <- extract_mediation(
  med_model2_std,
  out_model2_std,
  med2_std,
  treat = "fbsv_total",
  mediator = "ptq_total"
)

tab2


### FBSV Log-Transformed Analysis ###

EDData <- EDData %>%
  mutate(fbsv_log = log(fbsv_total + 1))

# ----------------------------
# Helper functions
# ----------------------------

set.seed(1234)

# Significance code
sig_code <- function(p_value) {
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

# Extract path info from lm model
extract_path <- function(model, var_name) {
  
  coef_info <- summary(model)$coefficients[var_name, ]
  
  estimate <- coef_info["Estimate"]
  se <- coef_info["Std. Error"]
  p_val <- coef_info["Pr(>|t|)"]
  ci <- confint(model, parm = var_name)
  
  data.frame(
    SE = round(se, 3),
    Lower = round(ci[1], 3),
    Upper = round(ci[2], 3),
    B = round(estimate, 3),
    p = format_pval(p_val),
    Signif = sig_code(p_val)
  )
}

# ----------------------------
# Full table function
# ----------------------------

make_full_table <- function(model.m, model.y, med.out, treat, mediator) {
  
  # Path coefficients
  path_a <- extract_path(model.m, treat)
  path_b <- extract_path(model.y, mediator)
  
  s <- summary(med.out)
  
  # -------------------
  # Indirect (ACME)
  # -------------------
  est_ind <- if (!is.null(s$d0)) s$d0 else NA
  se_ind  <- if (!is.null(s$d0.sims)) sd(s$d0.sims) else NA
  p_ind <- s$d0.p
  ci_ind  <- if (!is.null(s$d0.ci)) s$d0.ci else c(NA, NA)
  beta_ind <- est_ind   # already standardized
  
  # -------------------
  # ADE (direct effect)
  # -------------------
  est_ade <- if (!is.null(s$z0)) s$z0 else NA
  se_ade  <- if (!is.null(s$z0.sims)) sd(s$z0.sims) else NA
  p_ade <- s$z0.p
  ci_ade  <- if (!is.null(s$z0.ci)) s$z0.ci else c(NA, NA)
  beta_ade <- est_ade   # already standardized
  
  # -------------------
  # Total effect
  # -------------------
  est_tot <- if (!is.null(s$tau.coef)) s$tau.coef else NA
  se_tot  <- if (!is.null(s$tau.sims)) sd(s$tau.sims) else NA
  p_tot <- s$tau.p
  ci_tot <- if (!is.null(s$tau.ci)) s$tau.ci else c(NA, NA)
  beta_tot <- est_tot   # already standardized
  
  # -------------------
  # Combine into table
  # -------------------
  
  data.frame(
    Effect = c("Indirect (ACME)", "Component a", "Component b", "ADE", "Total"),
    
    B = round(
      c(beta_ind, path_a$B, path_b$B, beta_ade, beta_tot), 3
    ),
    
    SE = round(
      c(se_ind, path_a$SE, path_b$SE, se_ade, se_tot), 3
    ),
    
    Lower = round(
      c(ci_ind[1], path_a$Lower, path_b$Lower, ci_ade[1], ci_tot[1]), 3
    ),
    
    Upper = round(
      c(ci_ind[2], path_a$Upper, path_b$Upper, ci_ade[2], ci_tot[2]), 3
    ),
    
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
    )
  )
}

# ----------------------------
# Standardize variables
# ----------------------------

EDData_log_std <- as.data.frame(
  scale(
    EDData[, c(
      "fbsv_log",
      "rses_sc_score",
      "ptq_total",
      "edeqs_total",
      "age"
    )]
  )
)

# ----------------------------
# Fit standardized models
# ----------------------------

# Mediator models
logmodel.m1.std <- lm(
  rses_sc_score ~ fbsv_log + age,
  data = EDData_log_std
)

logmodel.m2.std <- lm(
  ptq_total ~ fbsv_log + age,
  data = EDData_log_std
)

# Outcome model
logmodel.y.std <- lm(
  edeqs_total ~ fbsv_log + rses_sc_score + ptq_total + age,
  data = EDData_log_std
)

# ----------------------------
# Mediation analyses
# ----------------------------

logmed.out1.std <- mediate(
  logmodel.m1.std,
  logmodel.y.std,
  treat = "fbsv_log",
  mediator = "rses_sc_score",
  boot = TRUE,
  sims = 5000
)

logmed.out2.std <- mediate(
  logmodel.m2.std,
  logmodel.y.std,
  treat = "fbsv_log",
  mediator = "ptq_total",
  boot = TRUE,
  sims = 5000
)

# ----------------------------
# Generate tables
# ----------------------------

logtab1 <- make_full_table(
  model.m = logmodel.m1.std,
  model.y = logmodel.y.std,
  med.out = logmed.out1.std,
  treat = "fbsv_log",
  mediator = "rses_sc_score"
)

logtab2 <- make_full_table(
  model.m = logmodel.m2.std,
  model.y = logmodel.y.std,
  med.out = logmed.out2.std,
  treat = "fbsv_log",
  mediator = "ptq_total"
)

# View results
logtab1
logtab2



# Influence Diagnostics
# Self-esteem mediator
m_mod <- lm(rses_sc_score ~ fbsv_total + ptq_total + age, data = EDData)
y_mod <- lm(edeqs_total ~ fbsv_total + rses_sc_score + ptq_total + age, data = EDData)

plot(m_mod, which = 4)  # Cook's distance
plot(y_mod, which = 4)

which(cooks.distance(m_mod) > 4 / nrow(EDData))
which(cooks.distance(y_mod) > 4 / nrow(EDData))

sort(cooks.distance(m_mod), decreasing = TRUE)[1:10]
sort(cooks.distance(y_mod), decreasing = TRUE)[1:10]

sum(cooks.distance(m_mod) > 0.5, na.rm = TRUE)
sum(cooks.distance(y_mod) > 0.5, na.rm = TRUE)

sum(cooks.distance(m_mod) > 1, na.rm = TRUE)
sum(cooks.distance(y_mod) > 1, na.rm = TRUE)

# Rumination mediator
m_mod2 <- lm(ptq_total ~ fbsv_total + rses_sc_score + age, data = EDData)
y_mod2 <- lm(edeqs_total ~ fbsv_total + ptq_total + rses_sc_score + age, data = EDData)

plot(m_mod2, which = 4)  # Cook's distance
plot(y_mod2, which = 4)

which(cooks.distance(m_mod2) > 4 / nrow(EDData))
which(cooks.distance(y_mod2) > 4 / nrow(EDData))

sort(cooks.distance(m_mod2), decreasing = TRUE)[1:10]
sort(cooks.distance(y_mod2), decreasing = TRUE)[1:10]

sum(cooks.distance(m_mod2) > 0.5, na.rm = TRUE)
sum(cooks.distance(y_mod2) > 0.5, na.rm = TRUE)

sum(cooks.distance(m_mod2) > 1, na.rm = TRUE)
sum(cooks.distance(y_mod2) > 1, na.rm = TRUE)