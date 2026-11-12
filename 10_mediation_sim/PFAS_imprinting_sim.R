# =============================================================================
# Stage 10 | PFAS_imprinting_sim.R     one array element = one simulated study
# Null-scenario simulation (n = 100 samples, 1000 isoforms) of the planned
# PFAS -> imprinted isoform -> outcome (SRS) mediation analysis; compares the
# type-I error of SEM (lavaan, bootstrap), CoCA, PGC and an adjusted linear
# model. Run r uses seed r.
# Out: SIM_DIR/run_<r>.RData  (combined_type1_errors, sem_df, coca_df, pgc_df,
#                              adjusted_lm_results)
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lavaan); library(boot); library(parallel)
})
run_id <- as.integer(Sys.getenv("IDX", "1"))
hash <- sprintf("run_%03d", run_id)
RNGkind("L'Ecuyer-CMRG")
set.seed(run_id)

# Simulate Dataset

## Define Samples and Isoforms

n_samples <- 100
n_isoforms <- 1000
simulation_scenario <- "null"
n_boot <- 1000

# Fetal sex (0 = female, 1 = male)
fetal_sex <- rbinom(n_samples, size = 1, prob = 0.5)

# Maternal PFAS exposure (ng/mL, log-normal)
# Common sdlog values for environmental contaminants in humans can range from 0.8 to 1.5.
PFAS <- rlnorm(n_samples, meanlog = log(3), sdlog = 0.9)

# Latent confounder U, correlated with PFAS
U <- rnorm(n_samples, mean = scale(PFAS), sd = 1)

## Simulate CpG Methylation from Shared Confounding Signal

n_cpg_per_isoform <- sample(1:5, size = n_isoforms, replace = TRUE)
maternal_meth_list <- lapply(n_cpg_per_isoform, function(n_cpg) {
  matrix(
    plogis(matrix(rnorm(n_samples * n_cpg, mean = rep(U, each = n_cpg), sd = 0.5), nrow = n_samples)),
    nrow = n_samples
  )
})

## Simulate Isoform-Level Expression Ratios

# Start with baseline biallelic expression + noise (maternal expression ratio)
# Using Beta distribution for ratios (bounded 0-1), which is more appropriate for proportions.
# shape1 = 50, shape2 = 50 gives a mean of 0.5 and a reasonable, but not extremely narrow, variance.
expr_ratios <- matrix(rbeta(n_samples * n_isoforms, shape1 = 50, shape2 = 50), nrow = n_samples)

# Define 50% of isoforms as true mediators
mediated_isoforms <- sample(1:n_isoforms, size = n_isoforms * 0.5)

# PFAS influences imprinting in these isoforms
for (i in mediated_isoforms) {
  effect_mean <- 0.05 # mean effect
  if (simulation_scenario == "null") {
    effect_mean <- 0 # Set mean effect to 0 for null scenario
  }
  effect <- scale(PFAS) * rnorm(1, mean = effect_mean, sd = 0.01) # Use the adjusted mean
  expr_ratios[, i] <- expr_ratios[, i] + effect
  
  # Clamp values to ensure they stay within [0,1]
  expr_ratios[, i] <- pmax(0, pmin(1, expr_ratios[, i])) 
}

## Simulate SRS Scores Influenced by Disrupted Imprinting

# Create SRS score from expression of mediated isoforms + confounder + noise
srs_effect_mean <- 2 # mean effect
if (simulation_scenario == "null") {
  srs_effect_mean <- 0 # Set mean effect to 0 for null scenario
}
SRS_contrib <- rowSums(scale(expr_ratios[, mediated_isoforms]) %*%
                         matrix(rnorm(length(mediated_isoforms), mean = srs_effect_mean, sd = 0.5), ncol = 1)) # Use the adjusted mean
SRS <- scale(SRS_contrib + scale(U) * 0.5 + rnorm(n_samples, 0, 1))

## Include Confounding Effects

# Select 50% of mediated isoforms to be confounded by methylation
confounded_idx <- sample(mediated_isoforms, size = length(mediated_isoforms) * 0.5)
confounding_contrib <- rep(0, n_samples)

for (i in confounded_idx) {
  meth_mean <- rowMeans(maternal_meth_list[[i]])
  effect_meth_expr <- scale(meth_mean) * rnorm(1, mean = 0.2, sd = 0.05)
  expr_ratios[, i] <- expr_ratios[, i] + effect_meth_expr
  expr_ratios[, i] <- pmax(0, pmin(1, expr_ratios[, i]))
  
  confounding_contrib <- confounding_contrib + scale(meth_mean) * rnorm(1, mean = 1.5, sd = 0.3)
}

# Inject confounding into SRS
SRS <- scale(SRS + confounding_contrib)

## Finalize and Format Output

expr_df <- as.data.frame(expr_ratios)
colnames(expr_df) <- paste0("isoform_", 1:n_isoforms)

metadata_df <- data.frame(
  sample = 1:n_samples,
  SRS = SRS,
  PFAS = PFAS,
  fetal_sex = fetal_sex
)

meth_long <- do.call(rbind, lapply(1:n_isoforms, function(i) {
  meth_mat <- maternal_meth_list[[i]]
  n_cpg <- ncol(meth_mat)
  data.frame(
    sample = rep(1:n_samples, times = n_cpg),
    isoform = paste0("isoform_", i),
    cpg = rep(1:n_cpg, each = n_samples),
    meth = as.vector(meth_mat)
  )
}))



# Parallelized Analysis Across Isoforms

# Predefine names for output list
isoform_names <- paste0("isoform_", 1:n_isoforms)

# Parallel processing (forked workers; reproducible for a given NCPU)
results_list <- mclapply(seq_len(n_isoforms), mc.cores = NCPU, FUN = function(i) {
  iso_name <- isoform_names[i]
  
  df <- cbind(metadata_df, expr = expr_df[[iso_name]])
  
  meth_avg <- meth_long %>%
    filter(isoform == iso_name) %>%
    group_by(sample) %>%
    summarise(m_meth = mean(meth), .groups = "drop")
  
  df <- merge(df, meth_avg, by = "sample")
  
  # Lavaan model
  model <- '
    expr ~ a1*PFAS + a2*fetal_sex + a3*PFAS:fetal_sex
    SRS ~ b1*expr + b2*fetal_sex + b3*expr:fetal_sex + c*PFAS
    ind_female := a1 * b1
    ind_male := (a1 + a3) * (b1 + b3)
    dir := c
    total_female := dir + ind_female
    total_male := dir + ind_male
  '
  
  lavaan_fit <- tryCatch({
    sem(model, data = df, se = "bootstrap", bootstrap = n_boot)
  }, error = function(e) NULL)
  
  # COCA
  coca_model <- lm(m_meth ~ PFAS + SRS * fetal_sex, data = df)
  beta_A <- coef(coca_model)["PFAS"]
  beta_Y <- coef(coca_model)["SRS"]
  psi_coca <- -beta_A / beta_Y
  
  coca_fn <- function(data, indices) {
    d <- data[indices, ]
    fit <- lm(m_meth ~ PFAS + SRS * fetal_sex, data = d)
    -coef(fit)["PFAS"] / coef(fit)["SRS"]
  }
  
  coca_boot <- boot(df, coca_fn, R = n_boot)
  coca_ci <- tryCatch(boot.ci(coca_boot, type = "perc"), error = function(e) NULL)
  
  # PGC
  bridge_model <- lm(m_meth ~ expr + SRS * fetal_sex, data = df)
  df$Z_hat <- predict(bridge_model)
  df$Y_adj <- df$SRS - (df$Z_hat - mean(df$Z_hat))
  
  pgc_model <- lm(Y_adj ~ PFAS * fetal_sex, data = df)
  
  pgc_fn <- function(data, indices) {
    d <- data[indices, ]
    bridge <- lm(m_meth ~ expr + SRS * fetal_sex, data = d)
    Z_hat <- predict(bridge)
    Y_adj <- d$SRS - (Z_hat - mean(Z_hat))
    pgc_fit <- lm(Y_adj ~ PFAS * fetal_sex, data = cbind(d, Y_adj))
    coef(pgc_fit)["PFAS"]
  }
  
  pgc_boot <- boot(df, pgc_fn, R = n_boot)
  pgc_ci <- tryCatch(boot.ci(pgc_boot, type = "perc"), error = function(e) NULL)
  
  # Return result for this isoform
  list(
    lavaan_fit = lavaan_fit,
    COCA = list(
      model = coca_model,
      boot = coca_boot,
      ci = coca_ci
    ),
    PGC = list(
      model = pgc_model,
      boot = pgc_boot,
      ci = pgc_ci
    )
  )
})
bad <- vapply(results_list, inherits, logical(1), "try-error")
if (any(bad)) stop(sum(bad), " isoforms failed: ", results_list[[which(bad)[1]]])

# Set names to match isoform identifiers
names(results_list) <- isoform_names


# Initialize results data.frames

sem_df <- data.frame()
sem_pval_df <- data.frame()
coca_df <- data.frame()
pgc_df <- data.frame()

iso_names <- names(results_list)

for (i in seq_along(results_list)) {
  iso <- iso_names[i]
  res <- results_list[[i]]
  
  # --- SEM path p-values
  fit <- res$lavaan_fit
  if (!is.null(fit)) {
    pe <- tryCatch(parameterEstimates(fit), error = function(e) NULL)
    if (!is.null(pe)) {
      selected_paths <- pe %>%
        filter(label %in% c("a1", "a3", "b1", "b3", "c",
                            "ind_female", "ind_male", "total_female", "total_male"))
      sem_df_out <- selected_paths
      sem_df_out$isoform <- iso
      sem_df <- bind_rows(sem_df,sem_df_out)
      selected_paths <- selected_paths %>% select(label, pvalue)
      if (nrow(selected_paths) > 0) {
        sem_pval_df <- bind_rows(sem_pval_df,
                                 data.frame(
                                   isoform = iso,
                                   path = selected_paths$label,
                                   pval = selected_paths$pvalue
                                 ))
      }
    }
  }
  
  # --- COCA model: extract PFAS, fetal_sex, and interaction p-values and CI
  coca_fit <- results_list[[i]][["COCA"]][["model"]]
  
  if (!is.null(coca_fit)) {
    coef_table <- coef(summary(coca_fit))
    confint_table <- tryCatch(confint(coca_fit), error = function(e) NULL)
    
    pval_pfas  <- coef_table["PFAS", "Pr(>|t|)"]
    pval_sex   <- if ("fetal_sex" %in% rownames(coef_table)) coef_table["fetal_sex", "Pr(>|t|)"] else NA
    pval_int   <- if ("SRS:fetal_sex" %in% rownames(coef_table)) coef_table["SRS:fetal_sex", "Pr(>|t|)"] else NA
    
    ci_lower <- ci_upper <- NA
    if (!is.null(confint_table) && "PFAS" %in% rownames(confint_table)) {
      ci_lower <- confint_table["PFAS", 1]
      ci_upper <- confint_table["PFAS", 2]
    }
    
    signif <- !(is.na(ci_lower) || (ci_lower < 0 & ci_upper > 0))  # TRUE if CI excludes 0
    
    coca_df <- bind_rows(coca_df, data.frame(
      isoform = iso,
      pval_PFAS = pval_pfas,
      pval_fetal_sex = pval_sex,
      pval_PFAS_fetal_sex = pval_int,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      significant = signif
    ))
  }
  
  # --- PGC model: extract PFAS, fetal_sex, and interaction p-values
  pgc_model <- tryCatch(res$PGC$model, error = function(e) NULL)
  pgc_boot  <- tryCatch(res$PGC$boot, error = function(e) NULL)
  
  if (!is.null(pgc_model) && !is.null(pgc_boot)) {
    coef_pg <- coef(summary(pgc_model))
    
    pval_pfas <- coef_pg["PFAS", "Pr(>|t|)"]
    pval_sex  <- if ("fetal_sex" %in% rownames(coef_pg)) coef_pg["fetal_sex", "Pr(>|t|)"] else NA
    pval_int  <- if ("PFAS:fetal_sex" %in% rownames(coef_pg)) coef_pg["PFAS:fetal_sex", "Pr(>|t|)"] else NA
    
    boot_estimates <- pgc_boot$t
    signif <- NA
    ci_lower <- ci_upper <- NA
    if (!is.null(boot_estimates) && length(boot_estimates) > 1) {
      ci_bounds <- quantile(boot_estimates, probs = c(0.025, 0.975), na.rm = TRUE)
      ci_lower <- ci_bounds[1]
      ci_upper <- ci_bounds[2]
      signif <- !(ci_lower < 0 & ci_upper > 0)
    }
    
    pgc_df <- bind_rows(pgc_df, data.frame(
      isoform = iso,
      pval_PFAS = pval_pfas,
      pval_fetal_sex = pval_sex,
      pval_PFAS_fetal_sex = pval_int,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      significant = signif
    ))
  }
}


# Regression adjustment models with control outcome variable as covariate

adjusted_lm_results <- data.frame(
  isoform = 1:n_isoforms,
  coef_expr = NA,           p_expr = NA,           sig_ci_expr = NA,
  coef_sex = NA,            p_sex = NA,            sig_ci_sex = NA,
  coef_interaction = NA,    p_interaction = NA,    sig_ci_interaction = NA
)

for (i in 1:n_isoforms) {
  meth_proxy <- rowMeans(maternal_meth_list[[i]])
  fit <- lm(SRS ~ PFAS + expr_ratios[, i] * fetal_sex + meth_proxy)
  confint_fit <- confint(fit)
  coef_names <- names(coef(fit))
  
  expr_term <- "expr_ratios[, i]"
  sex_term  <- "fetal_sex"
  int_term  <- "expr_ratios[, i]:fetal_sex"
  
  # Expression term
  if (expr_term %in% coef_names) {
    idx <- which(coef_names == expr_term)
    adjusted_lm_results$coef_expr[i] <- coef(fit)[idx]
    adjusted_lm_results$p_expr[i] <- summary(fit)$coefficients[idx, "Pr(>|t|)"]
    adjusted_lm_results$sig_ci_expr[i] <- !(0 >= confint_fit[idx, 1] & 0 <= confint_fit[idx, 2])
  }
  
  # Fetal sex term
  if (sex_term %in% coef_names) {
    idx <- which(coef_names == sex_term)
    adjusted_lm_results$coef_sex[i] <- coef(fit)[idx]
    adjusted_lm_results$p_sex[i] <- summary(fit)$coefficients[idx, "Pr(>|t|)"]
    adjusted_lm_results$sig_ci_sex[i] <- !(0 >= confint_fit[idx, 1] & 0 <= confint_fit[idx, 2])
  }
  
  # Interaction term
  if (int_term %in% coef_names) {
    idx <- which(coef_names == int_term)
    adjusted_lm_results$coef_interaction[i] <- coef(fit)[idx]
    adjusted_lm_results$p_interaction[i] <- summary(fit)$coefficients[idx, "Pr(>|t|)"]
    adjusted_lm_results$sig_ci_interaction[i] <- !(0 >= confint_fit[idx, 1] & 0 <= confint_fit[idx, 2])
  }
}





# For a null simulation, we expect these rates to be close to a = 0.05
# for true null effects, indicating correct Type I error control.
# Rates significantly above 0.05 suggest inflated Type I error.

type1_error_results <- list()

# 1. SEM Type I Error
# Focusing on the mediation paths (a1, b1) and indirect effects (ind_female, ind_male)
# which are set to null in the simulation scenario.
sem_type1 <- sem_pval_df %>%
  filter(path %in% c("a1", "b1", "ind_female", "ind_male")) %>%
  group_by(path) %>%
  summarise(
    Type1_Error_Rate = mean(pval < 0.05, na.rm = TRUE),
    .groups = "drop"
  )
type1_error_results$SEM <- sem_type1

# 2. COCA Type I Error
# For COCA, 'significant' is determined by whether the 95% CI for the PFAS effect excludes 0.
# For other terms (fetal_sex, SRS:fetal_sex interaction), we use a p-value threshold.
coca_type1 <- coca_df %>%
  summarise(
    Type1_Error_Rate_PFAS = mean(significant, na.rm = TRUE), # Uses the 'significant' flag from CI
    Type1_Error_Rate_fetal_sex = mean(pval_fetal_sex < 0.05, na.rm = TRUE),
    Type1_Error_Rate_SRS_fetal_sex_int = mean(pval_PFAS_fetal_sex < 0.05, na.rm = TRUE) # This is for the SRS:fetal_sex interaction in the COCA model
  ) %>%
  pivot_longer(everything(), names_to = "path", values_to = "Type1_Error_Rate") %>%
  mutate(path = gsub("Type1_Error_Rate_", "", path))
type1_error_results$COCA <- coca_type1

# 3. PGC Type I Error
# For PGC, 'significant' for PFAS is also based on CI exclusion of 0.
# For other terms (fetal_sex, PFAS:fetal_sex interaction), p-value threshold is used.
pgc_type1 <- pgc_df %>%
  summarise(
    Type1_Error_Rate_PFAS = mean(significant, na.rm = TRUE), # Uses the 'significant' flag from CI
    Type1_Error_Rate_fetal_sex = mean(pval_fetal_sex < 0.05, na.rm = TRUE),
    Type1_Error_Rate_PFAS_fetal_sex_int = mean(pval_PFAS_fetal_sex < 0.05, na.rm = TRUE) # This is for the PFAS:fetal_sex interaction in the PGC model
  ) %>%
  pivot_longer(everything(), names_to = "path", values_to = "Type1_Error_Rate") %>%
  mutate(path = gsub("Type1_Error_Rate_", "", path))
type1_error_results$PGC <- pgc_type1

# 4. Regression Adjustment Type I Error
# 'sig_ci_*' flags are based on CI excluding 0.
adj_type1 <- adjusted_lm_results %>%
  summarise(
    Type1_Error_Rate_expr = mean(sig_ci_expr, na.rm = TRUE),
    Type1_Error_Rate_sex = mean(sig_ci_sex, na.rm = TRUE),
    Type1_Error_Rate_expr_fetal_sex_int = mean(sig_ci_interaction, na.rm = TRUE) # This is for the expr:fetal_sex interaction in the direct adjustment model
  ) %>%
  pivot_longer(everything(), names_to = "path", values_to = "Type1_Error_Rate") %>%
  mutate(path = gsub("Type1_Error_Rate_", "", path))
type1_error_results$DirectAdj <- adj_type1


# Combine all tibbles into a single data frame
combined_type1_errors <- bind_rows(
  type1_error_results$SEM %>% mutate(Method = "SEM"),
  type1_error_results$COCA %>% mutate(Method = "COCA"),
  type1_error_results$PGC %>% mutate(Method = "PGC"),
  type1_error_results$DirectAdj %>% mutate(Method = "Direct Adjustment")
) %>%
  # Reorder columns for better presentation
  select(Method, path, Type1_Error_Rate) %>%
  # Rename paths for clarity if desired (e.g., "a1" to "PFAS -> Expression")
  mutate(
    path = case_when(
      path == "a1" ~ "PFAS -> Expression (a1)",
      path == "b1" ~ "Expression -> SRS (b1)",
      path == "ind_female" ~ "Indirect Effect (Female)",
      path == "ind_male" ~ "Indirect Effect (Male)",
      path == "PFAS" ~ "PFAS Main Effect",
      path == "fetal_sex" ~ "Fetal Sex Main Effect",
      path == "SRS_fetal_sex_int" ~ "SRS:Fetal Sex Interaction",
      path == "PFAS_fetal_sex_int" ~ "PFAS:Fetal Sex Interaction",
      path == "expr" ~ "Expression Main Effect",
      path == "sex" ~ "Sex Main Effect",
      path == "expr_fetal_sex_int" ~ "Expression:Fetal Sex Interaction",
      TRUE ~ path # Keep others as is
    )
  )


# Prepare output for writing

combined_type1_errors <- data.frame(combined_type1_errors)
combined_type1_errors$run <- hash
sem_df$run <- hash
coca_df$run <- hash
row.names(pgc_df) <- iso_names
row.names(pgc_df) <- NULL
pgc_df$run <- hash
adjusted_lm_results$run <- hash


# Save results to file
dir.create(P("SIM_DIR"), recursive = TRUE, showWarnings = FALSE)
save(list = c("combined_type1_errors", "sem_df", "coca_df", "pgc_df", "adjusted_lm_results"),
     file = file.path(P("SIM_DIR"), paste0(hash, ".RData")))
