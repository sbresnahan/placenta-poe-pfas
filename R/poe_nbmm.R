# =============================================================================
# R/poe_nbmm.R -- cohort-level parent-of-origin NB-GLMMs (glmmTMB)
#
# Baseline (stage 07), per feature (gene or transcript):
#   count ~ allele + sex + ancestry + (1|SNPpos) + (1|sample_id),  offset log(size factor)
#   beta = log(paternal / maternal) expression ratio ("allelepaternal")
#
# Interaction (stage 08), per feature and condition X:
#   count ~ allele * X + ancestry + sex + (1|SNPpos) + (1|sample_id)
#   beta_int = change in the POE log-ratio per unit of X
#   continuous X : POE at X = 0 (beta_POE_at_ref), at +1 step (step = SD of X),
#                  at the mean and mean + 1 SD of the modelled samples ("subset"),
#                  and along dose grids -8..8 (delta-method 95% CIs)
#   binary X     : POE in the reference group (cond0) and the other group (cond1)
#
# Family nbinom2, nlminb (rel.tol 1e-5, 300 iterations). Counts are rounded.
# Size factors: median-of-ratios on maternal + paternal counts (all features).
# =============================================================================

POE_FIT_DEFAULTS <- list(min_snps = 3, min_samples = 15, min_total_counts = 50,
                         rel_tol = 1e-5, max_iter = 300, optimizer = "nlminb",
                         dose_range = -8:8, n_boot = 100)

# ---- matrices --------------------------------------------------------------

#' Split a POE count matrix (rows feature:SNP, cols <sample>_maternal/_paternal)
poe_split_matrix <- function(x) {
  x <- as.matrix(x); storage.mode(x) <- "double"
  rn <- rownames(x)
  nm <- colnames(x)
  samp <- sub("_(maternal|paternal)$", "", nm)
  samples <- intersect(samp[grepl("_maternal$", nm)], samp[grepl("_paternal$", nm)])
  M <- x[, paste0(samples, "_maternal"), drop = FALSE]
  P <- x[, paste0(samples, "_paternal"), drop = FALSE]
  colnames(M) <- colnames(P) <- samples
  list(M = M, P = P, samples = samples,
       feature_id = sub("^(.*):[^:]+$", "\\1", rn),
       snppos     = sub("^.*:([^:]+)$", "\\1", rn))
}

#' Median-of-ratios size factors from maternal + paternal counts
poe_size_factors <- function(M, P) {
  lib <- M + P
  geo <- exp(rowMeans(log(lib + 1), na.rm = TRUE)) - 1
  valid <- is.finite(geo) & geo > 0
  if (sum(valid) < 100) {
    keep <- rowSums(lib) > 0
    geo <- rep(NA_real_, nrow(lib))
    geo[keep] <- exp(rowMeans(log(lib[keep, , drop = FALSE] + 0.5), na.rm = TRUE)) - 0.5
    valid <- is.finite(geo) & geo > 0
  }
  ratios <- lib[valid, , drop = FALSE] / geo[valid]
  sf <- apply(ratios, 2, function(v) median(v[v > 0], na.rm = TRUE))
  sf <- sf / exp(mean(log(sf[sf > 0]), na.rm = TRUE))
  sf[!is.finite(sf) | sf <= 0] <- 1
  setNames(as.numeric(sf), colnames(M))
}

#' Features passing the count filters (computed on maternal + paternal counts)
poe_keep_features <- function(M, P, feature_id, min_snps, min_samples, min_total_counts) {
  lib <- M + P
  tot <- rowsum(rowSums(lib), feature_id)[, 1]
  n_samp <- rowSums(rowsum((lib > 0) * 1, feature_id) > 0)
  n_snp <- table(feature_id)[names(tot)]
  names(tot)[tot >= min_total_counts & n_samp >= min_samples & as.numeric(n_snp) >= min_snps]
}

#' Long table (one row per feature x SNP x sample x allele) for a set of rows
poe_long <- function(M, P, feature_id, snppos, size_factor, rows) {
  nR <- length(rows); S <- colnames(M)
  base <- data.table::data.table(feature_id = rep(feature_id[rows], times = length(S)),
                                 SNPpos = rep(snppos[rows], times = length(S)),
                                 sample_id = rep(S, each = nR))
  dt <- data.table::rbindlist(list(
    base[, `:=`(allele = "maternal", count = as.vector(M[rows, , drop = FALSE]))][],
    data.table::copy(base)[, `:=`(allele = "paternal", count = as.vector(P[rows, , drop = FALSE]))][]))
  dt[, offset := log(size_factor[sample_id])]
  dt[is.finite(count) & count >= 0 & is.finite(offset)][, count := round(count)][]
}

# ---- fitting ----------------------------------------------------------------

.poe_control <- function(opt) glmmTMB::glmmTMBControl(
  optimizer = match.fun(opt$optimizer), parallel = 1L,
  optCtrl = list(rel.tol = opt$rel_tol, eval.max = opt$max_iter, iter.max = opt$max_iter))

.feature_stats <- function(sub) {
  by_s <- sub[, .(n = sum(count)), by = sample_id]
  list(n_samples = sum(by_s$n > 0), n_snps = data.table::uniqueN(sub$SNPpos), total_counts = sum(sub$count))
}

.prep_factors <- function(sub) {
  sub <- data.table::copy(sub)
  sub[, `:=`(allele = factor(allele, levels = c("maternal", "paternal")),
             sample_id = factor(sample_id), SNPpos = factor(SNPpos),
             sex = factor(sex), ancestry = factor(ancestry))]
  sub
}

#' Baseline POE model for one feature (sub: long table with sex, ancestry)
fit_poe_baseline_one <- function(sub, opt = POE_FIT_DEFAULTS) {
  st <- .feature_stats(sub)
  out <- data.frame(feature_id = sub$feature_id[1], beta = NA_real_, se = NA_real_, z = NA_real_,
                    p = NA_real_, direction = NA_character_, n_samples = st$n_samples,
                    n_snps = st$n_snps, total_counts = st$total_counts, stringsAsFactors = FALSE)
  if (st$n_samples < opt$min_samples || st$total_counts < opt$min_total_counts || st$n_snps < opt$min_snps)
    return(out)
  fit <- try(glmmTMB::glmmTMB(count ~ allele + sex + ancestry + (1 | SNPpos) + (1 | sample_id),
                              family = glmmTMB::nbinom2(), offset = offset,
                              data = .prep_factors(sub), control = .poe_control(opt)), silent = TRUE)
  if (inherits(fit, "try-error")) return(out)
  co <- try(summary(fit)$coefficients$cond, silent = TRUE)
  if (inherits(co, "try-error") || !"allelepaternal" %in% rownames(co)) return(out)
  out[c("beta", "se", "z", "p")] <- as.list(co["allelepaternal", 1:4])
  out$direction <- if (is.na(out$beta)) NA_character_ else if (out$beta > 0) "paternal" else "maternal"
  out
}

.rmvnorm2 <- function(n, mu, Sigma) {
  L <- try(chol(Sigma), silent = TRUE); eps <- 1e-8
  while (inherits(L, "try-error") && eps <= 1e-2) { L <- try(chol(Sigma + diag(eps, nrow(Sigma))), silent = TRUE); eps <- eps * 10 }
  if (inherits(L, "try-error")) L <- diag(length(mu))
  matrix(stats::rnorm(n * length(mu)), n) %*% L + matrix(mu, n, length(mu), byrow = TRUE)
}

.dose_cols <- function(prefix, dose_range) {
  b <- paste0(prefix, ifelse(dose_range < 0, paste0("m", abs(dose_range)), dose_range))
  c(b, paste0(b, "_ci_lo"), paste0(b, "_ci_hi"))
}

#' Interaction model for one feature.
#' sub has column .cond; cinfo = list(type, step, A = subset mean, B = subset SD)
fit_poe_interaction_one <- function(sub, cinfo, opt = POE_FIT_DEFAULTS) {
  st <- .feature_stats(sub)
  num <- c("beta_int", "se_int", "z_int", "p_int",
           "beta_POE_cond0", "se_POE_cond0", "z_POE_cond0", "p_POE_cond0",
           "beta_POE_cond1", "se_POE_cond1", "z_POE_cond1", "p_POE_cond1",
           "beta_POE_at_ref", "se_POE_at_ref",
           "se_beta_POE_after", "beta_POE_after_ci_lo", "beta_POE_after_ci_hi",
           "abs_beta_after", "se_abs_beta_after", "abs_beta_after_ci_lo", "abs_beta_after_ci_hi",
           "se_delta_abs_beta",
           "beta_POE_subset", "se_beta_POE_subset", "beta_POE_subset_ci_lo", "beta_POE_subset_ci_hi",
           "beta_POE_after_subset", "se_beta_POE_after_subset",
           "beta_POE_after_subset_ci_lo", "beta_POE_after_subset_ci_hi",
           .dose_cols("beta_POE_at_dose_", opt$dose_range),
           .dose_cols("beta_POE_subset_at_dose_", opt$dose_range))
  out <- c(list(feature_id = sub$feature_id[1], cond_type = cinfo$type,
                n_samples = st$n_samples, n_snps = st$n_snps, total_counts = st$total_counts),
           setNames(as.list(rep(NA_real_, length(num))), num),
           list(subset_mean_cond = cinfo$A, subset_sd_cond = cinfo$B,
                vcov_source = NA_character_, boot_status = NA_character_))
  done <- function(status) { out$boot_status <- status; as.data.frame(out, stringsAsFactors = FALSE) }
  if (st$n_samples < opt$min_samples || st$total_counts < opt$min_total_counts || st$n_snps < opt$min_snps)
    return(done("filtered_pre_fit"))

  fit <- try(glmmTMB::glmmTMB(count ~ allele * .cond + ancestry + sex + (1 | SNPpos) + (1 | sample_id),
                              family = glmmTMB::nbinom2(), offset = offset,
                              data = .prep_factors(sub), control = .poe_control(opt)), silent = TRUE)
  if (inherits(fit, "try-error")) return(done("fit_error"))
  co <- try(summary(fit)$coefficients$cond, silent = TRUE)
  if (inherits(co, "try-error") || is.null(co)) return(done("summary_error"))
  ka <- "allelepaternal"; ki <- "allelepaternal:.cond"
  if (!ka %in% rownames(co)) return(done("no_allele_coef"))
  if (!ki %in% rownames(co)) return(done("no_int_coef"))
  out[c("beta_int", "se_int", "z_int", "p_int")] <- as.list(co[ki, 1:4])
  V <- try(as.matrix(stats::vcov(fit)$cond), silent = TRUE)
  if (inherits(V, "try-error") || !all(c(ka, ki) %in% rownames(V))) return(done("no_vcov"))
  out$vcov_source <- "vcov(fit)$cond"
  V <- V[c(ka, ki), c(ka, ki)]
  if (!all(is.finite(V))) return(done("vcov_nonfinite"))
  b0 <- co[ka, 1]; bi <- co[ki, 1]; z <- stats::qnorm(0.975)
  lin <- function(x) {                       # beta_POE at condition value x, delta-method SE
    v <- V[1, 1] + x^2 * V[2, 2] + 2 * x * V[1, 2]
    c(est = b0 + bi * x, se = if (is.finite(v) && v >= 0) sqrt(v) else NA_real_)
  }

  if (cinfo$type == "binary") {
    out[c("beta_POE_cond0", "se_POE_cond0", "z_POE_cond0", "p_POE_cond0")] <- as.list(co[ka, 1:4])
    c1 <- lin(1)
    out$beta_POE_cond1 <- c1[["est"]]; out$se_POE_cond1 <- c1[["se"]]
    out$z_POE_cond1 <- c1[["est"]] / c1[["se"]]
    out$p_POE_cond1 <- 2 * stats::pnorm(-abs(out$z_POE_cond1))
    ref <- c(est = b0, se = co[ka, 2])
  } else {
    out$beta_POE_at_ref <- b0; out$se_POE_at_ref <- co[ka, 2]
  }

  # POE one step above the reference (continuous: +1 SD of X; binary: group 1),
  # with parametric draws for |beta| summaries
  step <- cinfo$step
  draws <- try(.rmvnorm2(max(10L, as.integer(opt$n_boot)), c(b0, bi), V), silent = TRUE)
  if (inherits(draws, "try-error")) return(done("boot_failed"))
  d0 <- draws[, 1]; d1 <- draws[, 1] + draws[, 2] * step
  after <- b0 + bi * step
  out$se_delta_abs_beta <- stats::sd(abs(d1) - abs(d0))
  out$abs_beta_after <- abs(after)
  out$se_abs_beta_after <- stats::sd(abs(d1))
  out$abs_beta_after_ci_lo <- max(0, abs(after) - z * out$se_abs_beta_after)
  out$abs_beta_after_ci_hi <- abs(after) + z * out$se_abs_beta_after
  out$se_beta_POE_after <- stats::sd(d1)
  out$beta_POE_after_ci_lo <- after - z * out$se_beta_POE_after
  out$beta_POE_after_ci_hi <- after + z * out$se_beta_POE_after

  if (cinfo$type != "binary") {
    set_ci <- function(prefix, x) {
      e <- lin(x)
      out[[prefix]] <<- e[["est"]]
      out[[paste0(prefix, "_ci_lo")]] <<- e[["est"]] - z * e[["se"]]
      out[[paste0(prefix, "_ci_hi")]] <<- e[["est"]] + z * e[["se"]]
      e
    }
    e <- set_ci("beta_POE_subset", cinfo$A);            out$se_beta_POE_subset <- e[["se"]]
    e <- set_ci("beta_POE_after_subset", cinfo$A + cinfo$B); out$se_beta_POE_after_subset <- e[["se"]]
    for (d in opt$dose_range) {
      suf <- if (d < 0) paste0("m", abs(d)) else d
      set_ci(paste0("beta_POE_at_dose_", suf), d)                          # X = d
      set_ci(paste0("beta_POE_subset_at_dose_", suf), cinfo$A + d * cinfo$B)  # d SDs from the subset mean
    }
  }
  done("ok")
}

# ---- combining ---------------------------------------------------------------

#' Baseline: BH across features + derived columns
poe_finish_baseline <- function(res) {
  res <- as.data.frame(res)
  res$q <- NA_real_
  ok <- is.finite(res$p)
  res$q[ok] <- stats::p.adjust(res$p[ok], method = "BH")
  res$log2_ratio_p_over_m <- res$beta / log(2)
  res$maternal_prop_imp <- ifelse(is.finite(res$beta), 1 / (1 + exp(res$beta)), NA_real_)
  res <- res[order(res$q, res$p, na.last = TRUE), ]
  rownames(res) <- NULL
  res
}

#' Interaction: BH across features within a condition + derived columns
poe_finish_interaction <- function(res, cinfo, level, condition, opt = POE_FIT_DEFAULTS) {
  res <- as.data.frame(res)
  fin <- function(x) ifelse(is.finite(x), x, NA_real_)
  res$q_int <- NA_real_
  ok <- is.finite(res$p_int)
  res$q_int[ok] <- stats::p.adjust(res$p_int[ok], method = "BH")
  res$log2_POE_slope_per_SD <- fin(res$beta_int) / log(2)
  res$log2_POE_cond0 <- fin(res$beta_POE_cond0) / log(2)
  res$log2_POE_cond1 <- fin(res$beta_POE_cond1) / log(2)
  res$paternal_prop_cond0 <- stats::plogis(fin(res$beta_POE_cond0))
  res$paternal_prop_cond1 <- stats::plogis(fin(res$beta_POE_cond1))
  if (cinfo$type == "continuous") {
    res$beta_POE_before <- res$beta_POE_at_ref
    res$beta_POE_after  <- res$beta_POE_at_ref + res$beta_int * cinfo$step
  } else {
    res$beta_POE_before <- res$beta_POE_cond0
    res$beta_POE_after  <- res$beta_POE_cond1
  }
  b0 <- fin(res$beta_POE_before); b1 <- fin(res$beta_POE_after)
  res$x_enhance <- abs(b1) - abs(b0)
  res$delta_abs_beta <- abs(res$x_enhance)
  res$abs_beta_int <- abs(res$beta_int)
  res$paternal_prop_before <- stats::plogis(b0)
  res$paternal_prop_after  <- stats::plogis(b1)
  res$log2_POE_before <- b0 / log(2)
  res$log2_POE_after  <- b1 / log(2)
  b0s <- fin(res$beta_POE_subset); b1s <- fin(res$beta_POE_after_subset)
  res$paternal_prop_subset <- stats::plogis(b0s)
  res$paternal_prop_after_subset <- stats::plogis(b1s)
  res$log2_POE_subset <- b0s / log(2)
  res$log2_POE_after_subset <- b1s / log(2)
  res$x_enhance_subset <- abs(b1s) - abs(b0s)
  res$pass_pre_counts <- with(res, n_samples >= opt$min_samples & n_snps >= opt$min_snps &
                                     total_counts >= opt$min_total_counts)
  res$pass_diff_strict <- is.finite(res$q_int) & res$q_int <= 0.05
  res$level <- level; res$condition <- condition; res$cond_type <- cinfo$type
  res <- res[order(res$q_int, res$p_int, na.last = TRUE), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "condition") <- condition; attr(res, "cond_type") <- cinfo$type
  attr(res, "delta_step_for_after") <- cinfo$step; attr(res, "level") <- level
  attr(res, "cond_levels") <- cinfo$levels
  attr(res, "min_snps") <- opt$min_snps; attr(res, "min_samples") <- opt$min_samples
  attr(res, "min_total_counts") <- opt$min_total_counts
  res
}

# ---- condition coding ------------------------------------------------------------

#' Code one metadata column for the interaction model
#' continuous: used as is, step = SD; binary (non-numeric, or numeric with two
#' values): 0 = first sorted level (reference), 1 = second, step = 1
poe_code_condition <- function(x) {
  u <- sort(unique(stats::na.omit(as.character(x))))
  if (!is.numeric(x) || length(u) == 2L) {
    if (length(u) != 2L) stop("non-numeric condition must have exactly 2 levels; found ", length(u))
    if (is.factor(x)) u <- intersect(levels(x), u)
    v <- ifelse(is.na(x), NA_real_, as.numeric(as.character(x) == u[2]))
    return(list(x = v, type = "binary", step = 1, levels = u))
  }
  list(x = as.numeric(x), type = "continuous", step = stats::sd(x, na.rm = TRUE), levels = NULL)
}

# ---- batches -----------------------------------------------------------------
# Array element i of the fit steps maps to gene batch i (i <= N_BATCHES_GENE)
# or transcript batch i - N_BATCHES_GENE.

poe_n_batches <- function() c(gene = as.integer(P("N_BATCHES_GENE")),
                              transcript = as.integer(P("N_BATCHES_TX")))

poe_batch_of <- function(idx) {
  nb <- poe_n_batches()
  if (idx < 1 || idx > sum(nb)) stop("batch index ", idx, " outside 1..", sum(nb))
  if (idx <= nb[["gene"]]) list(level = "gene", k = idx) else list(level = "transcript", k = idx - nb[["gene"]])
}

poe_batch_input  <- function(level, k) file.path(P("BATCH_DIR"), sprintf("%s_batch_%03d.rds", level, k))
poe_batch_result <- function(model, level, k)
  file.path(P("BATCH_DIR"), "results", sprintf("%s_%s_%03d.rds", model, level, k))

#' Assign features to n batches with similar numbers of rows (largest first, snake order)
poe_assign_batches <- function(rows_per_feature, n) {
  o <- order(rows_per_feature, decreasing = TRUE)
  m <- length(o)
  snake <- rep(c(seq_len(n), rev(seq_len(n))), length.out = m)
  setNames(snake, names(rows_per_feature)[o])
}

#' Deterministic per-feature seed (independent of core count and scheduling)
poe_seed <- function(feature_id) {
  b <- utf8ToInt(feature_id)
  sum(b * seq_along(b)) %% 2147483L + 1L
}
