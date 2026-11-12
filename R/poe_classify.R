# =============================================================================
# R/poe_classify.R -- significance tiers for the NB-GLMM results
# classify_poe_nbmm():      baseline POE (BH q and 95% CI of beta vs |beta| threshold)
# classify_poe_nbmm_diff(): allele x condition interaction (q_int and |beta_int|)
# =============================================================================

classify_poe_nbmm <- function(
    res,
    val_thresh       = 0.10,    # BH q threshold
    min_abs_beta     = log(1.5),
    ci               = 0.95,
    min_snps         = 0L,
    min_total_counts = 0L,
    min_sample_prop  = 0,
    require_CI_low   = FALSE
) {
  # Working columns
  stopifnot(is.data.frame(res), all(c("beta","se","q") %in% names(res)))
  beta_col <- "beta"
  se_col <- "se"
  val_col <- "q"
  p_col <- "p"
  z_col <- "z"

  
  # Optional columns; if missing, create permissive placeholders
  if (!"n_snps"       %in% names(res)) res$n_snps       <- Inf
  if (!"n_samples"    %in% names(res)) res$n_samples    <- Inf
  if (!"total_counts" %in% names(res)) res$total_counts <- Inf
  if (!p_col          %in% names(res)) res[[p_col]]     <- NA_real_
  if (!z_col          %in% names(res)) res[[z_col]]     <- NA_real_
  
  # Extract the working columns
  beta_work <- res[[beta_col]]
  se_work <- res[[se_col]]
  val_work <- res[[val_col]]
  p_work <- res[[p_col]]
  z_work <- res[[z_col]]
  
  # helpers
  inv_logit <- function(x) 1/(1+exp(-x))
  
  # Total matched samples (if recorded in size_factors attribute)
  sf    <- attr(res, "size_factors")
  n_all <- if (!is.null(sf)) length(sf) else suppressWarnings(max(res$n_samples, na.rm = TRUE))
  if (!is.finite(n_all) || n_all <= 0) n_all <- NA_integer_
  
  # z* for CI
  zcrit <- qnorm(1 - (1 - ci)/2)
  
  # Prefer SE back-calculated from z when available (keeps CI consistent with the test)
  se_from_z <- ifelse(is.finite(z_work) & z_work != 0, abs(beta_work) / abs(z_work), NA_real_)
  se_use    <- ifelse(is.finite(se_from_z), se_from_z, se_work)
  se_use[!is.finite(se_use) | se_use <= 0] <- NA_real_
  
  # CI endpoints
  beta_lo <- beta_work - zcrit * se_use
  beta_hi <- beta_work + zcrit * se_use
  
  # Store CI endpoints
  ci_lo_col <- "beta_lo"
  ci_hi_col <- "beta_hi"
  res[[ci_lo_col]] <- beta_lo
  res[[ci_hi_col]] <- beta_hi
  
  # Threshold checks using the working beta values
  abs_beta <- abs(beta_work)
  ci_outside_threshold <- (beta_lo > min_abs_beta) | (beta_hi < -min_abs_beta)
  
  # Store derived values
  abs_beta_col <- "abs_beta"
  res[[abs_beta_col]] <- abs_beta
  res$ci_outside_threshold <- ci_outside_threshold
  
  # Convenience: implied maternal fractions from beta values (for display)
  maternal_col <- "maternal_frac"
  maternal_lo_col <- "maternal_frac_lo"
  maternal_hi_col <- "maternal_frac_hi"
  
  res[[maternal_col]] <- 1/(1+exp(beta_work))      # plogis(-beta)
  res[[maternal_lo_col]] <- 1/(1+exp(beta_lo))
  res[[maternal_hi_col]] <- 1/(1+exp(beta_hi))
  
  # Core numeric sanity
  ok_core <- is.finite(beta_work) & is.finite(se_use) & is.finite(val_work)
  
  # Information thresholds
  min_ns <- if (!is.na(n_all)) ceiling(min_sample_prop * n_all) else 0L
  info_ok <- (res$n_snps >= min_snps) &
    (res$total_counts >= min_total_counts) &
    (res$n_samples >= min_ns)
  
  # Final significance logic
  meets_val <- ok_core & (val_work <= val_thresh)  # Use val_work and val_thresh
  if (require_CI_low) {
    meets_thr <- ok_core & ci_outside_threshold
  } else {
    meets_thr <- ok_core & (abs_beta >= min_abs_beta)
  }
  sig_idx <- info_ok & meets_val & meets_thr
  
  # Order helper - use the working values for ordering
  ord_fun <- function(d) {
    if (!nrow(d)) return(d)
    d[order(d[[val_col]], -d[[abs_beta_col]], d[[p_col]], na.last = TRUE), , drop = FALSE]
  }
  
  # Diagnostics columns
  res$meets_val            <- meets_val
  res$meets_threshold      <- meets_thr
  res$used_ci_requirement  <- require_CI_low
  
  out <- list(
    significant     = ord_fun(res[sig_idx, , drop = FALSE]),
    not_significant = ord_fun(res[!sig_idx, , drop = FALSE])
  )
  attr(out$significant, "size_factors")     <- sf
  attr(out$not_significant, "size_factors") <- sf
  out
}



classify_poe_nbmm_diff <- function(
    res,                       # output of fit_poe_nbmm_by_condition()
    # FDR thresholds (interaction only)
    fdr_strict = 0.05,
    fdr_loose  = 0.10,
    # Effect-size threshold: log fold-change on beta_int (absolute)
    # e.g., log(2) = 2-fold change
    min_abs_lfc_int = log(2)
) {
  stopifnot(is.data.frame(res),
            all(c("beta_int","q_int") %in% names(res)))
  
  # ---- convenience columns (not used for thresholds) ----
  if (!("x_enhance" %in% names(res)) || !("delta_abs_beta" %in% names(res))) {
    x_enh <- rep(NA_real_, nrow(res))
    dabs  <- rep(NA_real_, nrow(res))
    if (all(c("beta_POE_before","beta_POE_after") %in% names(res))) {
      b0 <- res$beta_POE_before; b1 <- res$beta_POE_after
      x_enh <- abs(b1) - abs(b0); dabs <- abs(x_enh)
    } else {
      is_cont <- if ("cond_type" %in% names(res)) res$cond_type == "continuous" else FALSE
      is_bin  <- if ("cond_type" %in% names(res)) res$cond_type == "binary"     else FALSE
      if (any(is_cont, na.rm = TRUE) && all(c("beta_POE_at_ref","beta_int") %in% names(res))) {
        b0 <- res$beta_POE_at_ref; b1 <- b0 + res$beta_int * (attr(res, "delta_step_for_after") %||% 1)
        x_enh[is_cont] <- abs(b1[is_cont]) - abs(b0[is_cont])
        dabs[is_cont]  <- abs(x_enh[is_cont])
      }
      if (any(is_bin, na.rm = TRUE) && all(c("beta_POE_cond0","beta_POE_cond1") %in% names(res))) {
        b0 <- res$beta_POE_cond0; b1 <- res$beta_POE_cond1
        x_enh[is_bin] <- abs(b1[is_bin]) - abs(b0[is_bin])
        dabs[is_bin]  <- abs(x_enh[is_bin])
      }
    }
    if (!"x_enhance" %in% names(res))      res$x_enhance      <- x_enh
    if (!"delta_abs_beta" %in% names(res)) res$delta_abs_beta <- dabs
  }
  if (!("se_delta_abs_beta" %in% names(res))) res$se_delta_abs_beta <- NA_real_
  
  # ---- absolute log fold change on interaction ----
  res$abs_beta_int <- abs(res$beta_int)              # treat beta_int as log fold-change
  res$lfc_int_ok   <- is.finite(res$abs_beta_int) & (res$abs_beta_int >= min_abs_lfc_int)
  
  # ---- FDR gates (interaction only) ----
  if (!"q_int" %in% names(res)) res$q_int <- NA_real_
  fdr_ok_diff_strict <- is.finite(res$q_int) & (res$q_int <= fdr_strict)
  fdr_ok_diff_loose  <- is.finite(res$q_int) & (res$q_int <= fdr_loose)
  
  # ---- Final passes (ONLY q + log fold-change on beta_int) ----
  res$pass_diff_strict <- fdr_ok_diff_strict & res$lfc_int_ok
  res$pass_diff_loose  <- fdr_ok_diff_loose  & res$lfc_int_ok
  
  # Aliases
  res$pass_strict <- res$pass_diff_strict
  res$pass_loose  <- res$pass_diff_loose
  
  # ---- Order & return ----
  ord_q <- ifelse(is.finite(res$q_int), res$q_int,
                  ifelse(is.finite(res$p_int), res$p_int, NA_real_))
  res <- res[order(ord_q, -res$abs_beta_int, na.last = TRUE), ]
  rownames(res) <- NULL
  res
}

