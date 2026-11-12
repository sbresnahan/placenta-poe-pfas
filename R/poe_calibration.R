# =============================================================================
# R/poe_calibration.R -- model diagnostics
# Observed vs predicted paternal proportion per feature (calibration) and QQ
# plots of baseline and interaction p-values.
# =============================================================================

plot_calibration_observed_vs_predicted <- function(
    res,
    poMat,
    ensembl,                 # data.frame(gene_id, symbol)
    tx2g = NULL,             # transcript->gene map if level="transcript"
    features = NULL,         # optional vector of gene IDs, tx IDs, or symbols
    min_total = 10,          # per-sample minimum total counts (M+P) to keep sample for pooling
    level = c("gene","transcript"),
    smooth_method = "loess" # "loess", "gam", or "lm"
) {
  level <- match.arg(level, c("gene","transcript"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2')")
  
  # --- helpers ---
  inv_logit <- function(x) 1/(1+exp(-x))
  feature_from_rownm <- function(rn) sub(":([^:]*)$", "", rn)  # drop after LAST colon
  has_cols <- function(df, cols) all(cols %in% names(df))
  
  convert_to_symbols <- function(ids, ensembl, tx2g = NULL, level = c("gene","transcript")) {
    level <- match.arg(level)
    if (is.null(ensembl) || !has_cols(ensembl, c("gene_id","symbol"))) return(rep(NA_character_, length(ids)))
    if (level == "gene") {
      gid <- sub("\\.\\d+$","", ids)
      sym <- ensembl$symbol[match(gid, sub("\\.\\d+$","", ensembl$gene_id))]
      return(sym)
    } else {
      if (is.null(tx2g) || !has_cols(tx2g, c("transcript_id","gene_id"))) return(rep(NA_character_, length(ids)))
      gid <- tx2g$gene_id[match(ids, tx2g$transcript_id)]
      gid <- sub("\\.\\d+$","", gid)
      sym <- ensembl$symbol[match(gid, sub("\\.\\d+$","", ensembl$gene_id))]
      return(sym)
    }
  }
  
  # --- baseline POE beta column ---
  beta_col <- if ("beta" %in% names(res)) {
    "beta"
  } else if ("beta_POE_at_ref" %in% names(res)) {
    "beta_POE_at_ref"
  } else if ("beta_POE" %in% names(res)) {
    "beta_POE"
  } else {
    stop("no beta column found. Expected one of: beta, beta_POE_at_ref, beta_POE")
  }

  
  if (!has_cols(res, c("feature_id", beta_col))) {
    stop("`res` must include columns: feature_id and ", beta_col)
  }
  
  # --- optional feature subsetting (IDs and/or symbols) ---
  feature_filter_ids <- NULL
  if (!is.null(features)) {
    symbols_in <- if (!is.null(ensembl)) intersect(features, unique(ensembl$symbol)) else character(0)
    ids_from_symbols <- if (length(symbols_in) && !is.null(ensembl)) {
      out <- ensembl$gene_id[match(symbols_in, ensembl$symbol)]
      out[!is.na(out)]
    } else character(0)
    direct_ids <- setdiff(features, symbols_in)
    
    if (level == "transcript" && length(ids_from_symbols) && !is.null(tx2g)) {
      tx_keep <- tx2g$transcript_id[tx2g$gene_id %in% sub("\\.\\d+$","", ids_from_symbols)]
      feature_filter_ids <- unique(c(direct_ids, tx_keep))
    } else {
      feature_filter_ids <- unique(c(direct_ids, ids_from_symbols))
    }
    if (!length(feature_filter_ids)) stop("No features matched `features` after mapping.")
  }
  
  # --- restrict result rows ---
  res_use <- if (is.null(feature_filter_ids)) res else res[res$feature_id %in% feature_filter_ids, , drop = FALSE]
  if (!nrow(res_use)) stop("No rows left in `res` after subsetting.")
  
  # --- parse poMat into maternal/paternal count matrices ---
  rn <- rownames(poMat)
  if (is.null(rn)) stop("`poMat` rownames must be 'feature_id:SNPpos'.")
  feat_ids_all <- feature_from_rownm(rn)
  
  col_m <- grep("_maternal$", colnames(poMat), value = TRUE)
  if (!length(col_m)) stop("No '*_maternal' columns in `poMat`.")
  col_p <- sub("_maternal$", "_paternal", col_m)
  if (!all(col_p %in% colnames(poMat))) stop("Missing matching '*_paternal' columns in `poMat`.")
  
  # optional row filter to the requested features
  if (!is.null(feature_filter_ids)) {
    keep_row <- feat_ids_all %in% feature_filter_ids
    if (!any(keep_row)) stop("No poMat rows match requested `features`.")
    poMat <- poMat[keep_row, , drop = FALSE]
    feat_ids_all <- feat_ids_all[keep_row]
  }
  
  # matrices, standardized sample IDs as column names
  M_mat <- as.matrix(poMat[, col_m, drop = FALSE])
  P_mat <- as.matrix(poMat[, col_p, drop = FALSE])
  colnames(M_mat) <- sub("_maternal$", "", colnames(M_mat))
  colnames(P_mat) <- sub("_paternal$",  "", colnames(P_mat))
  
  # --- aggregate SNP rows -> feature-level totals per sample (pooled = sum) ---
  M_byfeat <- rowsum(M_mat, group = feat_ids_all, reorder = FALSE, na.rm = TRUE)
  P_byfeat <- rowsum(P_mat, group = feat_ids_all, reorder = FALSE, na.rm = TRUE)
  
  Mm <- as.data.frame(as.table(M_byfeat)); names(Mm) <- c("feature_id","sample_id","M")
  Pm <- as.data.frame(as.table(P_byfeat)); names(Pm) <- c("feature_id","sample_id","P")
  
  obs_long <- merge(Mm, Pm, by = c("feature_id","sample_id"), all = TRUE)
  obs_long$M[is.na(obs_long$M)] <- 0
  obs_long$P[is.na(obs_long$P)] <- 0
  obs_long$total <- obs_long$M + obs_long$P
  obs_long <- obs_long[is.finite(obs_long$total) & obs_long$total >= min_total, , drop = FALSE]
  if (!nrow(obs_long)) stop("No observed pairs after filtering (check `min_total`).")
  
  # pooled observed per feature (sum across kept samples) + n_samples contributing
  pooled_obs <- aggregate(cbind(P = P, total = total) ~ feature_id, data = obs_long, FUN = sum)
  pooled_obs$p_obs <- with(pooled_obs, ifelse(total > 0, P/total, NA_real_))
  ns <- aggregate(sample_id ~ feature_id, data = obs_long, FUN = function(x) length(unique(x)))
  names(ns)[2] <- "n_samples"
  
  # predicted per feature from baseline beta
  pred <- res_use[, c("feature_id", beta_col)]
  names(pred)[2] <- "beta"
  pred$p_hat <- inv_logit(pred$beta)
  
  # merge + symbols
  pooled <- merge(pooled_obs, pred, by = "feature_id", all.x = TRUE)
  pooled <- merge(pooled, ns, by = "feature_id", all.x = TRUE)
  pooled$symbol <- convert_to_symbols(pooled$feature_id, ensembl, tx2g = tx2g, level = level)
  
  pooled <- pooled[is.finite(pooled$p_obs) & is.finite(pooled$p_hat), , drop = FALSE]
  if (!nrow(pooled)) stop("No finite observed/predicted pairs to plot.")
  
  # --- R^2 for display (cor between p_hat and p_obs) ---
  R2_display <- if (nrow(pooled) > 1) (cor(pooled$p_hat, pooled$p_obs, use = "complete.obs"))^2 else NA_real_
  
  # --- plot ---
  library(ggplot2)
  calib_label <- ""
  ttl <- sprintf("Calibration: Predicted vs Observed (baseline, %s-level; pooled)%s", level, calib_label)
  
  p <- ggplot(pooled, aes(x = p_hat, y = p_obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.7, color = "black") +
    # geom_point(aes(color = n_samples), alpha = 0.85, size = 2.5) +
    geom_point(color="grey40", alpha = 0.85, size = 2.5) +
    geom_smooth(
      inherit.aes = FALSE,
      aes(x = p_hat, y = p_obs, group = 1),
      method = smooth_method, se = TRUE, linewidth = 1, color = "red"
    ) +
    # scale_color_gradient(name = "Samples contributing", low = "skyblue", high = "darkblue") +
    labs(
      title = ttl,
      subtitle = if (!is.na(R2_display)) sprintf("R² (p_hat vs p_obs) = %.3f", R2_display) else NULL,
      x = "Predicted paternal proportion (inv_logit(beta))",
      y = "Observed paternal proportion (pooled across samples)"
    ) +
    coord_equal(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank())
  
  # also return per-sample observed for diagnostics
  obs_long$p_obs_sample <- with(obs_long, ifelse(total > 0, P/total, NA_real_))
  
  # --- return ---
  out <- list(
    plot = p,
    pooled = pooled[order(pooled$feature_id), ],
    sample_level = obs_long[order(obs_long$feature_id, obs_long$sample_id),
                            c("feature_id","sample_id","M","P","total","p_obs_sample")],
    beta_column_used = beta_col
  )
  
  return(out)
}



# QQ plot function
make_glm_qq <- function(
    res,
    mode = c("p", "z"),        # "p" = p-value QQ on -log10 scale; "z" = z-score QQ vs Normal
    p_col = "p",               # column in `res` with p-values
    z_col = "z",               # column in `res` with z-scores
    title = NULL,
    subtitle = NULL,
    point_alpha = 0.6,
    point_size = 1.2,
    ci = 0.95                  # confidence band level (for p-QQ)
) {
  stopifnot(is.data.frame(res))
  mode <- match.arg(mode)
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install.packages('ggplot2')")
  if (!requireNamespace("stats", quietly = TRUE))
    stop("Base 'stats' not available?")
  
  if (mode == "p") {
    p <- res[[p_col]]
    p <- p[is.finite(p) & p >= 0]
    if (!length(p)) stop("No finite p-values found in column '", p_col, "'.")
    # protect against p=0
    p[p == 0] <- .Machine$double.xmin
    
    n <- length(p)
    p_obs <- sort(p)
    # Expected order-statistics of Uniform(0,1) ≈ (i - 0.5)/n
    i <- seq_len(n)
    p_exp <- (i - 0.5) / n
    
    # 95% CI from Beta(i, n-i+1)
    alpha <- (1 - ci) / 2
    lo <- stats::qbeta(alpha, i, n - i + 1)
    hi <- stats::qbeta(1 - alpha, i, n - i + 1)
    
    df <- data.frame(
      exp = -log10(p_exp),
      obs = -log10(p_obs),
      lo  = -log10(hi),  # note: hi on p-scale is lower on -log10 scale
      hi  = -log10(lo)
    )
    
    library(ggplot2)
    gg <- ggplot(df, aes(x = exp, y = obs)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
      geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = 2) +
      geom_point(alpha = point_alpha, size = point_size) +
      labs(
        title    = title %||% "QQ plot of p-values",
        subtitle = subtitle %||% sprintf("n = %d tests; %d%% CI band", n, round(ci*100)),
        x = expression(Expected~~-log[10](p)),
        y = expression(Observed~~-log[10](p))
      ) +
      theme_minimal(base_size = 12)
    
    return(gg)
  } else {
    z <- res[[z_col]]
    z <- z[is.finite(z)]
    if (!length(z)) stop("No finite z-scores found in column '", z_col, "'.")
    n <- length(z)
    z <- sort(z)
    # expected quantiles from standard normal
    q <- stats::qnorm(stats::ppoints(n))
    
    # genomic inflation (optional diagnostic)
    lambda <- median(z^2, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
    
    df <- data.frame(exp = q, obs = z)
    
    library(ggplot2)
    gg <- ggplot(df, aes(x = exp, y = obs)) +
      geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = 2) +
      geom_point(alpha = point_alpha, size = point_size) +
      labs(
        title    = title %||% "QQ plot of z-scores",
        subtitle = subtitle %||% sprintf("n = %d tests; inflation λ = %.3f", n, lambda),
        x = "Expected Normal quantile",
        y = "Observed Normal quantile"
      ) +
      theme_minimal(base_size = 12)
    
    return(gg)
  }
}





# QQ plotting
make_nb_int_qq <- function(
    res,
    mode   = c("p","z"),                 # "p" for -log10(p) QQ, "z" for Normal QQ
    which  = c("interaction","custom"),
    p_col  = NULL,                       # used if which = "custom" or to override defaults
    z_col  = NULL,
    title = NULL,
    subtitle = NULL,
    point_alpha = 0.6,
    point_size  = 1.2,
    ci = 0.95
) {
  stopifnot(is.data.frame(res))
  mode  <- match.arg(mode)
  which <- match.arg(which)
  
  # tiny helper (avoid external dep)
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # Map to your data structure
  if (which == "interaction") {
    p_col <- p_col %||% "p_int"
    z_col <- z_col %||% "z_int"
  } else { # custom
    if (mode == "p" && is.null(p_col)) stop("For which='custom', provide p_col.")
    if (mode == "z" && is.null(z_col)) stop("For which='custom', provide z_col.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Please install.packages('ggplot2')")
  
  if (mode == "p") {
    if (!p_col %in% names(res)) stop("Column '", p_col, "' not found in res.")
    p <- res[[p_col]]
    p <- p[is.finite(p) & p >= 0]
    if (!length(p)) stop("No finite p-values found in column '", p_col, "'.")
    p[p == 0] <- .Machine$double.xmin
    
    n <- length(p)
    p_obs <- sort(p)
    i <- seq_len(n)
    p_exp <- (i - 0.5) / n
    
    alpha <- (1 - ci) / 2
    lo <- stats::qbeta(alpha, i, n - i + 1)
    hi <- stats::qbeta(1 - alpha, i, n - i + 1)
    
    df <- data.frame(
      exp = -log10(p_exp),
      obs = -log10(p_obs),
      lo  = -log10(hi),
      hi  = -log10(lo)
    )
    
    ggplot2::ggplot(df, ggplot2::aes(x = exp, y = obs)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey85") +
      ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = 2) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::labs(
        title    = title %||% sprintf("QQ plot of p-values (%s)", which),
        subtitle = subtitle %||% sprintf("n = %d tests; %d%% CI band", n, round(ci*100)),
        x = expression(Expected~~-log[10](p)),
        y = expression(Observed~~-log[10](p))
      ) +
      ggplot2::theme_minimal(base_size = 12)
    
  } else { # mode == "z"
    if (!z_col %in% names(res)) stop("Column '", z_col, "' not found in res.")
    z <- res[[z_col]]
    z <- z[is.finite(z)]
    if (!length(z)) stop("No finite z-scores found in column '", z_col, "'.")
    z <- sort(z)
    n <- length(z)
    q <- stats::qnorm(stats::ppoints(n))
    lambda <- stats::median(z^2, na.rm = TRUE) / stats::qchisq(0.5, df = 1)
    
    df <- data.frame(exp = q, obs = z)
    
    ggplot2::ggplot(df, ggplot2::aes(x = exp, y = obs)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linewidth = 0.5, linetype = 2) +
      ggplot2::geom_point(alpha = point_alpha, size = point_size) +
      ggplot2::labs(
        title    = title %||% sprintf("QQ plot of z-scores (%s)", which),
        subtitle = subtitle %||% sprintf("n = %d tests; inflation \u03BB = %.3f", n, lambda),
        x = "Expected Normal quantile",
        y = "Observed Normal quantile"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  }
}

