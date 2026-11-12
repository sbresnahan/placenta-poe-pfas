# =============================================================================
# R/poe_within_sample.R -- within-sample parent-of-origin tests
# Per-SNP exact binomial tests (maternal vs paternal counts), ACAT combination
# per feature within each sample, BH across features within the sample.
# =============================================================================

poe_tx_tests_hierarchical <- function(
    x,
    q_threshold = 0.10,
    p_adjust_method = "BH",
    acat_weights = c("equal","counts"),
    rounding = c("round","floor","ceiling","none"),
    n_workers = max(1, parallel::detectCores() - 1),
    BPPARAM = NULL,
    verbose = TRUE
) {
  acat_weights <- match.arg(acat_weights)
  rounding     <- match.arg(rounding)
  
  if (!requireNamespace("BiocParallel", quietly = TRUE))
    stop("Please install BiocParallel: BiocManager::install('BiocParallel')")
  
  # backend
  if (is.null(BPPARAM)) {
    if (.Platform$OS.type == "windows") {
      BPPARAM <- BiocParallel::SnowParam(workers = n_workers, type = "PSOCK", progressbar = verbose)
    } else {
      BPPARAM <- BiocParallel::MulticoreParam(workers = n_workers, progressbar = verbose)
    }
  }
  
  stopifnot(!is.null(rownames(x)))
  rn <- rownames(x)
  
  # 1) STRICT: keep only rows that end with :<digits>
  has_numeric_tail <- grepl(":\\d+$", rn)
  if (!all(has_numeric_tail)) {
    if (verbose) message("Filtered out ", sum(!has_numeric_tail), " rows without trailing :<digits>.")
    x  <- x[has_numeric_tail, , drop = FALSE]
    rn <- rownames(x)
  }
  
  # LAST-colon split
  feature_id <- sub("^(.*):[^:]+$", "\\1", rn)       # before last ':'
  snppos_str <- sub("^.*:([^:]+)$", "\\1", rn)       # after  last ':'
  suppressWarnings(SNPpos <- as.integer(snppos_str)) # now guaranteed numeric by filter above
  
  nm <- names(x)
  mat_idx <- grep("_maternal$", nm)
  pat_idx <- grep("_paternal$", nm)
  if (length(mat_idx) == 0L || length(pat_idx) == 0L)
    stop("Expected columns named like '<sample>_maternal' and '<sample>_paternal'.")
  
  sample_from_name <- function(v) sub("_(maternal|paternal)$", "", v)
  mat_samples <- sample_from_name(nm[mat_idx])
  pat_samples <- sample_from_name(nm[pat_idx])
  common_samples <- intersect(mat_samples, pat_samples)
  if (length(common_samples) == 0L)
    stop("No matched maternal/paternal sample column pairs found.")
  
  mat_use <- mat_idx[match(common_samples, mat_samples)]
  pat_use <- pat_idx[match(common_samples, pat_samples)]
  
  X <- as.matrix(x[, c(mat_use, pat_use), drop = FALSE])
  storage.mode(X) <- "double"
  M <- X[, seq_along(common_samples), drop = FALSE]
  P <- X[, length(common_samples) + seq_along(common_samples), drop = FALSE]
  colnames(M) <- colnames(P) <- common_samples
  
  Tot   <- M + P
  Ratio <- ifelse(Tot > 0, M / Tot, NA_real_)
  
  round_fun <- switch(rounding,
                      round   = base::round,
                      floor   = base::floor,
                      ceiling = base::ceiling,
                      none    = function(z) z)
  
  pvals_for_col <- function(m_col, t_col) {
    m_int <- round_fun(m_col); t_int <- round_fun(t_col)
    if (identical(round_fun, function(z) z)) {
      if (any(abs(m_int - m_col) > 1e-8 | abs(t_int - t_col) > 1e-8, na.rm = TRUE)) {
        m_int <- base::round(m_col); t_int <- base::round(t_col)
      }
    }
    vapply(seq_along(m_int), function(i) {
      mi <- m_int[i]; ti <- t_int[i]
      if (is.na(mi) || is.na(ti) || ti <= 0L || mi < 0L || mi > ti) return(NA_real_)
      stats::binom.test(mi, ti, p = 0.5, alternative = "two.sided")$p.value
    }, numeric(1))
  }
  
  acat <- function(pvals, weights = NULL) {
    p <- pvals[!is.na(pvals)]
    if (length(p) == 0L) return(NA_real_)
    eps_lo <- 1e-300; eps_hi <- 1 - 1e-16
    p[p <= 0] <- eps_lo
    p[p >= 1] <- eps_hi
    if (is.null(weights)) {
      w <- rep(1, length(p))
    } else {
      w <- weights[!is.na(pvals)]
      if (length(w) != length(p)) w <- rep(1, length(p))
    }
    w <- w / sum(w)
    t_sum <- sum(w * tan((0.5 - p) * pi))
    p_acat <- 0.5 - atan(t_sum) / pi
    max(min(p_acat, 1), 0)
  }
  
  if (verbose) message("Hierarchical testing with ACAT: per-SNP exact tests, ACAT within features, FDR across features…")
  
  per_sample_results <- BiocParallel::bplapply(seq_along(common_samples), function(j) {
    s <- common_samples[j]
    
    # 2) PER-SAMPLE GUARD: only rows that are usable in this sample
    keep_rows <- is.finite(SNPpos) & nzchar(feature_id) & (Tot[, j] > 0)
    if (!any(keep_rows)) {
      empty_tests <- data.frame(
        sample_id=character(0), feature_id=character(0), SNPpos=integer(0),
        total_counts=numeric(0), ratio=numeric(0), p_value=numeric(0),
        check.names = FALSE
      )
      empty_summary <- data.frame(
        sample_id=character(0), feature_id=character(0),
        ratio_summary=numeric(0), combined_p=numeric(0),
        combined_q=numeric(0), POE=logical(0), check.names = FALSE
      )
      return(list(tests = empty_tests, summary = empty_summary))
    }
    
    tj_feature <- feature_id[keep_rows]
    tj_SNPpos  <- SNPpos[keep_rows]
    m_j        <- M[keep_rows, j]
    t_j        <- Tot[keep_rows, j]
    r_j        <- Ratio[keep_rows, j]
    
    pvals <- pvals_for_col(m_j, t_j)
    
    tests_j <- data.frame(
      sample_id     = s,
      feature_id    = tj_feature,
      SNPpos        = tj_SNPpos,
      total_counts  = as.numeric(t_j),
      ratio         = as.numeric(r_j),
      p_value       = as.numeric(pvals),
      row.names = NULL, check.names = FALSE
    )
    
    # rollup counts for ratio_summary
    tests_j$maternal_counts <- tests_j$ratio * tests_j$total_counts
    
    agg_counts <- stats::aggregate(
      cbind(total_counts, maternal_counts) ~ sample_id + feature_id,
      data = tests_j, FUN = sum, na.rm = TRUE
    )
    ratio_summary <- with(agg_counts, ifelse(total_counts > 0, maternal_counts / total_counts, NA_real_))
    
    split_inds <- split(seq_len(nrow(tests_j)), tests_j$feature_id)
    weights_vec <- if (acat_weights == "counts") tests_j$total_counts else NULL
    
    combined_p_map <- vapply(split_inds, function(inds) {
      pv <- tests_j$p_value[inds]
      ww <- if (is.null(weights_vec)) NULL else weights_vec[inds]
      acat(pv, weights = ww)
    }, numeric(1))
    
    tr_order <- agg_counts$feature_id
    combined_p_vec <- unname(as.numeric(combined_p_map[tr_order]))
    combined_q_vec <- rep(NA_real_, length(combined_p_vec))
    ok <- is.finite(combined_p_vec)
    if (any(ok)) combined_q_vec[ok] <- p.adjust(combined_p_vec[ok], method = p_adjust_method)
    
    summary_j <- data.frame(
      sample_id     = agg_counts$sample_id,
      feature_id    = tr_order,
      ratio_summary = ratio_summary,
      combined_p    = combined_p_vec,
      combined_q    = combined_q_vec,
      POE           = combined_q_vec < q_threshold,
      row.names = NULL, check.names = FALSE
    )
    summary_j$POE[is.na(summary_j$combined_q)] <- FALSE
    
    list(
      tests   = tests_j[, c("sample_id","feature_id","SNPpos","total_counts","ratio","p_value")],
      summary = summary_j
    )
  }, BPPARAM = BPPARAM)
  
  # combine outputs
  tests_list   <- lapply(per_sample_results, `[[`, "tests")
  summary_list <- lapply(per_sample_results, `[[`, "summary")
  tests   <- do.call(rbind, tests_list)
  summary <- do.call(rbind, summary_list)
  
  # stable ordering
  if (nrow(tests))   tests   <- tests[order(tests$sample_id, tests$feature_id, tests$SNPpos), ]
  if (nrow(summary)) summary <- summary[order(summary$sample_id, summary$feature_id), ]
  
  list(tests = tests, summary = summary)
}

