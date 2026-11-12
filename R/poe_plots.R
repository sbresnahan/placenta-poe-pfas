# =============================================================================
# R/poe_plots.R -- figures for stages 07-08
# Volcano and forest plots (baseline and interaction), counts of strict
# interaction hits per PFAS, and dose-response curves of the allelic ratio.
# =============================================================================

poe_forest_nbmm <- function(
    res,
    ensembl, tx2g = NULL,
    tx = FALSE,
    top_n = Inf,
    beta.thresh = .405,                  # threshold on beta scale
    x_scale = c("beta","log","maternal_ratio","beta.ratio"),
    title = "Isoforms with significant POE"
) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(forcats)
  })
  
  x_scale <- match.arg(x_scale)
  
  # ----- gene symbols -----
  if (tx) {
    res <- res %>%
      left_join(tx2g, by = c("feature_id" = "transcript_id")) %>%
      mutate(gene_id = sub("\\..*$","", gene_id))
  } else {
    res$gene_id <- sub("\\..*$","", res$feature_id)
  }
  
  res <- res %>% left_join(ensembl, by="gene_id")
  
  # ----- order (|β| desc, q asc) -----
  ord <- res %>%
    arrange(desc(abs_beta), q, feature_id) %>%
    head(top_n)
  
  # ----- labeling -----
  ord <- ord %>%
    group_by(symbol) %>%
    mutate(
      dup_count = n(),
      dup_index = row_number(),
      label = case_when(
        is.na(symbol) | !nzchar(symbol) ~ feature_id,
        dup_count == 1 ~ symbol,
        dup_count > 1 ~ paste0(symbol, " (", dup_index, ")")
      )
    ) %>%
    ungroup()
  
  # ----- filter: drop CIs fully within exclusion band -----
  ord <- ord %>%
    filter(
      (beta_lo < -beta.thresh & beta_hi < -beta.thresh) |   # CI entirely below threshold
        (beta_lo >  beta.thresh & beta_hi >  beta.thresh) |   # CI entirely above threshold
        (beta_lo < -beta.thresh & beta_hi >  beta.thresh)     # CI spans across thresholds
    )
  
  plot_df <- ord %>%
    mutate(label = factor(label, levels = rev(unique(label))))
  
  # ----- transformations -----
  if (x_scale == "beta") {
    plot_df <- plot_df %>% mutate(x = beta, xmin = beta_lo, xmax = beta_hi)
    xlab <- expression(beta[POE])
    formatter <- function(x) sprintf("%.2f", x)
    vlines <- c(-beta.thresh, 0, beta.thresh)
    
    scale_x <- scale_x_continuous(name = xlab, labels = formatter)
    
  } else if (x_scale == "log") {
    plot_df <- plot_df %>%
      mutate(x = exp(beta), xmin = exp(beta_lo), xmax = exp(beta_hi))
    xlab <- "Fold change (exp(beta))"
    formatter <- function(x) sprintf("%.2f", x)
    vlines <- exp(c(-beta.thresh, 0, beta.thresh))
    
    scale_x <- scale_x_continuous(name = xlab, labels = formatter)
    
  } else if (x_scale == "maternal_ratio") {
    inv_logit <- function(b) 1 / (1 + exp(b))
    to_ratio <- function(p) paste0(round(p*100), ":", round((1-p)*100))
    
    plot_df <- plot_df %>%
      mutate(x = inv_logit(beta), xmin = inv_logit(beta_lo), xmax = inv_logit(beta_hi))
    xlab <- "Maternal allelic ratio"
    formatter <- to_ratio
    vlines <- inv_logit(c(-beta.thresh, 0, beta.thresh))
    
    scale_x <- scale_x_continuous(name = xlab, labels = formatter, limits = c(0,1))
    
  } else if (x_scale == "beta.ratio") {
    inv_logit <- function(b) 1 / (1 + exp(b))
    
    # Cap CIs at [-10, +10] to avoid infinite ratios
    plot_df <- plot_df %>%
      mutate(
        x    = beta,
        xmin = pmax(beta_lo, -10),
        xmax = pmin(beta_hi,  10)
      )
    xlab <- expression(beta[POE])
    formatter <- function(x) sprintf("%.2f", x)
    vlines <- c(-beta.thresh, 0, beta.thresh)
    
    # Secondary axis with full maternal ratio range
    beta_breaks <- c(-10, -2.5, 0, 2.5, 10)
    ratio_labels <- c("100:0", "92:8", "50:50", "8:92", "0:100")
    
    scale_x <- scale_x_continuous(
      name = xlab,
      labels = formatter,
      limits = c(-10, 10),
      sec.axis = sec_axis(
        ~ ., name = "Maternal allelic ratio",
        breaks = beta_breaks,
        labels = ratio_labels
      )
    )
  }
  
  # ----- plot -----
  ggplot(plot_df, aes(x = x, y = label)) +
    geom_vline(xintercept = vlines,
               linetype = c("dotted","dashed","dotted"),
               linewidth = 0.4, color = "grey40") +
    geom_errorbar(orientation = "y", aes(xmin = xmin, xmax = xmax), width = 0.25) +
    geom_point(size = 2.6, color = "black") +
    scale_x +
    labs(
      title = title,
      subtitle = paste0("Threshold on β = ±", beta.thresh, 
                        " (features with CI entirely inside excluded)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      axis.title.y = element_blank()
    )
}


plot_isoform_poe_volcano <- function(
    res_tx,
    tx2g,
    ensembl,
    thr_ratio,                 # majority ratio threshold, e.g. 0.60
    val_thresh = 0.10,         # BH q threshold
    n_top_axis_breaks = 7
){
  stopifnot(is.data.frame(res_tx))
  if (!(is.numeric(thr_ratio) && length(thr_ratio) == 1 && thr_ratio >= 0.5 && thr_ratio < 1))
    stop("thr_ratio must be in [0.5,1). Example: 0.60 for 60:40.")
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Install scales.")
  
  inv_logit <- function(x) 1/(1+exp(-x))
  logit     <- function(p) log(p/(1-p))
  clamp01   <- function(p, eps=1e-9) pmax(pmin(p,1-eps), eps)
  
  # Determine which columns to use
  # Required columns
  if (!all(c("beta", "q") %in% names(res_tx))) {
    stop("Missing required columns: beta, q")
  }
  beta_col <- "beta"
  val_col <- "q"
  se_col <- "se"
  z_col <- "z"

  
  # Extract working values
  beta_work <- res_tx[[beta_col]]
  val_work <- res_tx[[val_col]]
  
  # --- ensure CI endpoints ---
  ci_lo_col <- paste0(gsub("beta", "beta_lo", beta_col))
  ci_hi_col <- paste0(gsub("beta", "beta_hi", beta_col))
  
  if (!all(c(ci_lo_col, ci_hi_col) %in% names(res_tx))) {
    if (!se_col %in% names(res_tx)) stop("Need CI columns or SE to make CIs.")
    zcrit <- qnorm(0.975)
    se_from_z <- if (z_col %in% names(res_tx)) {
      z_work <- res_tx[[z_col]]
      ifelse(is.finite(z_work) & z_work != 0, abs(beta_work)/abs(z_work), NA_real_)
    } else {
      NA_real_
    }
    se_use <- ifelse(is.finite(se_from_z), se_from_z, res_tx[[se_col]])
    res_tx[[ci_lo_col]] <- beta_work - zcrit*se_use
    res_tx[[ci_hi_col]] <- beta_work + zcrit*se_use
  }
  
  # default direction
  if (!"direction" %in% names(res_tx)) {
    res_tx$direction <- ifelse(beta_work > 0, "paternal",
                               ifelse(beta_work < 0, "maternal", NA))
  }
  if (!"used_ci_requirement" %in% names(res_tx)) res_tx$used_ci_requirement <- FALSE
  
  # ----- Set x-axis range and breaks based on desired ratio coverage -----
  # Force standard ratio breakpoints: 0:100, 25:75, 50:50, 75:25, 100:0
  p_std <- c(0.001, 0.25, 0.5, 0.75, 0.999)  # avoid exact 0/1 for logit
  ratio_lbl <- c("0:100", "25:75", "50:50", "75:25", "100:0")
  x_ticks <- logit(p_std)
  
  # Set strict limits at the 0:100 and 100:0 boundaries
  x_limits <- c(x_ticks[1], x_ticks[5])
  lower_breaks <- x_ticks
  
  # ----- reference points (threshold distance from 0) -----
  beta50 <- 0                              # 50:50 always at origin in beta_POE space
  thr_shift <- abs(logit(thr_ratio))       # threshold distance from 0
  
  # Cap plotted values to stay strictly within the ratio display range
  res_tx$beta_shift    <- pmax(pmin(beta_work, x_limits[2]), x_limits[1])
  res_tx$beta_lo_shift <- pmax(pmin(res_tx[[ci_lo_col]], x_limits[2]), x_limits[1])
  res_tx$beta_hi_shift <- pmax(pmin(res_tx[[ci_hi_col]], x_limits[2]), x_limits[1])
  
  # significance & aesthetics - ONLY create sig_plot if it doesn't exist
  if (!"sig_plot" %in% names(res_tx)) {
    # sig if outside ±thr_shift band AND val <= cutoff
    res_tx$sig_plot <- (val_work <= val_thresh) &
      (beta_work <= -thr_shift | beta_work >= thr_shift)
  }
  
  # Create plotting aesthetics based on sig_plot
  res_tx$col_plot <- dplyr::case_when(
    res_tx$sig_plot & res_tx$direction == "paternal" ~ "blue",
    res_tx$sig_plot & res_tx$direction == "maternal" ~ "red",
    TRUE                                              ~ "grey"
  )
  res_tx$alpha_plot <- ifelse(res_tx$sig_plot, 1, 0.35)
  
  # ----- map tx -> gene/symbol -----
  tx2g_unique <- tx2g[!is.na(tx2g$transcript_id) & !duplicated(tx2g$transcript_id),
                      c("transcript_id","gene_id")]
  res_nb_tx_plot <- res_tx |>
    dplyr::left_join(tx2g_unique, by = c("feature_id"="transcript_id")) |>
    dplyr::mutate(gene_id = sub("\\..*$","", as.character(gene_id))) |>
    dplyr::left_join(ensembl, by = "gene_id")
  res_nb_tx_plot$symbol[is.na(res_nb_tx_plot$symbol)] <- res_nb_tx_plot$feature_id[is.na(res_nb_tx_plot$symbol)]
  
  # ----- subtitle -----
  ratio_label_txt <- sprintf("%d:%d", round(thr_ratio*100), round((1-thr_ratio)*100))
  subtxt <- bquote(
    beta[POE] < 0 ~ ": maternal bias; " ~
      beta[POE] > 0 ~ ": paternal bias | " ~
      "Threshold: |" * beta[POE] * "| > " * .(round(thr_shift, 3)) ~
      " (major allele " * .(ratio_label_txt) * ")"
  )
  xlab <- expression(beta[POE])

  
  # Split data for proper layering
  non_sig_data <- subset(res_nb_tx_plot, !sig_plot)
  sig_data <- subset(res_nb_tx_plot, sig_plot)
  
  g <- ggplot2::ggplot(res_nb_tx_plot, ggplot2::aes(x = beta_shift, y = -log10(.data[[val_col]]))) +
    # Plot non-significant points first (grey, behind)
    ggplot2::geom_errorbar(orientation = "y",
      data = non_sig_data,
      ggplot2::aes(xmin = beta_lo_shift, xmax = beta_hi_shift),
      color = "grey", alpha = 0.35, width = 0.1, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = non_sig_data,
      color = "grey", alpha = 0.35, size = 1.6, show.legend = FALSE
    ) +
    # Plot significant points on top (colored, in front)
    ggplot2::geom_errorbar(orientation = "y",
      data = sig_data,
      ggplot2::aes(xmin = beta_lo_shift, xmax = beta_hi_shift, color = col_plot),
      alpha = 1, width = 0.1, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = sig_data,
      ggplot2::aes(color = col_plot),
      alpha = 1, size = 1.6, show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    # guides: dashed center at 50:50 (x=0), dotted at ± threshold
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, color = "black") +
    ggplot2::geom_vline(xintercept = c(-thr_shift, thr_shift), linetype = "dotted", color = "black") +
    ggplot2::geom_hline(yintercept = -log10(val_thresh), linetype = "dashed", color = "black") +
    ggrepel::geom_label_repel(
      data = sig_data,
      ggplot2::aes(label = symbol, color = col_plot),
      alpha = 1,
      fill = "white", label.size = 0.25, label.r = grid::unit(0.1, "lines"),
      box.padding = 1, point.padding = 0, show.legend = FALSE, max.overlaps = 200
    ) +
    ggplot2::labs(
      title = "Isoform-level POE",
      subtitle = subtxt,
      x = xlab,
      y = bquote(-log[10] * "(q)")
    ) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      breaks = lower_breaks,
      labels = function(x) sprintf("%.1f", x),
      sec.axis = ggplot2::sec_axis(
        trans = ~ .,
        name  = "Allelic ratio (paternal : maternal)",
        breaks = x_ticks,
        labels = ratio_lbl
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray30"),
      axis.text.x.top = ggplot2::element_text(size = 10)
    )
  
  g
}


plot_gene_poe_volcano <- function(
    res_gene,
    ensembl,
    thr_ratio,                 # majority ratio threshold, e.g. 0.60
    val_thresh = 0.10,         # BH q threshold
    n_top_axis_breaks = 7
){
  stopifnot(is.data.frame(res_gene))
  if (!(is.numeric(thr_ratio) && length(thr_ratio) == 1 && thr_ratio >= 0.5 && thr_ratio < 1))
    stop("thr_ratio must be in [0.5,1). Example: 0.60 for 60:40.")
  
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Install scales.")
  
  inv_logit <- function(x) 1/(1+exp(-x))
  logit     <- function(p) log(p/(1-p))
  clamp01   <- function(p, eps=1e-9) pmax(pmin(p,1-eps), eps)
  
  # Determine which columns to use
  # Required columns
  if (!all(c("beta", "q") %in% names(res_gene))) {
    stop("Missing required columns: beta, q")
  }
  beta_col <- "beta"
  val_col <- "q"
  se_col <- "se"
  z_col <- "z"

  
  # Extract working values
  beta_work <- res_gene[[beta_col]]
  val_work <- res_gene[[val_col]]
  
  # --- ensure CI endpoints ---
  ci_lo_col <- paste0(gsub("beta", "beta_lo", beta_col))
  ci_hi_col <- paste0(gsub("beta", "beta_hi", beta_col))
  
  if (!all(c(ci_lo_col, ci_hi_col) %in% names(res_gene))) {
    if (!se_col %in% names(res_gene))
      stop("Need CI columns or SE to make CIs.")
    zcrit <- qnorm(0.975)
    se_from_z <- if (z_col %in% names(res_gene)) {
      z_work <- res_gene[[z_col]]
      ifelse(is.finite(z_work) & z_work != 0, abs(beta_work)/abs(z_work), NA_real_)
    } else {
      NA_real_
    }
    se_use <- ifelse(is.finite(se_from_z), se_from_z, res_gene[[se_col]])
    res_gene[[ci_lo_col]] <- beta_work - zcrit*se_use
    res_gene[[ci_hi_col]] <- beta_work + zcrit*se_use
  }
  
  # default direction
  if (!"direction" %in% names(res_gene)) {
    res_gene$direction <- ifelse(beta_work > 0, "paternal",
                                 ifelse(beta_work < 0, "maternal", NA))
  }
  if (!"used_ci_requirement" %in% names(res_gene)) res_gene$used_ci_requirement <- FALSE
  
  # ----- Set x-axis range and breaks based on desired ratio coverage -----
  # Force standard ratio breakpoints: 0:100, 25:75, 50:50, 75:25, 100:0
  p_std <- c(0.001, 0.25, 0.5, 0.75, 0.999)  # avoid exact 0/1 for logit
  ratio_lbl <- c("0:100", "25:75", "50:50", "75:25", "100:0")
  x_ticks <- logit(p_std)
  
  # Set strict limits at the 0:100 and 100:0 boundaries
  x_limits <- c(x_ticks[1], x_ticks[5])
  lower_breaks <- x_ticks
  
  # ----- reference points (threshold distance from 0) -----
  beta50 <- 0                              # 50:50 always at origin in beta_POE space
  thr_shift <- abs(logit(thr_ratio))       # threshold distance from 0
  
  # Cap plotted values to stay strictly within the ratio display range
  res_gene$beta_shift    <- pmax(pmin(beta_work, x_limits[2]), x_limits[1])
  res_gene$beta_lo_shift <- pmax(pmin(res_gene[[ci_lo_col]], x_limits[2]), x_limits[1])
  res_gene$beta_hi_shift <- pmax(pmin(res_gene[[ci_hi_col]], x_limits[2]), x_limits[1])
  
  # significance & aesthetics - ONLY create sig_plot if it doesn't exist
  if (!"sig_plot" %in% names(res_gene)) {
    # sig if outside ±thr_shift band AND val <= cutoff
    res_gene$sig_plot <- (val_work <= val_thresh) &
      (beta_work <= -thr_shift | beta_work >= thr_shift)
  }
  
  # Create plotting aesthetics based on sig_plot
  res_gene$col_plot <- dplyr::case_when(
    res_gene$sig_plot & res_gene$direction == "paternal" ~ "blue",
    res_gene$sig_plot & res_gene$direction == "maternal" ~ "red",
    TRUE                                                  ~ "grey"
  )
  res_gene$alpha_plot <- ifelse(res_gene$sig_plot, 1, 0.35)
  
  # ----- gene_id / symbol mapping -----
  res_nb_gene_plot <- res_gene |>
    dplyr::mutate(
      gene_id = sub("^(.*):[^:]+$", "\\1", feature_id),  # keep before last ":"
      gene_id = sub("\\..*$", "", gene_id)               # drop version
    ) |>
    dplyr::left_join(ensembl, by = "gene_id")
  res_nb_gene_plot$symbol[is.na(res_nb_gene_plot$symbol)] <- res_nb_gene_plot$feature_id[is.na(res_nb_gene_plot$symbol)]
  res_nb_gene_plot$symbol[res_nb_gene_plot$symbol == ""]  <- res_nb_gene_plot$feature_id[res_nb_gene_plot$symbol == ""]
  
  # ----- subtitle -----
  ratio_label_txt <- sprintf("%d:%d", round(thr_ratio*100), round((1-thr_ratio)*100))
  subtxt <- bquote(
    beta[POE] < 0 ~ ": maternal bias; " ~
      beta[POE] > 0 ~ ": paternal bias | " ~
      "Threshold: |" * beta[POE] * "| > " * .(round(thr_shift, 3)) ~
      " (major allele " * .(ratio_label_txt) * ")"
  )
  xlab <- expression(beta[POE])

  
  # Split data for proper layering
  non_sig_data <- subset(res_nb_gene_plot, !sig_plot)
  sig_data <- subset(res_nb_gene_plot, sig_plot)
  
  g <- ggplot2::ggplot(res_nb_gene_plot, ggplot2::aes(x = beta_shift, y = -log10(.data[[val_col]]))) +
    # Plot non-significant points first (grey, behind)
    ggplot2::geom_errorbar(orientation = "y",
      data = non_sig_data,
      ggplot2::aes(xmin = beta_lo_shift, xmax = beta_hi_shift),
      color = "grey", alpha = 0.35, width = 0.1, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = non_sig_data,
      color = "grey", alpha = 0.35, size = 1.6, show.legend = FALSE
    ) +
    # Plot significant points on top (colored, in front)
    ggplot2::geom_errorbar(orientation = "y",
      data = sig_data,
      ggplot2::aes(xmin = beta_lo_shift, xmax = beta_hi_shift, color = col_plot),
      alpha = 1, width = 0.1, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = sig_data,
      ggplot2::aes(color = col_plot),
      alpha = 1, size = 1.6, show.legend = FALSE
    ) +
    ggplot2::scale_color_identity() +
    # guides: dashed center at 50:50 (x=0), dotted at ± threshold
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, color = "black") +
    ggplot2::geom_vline(xintercept = c(-thr_shift, thr_shift), linetype = "dotted", color = "black") +
    ggplot2::geom_hline(yintercept = -log10(val_thresh), linetype = "dashed", color = "black") +
    ggrepel::geom_label_repel(
      data = sig_data,
      ggplot2::aes(label = symbol, color = col_plot),
      alpha = 1,
      fill = "white", label.size = 0.25, label.r = grid::unit(0.1, "lines"),
      box.padding = 1, point.padding = 0, show.legend = FALSE, max.overlaps = 200
    ) +
    ggplot2::labs(
      title = "Gene-level POE",
      subtitle = subtxt,
      x = xlab,
      y = bquote(-log[10] * "(q)")
    ) +
    ggplot2::scale_x_continuous(
      limits = x_limits,
      breaks = lower_breaks,
      labels = function(x) sprintf("%.1f", x),
      sec.axis = ggplot2::sec_axis(
        trans = ~ .,
        name  = "Allelic ratio (paternal : maternal)",
        breaks = x_ticks,
        labels = ratio_lbl
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray30"),
      axis.text.x.top = ggplot2::element_text(size = 10)
    )
  
  g
}




# ---- PLOTTING --------------------------------------------------------------
make_poe_diff_volcanos <- function(
    # Provide exactly one of the following pairs:
  res_tx   = NULL, tiers_tx  = NULL,      # transcript-level
  res_gene = NULL, tiers_gene = NULL,     # gene-level
  # Annotation
  tx2g = NULL,                            # required if transcript-level is provided
  ensembl,
  # Display controls
  exposure_label = NULL,
  fdr_hline = 0.10,
  min_abs_lfc_int = log(2),               # threshold on |beta_int| (log fold-change)
  label_tier = c("any","strict")
) {
  label_tier <- match.arg(label_tier)
  
  # --- choose level (exactly one) ---
  provided_tx   <- !is.null(res_tx)   && !is.null(tiers_tx)
  provided_gene <- !is.null(res_gene) && !is.null(tiers_gene)
  if (provided_tx == provided_gene) stop("Provide exactly one of (res_tx, tiers_tx) or (res_gene, tiers_gene).")
  level <- if (provided_tx) "tx" else "gene"
  if (level == "tx" && is.null(tx2g)) stop("tx2g is required when using transcript-level inputs.")
  if (missing(ensembl)) stop("ensembl annotation (gene_id -> symbol) must be supplied.")
  
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(ggrepel); library(grid)
  })
  
  # ---------- helpers ----------
  annotate_level <- function(res, tiers, level = c("tx","gene")) {
    level <- match.arg(level)
    res$gene_id <- NULL; res$symbol <- NULL
    
    # Use ONLY DIFF flags (enhance retired)
    res$diff_poe <- FALSE
    res$tier <- "ns"
    strict_sig <- is.logical(tiers$pass_diff_strict) & tiers$pass_diff_strict
    loose_sig  <- is.logical(tiers$pass_diff_loose)  & tiers$pass_diff_loose
    any_sig    <- strict_sig | loose_sig
    
    if (length(any_sig) == nrow(res)) {
      res$diff_poe[any_sig] <- TRUE
    } else if (length(any_sig) == 1L && isTRUE(any_sig)) {
      res$diff_poe[] <- TRUE
    }
    res$tier[strict_sig] <- "strict"
    res$tier[loose_sig & !strict_sig] <- "loose"
    
    if (level == "tx") {
      res <- dplyr::left_join(res, tx2g, by = c("feature_id" = "transcript_id"))
      res$gene_id <- sub("\\..*$", "", res$gene_id)
    } else {
      res$gene_id <- sub("\\..*$", "", res$feature_id)
    }
    dplyr::left_join(res, ensembl, by = "gene_id")
  }
  
  get_cond_type <- function(df) {
    u <- unique(na.omit(df$cond_type)); if (!length(u)) return(NA_character_)
    if (length(u) > 1L) u <- u[1L]; u
  }
  
  # ---------- DIFF volcano: gates on q_int and |beta_int| ----------
  make_diff_volcano <- function(df,
                                q_col  = "q_int",
                                title = "Differential POE (interaction)",
                                fdr_hline = 0.10,
                                min_abs_lfc_int = log(2),
                                label_selector = c("any","strict")) {
    label_selector <- match.arg(label_selector)
    cond_type <- get_cond_type(df)
    cond_phrase <- if (isTRUE(cond_type == "binary")) "group difference" else "per +1 SD exposure"
    
    # Gates
    abs_beta <- abs(df$beta_int)
    lfc_ok   <- is.finite(abs_beta) & (abs_beta >= min_abs_lfc_int)
    q_vals   <- df[[q_col]]
    pass_q   <- is.finite(q_vals) & (q_vals <= fdr_hline)
    pass     <- pass_q & lfc_ok
    
    # Colors
    col <- rep("grey70", nrow(df))
    col[pass & df$beta_int > 0] <- "steelblue3"   # shift toward paternal
    col[pass & df$beta_int < 0] <- "indianred3"   # shift toward maternal
    df$col_plot <- factor(col, levels = c("grey70","steelblue3","indianred3"))
    
    # Safe y
    q_safe <- pmax(q_vals, .Machine$double.xmin)
    df$.yqq <- -log10(q_safe)
    
    # Labels
    label_rows <- if (label_selector == "strict") {
      df[pass & df$tier == "strict", , drop = FALSE]
    } else {
      df[pass & df$diff_poe, , drop = FALSE]
    }
    
    # Plot
    p <- ggplot(df, aes(x = beta_int, y = .yqq)) +
      geom_point(aes(color = .data$col_plot), size = 1.2, alpha = 0.9) +
      scale_color_identity() +
      geom_errorbar(orientation = "y",
        aes(xmin = beta_int - 1.96*se_int, xmax = beta_int + 1.96*se_int, color = .data$col_plot),
        width = 0.03, alpha = 0.45
      ) +
      geom_vline(xintercept = c(-1, 1) * min_abs_lfc_int, linetype = "dotted") +
      geom_hline(yintercept = -log10(fdr_hline), linetype = "dashed") +
      geom_label_repel(
        data = label_rows,
        aes(x = beta_int, y = .yqq, label = symbol, color = .data$col_plot),
        fill = "white", label.size = 0.25, label.r = unit(0.08, "lines"),
        box.padding = 0.8, point.padding = 0.1,
        show.legend = FALSE, max.overlaps = 200
      ) +
      labs(
        title = title,
        subtitle = bquote(
          atop(
            beta[int] < 0 ~ ": shifts toward maternal; " ~ beta[int] > 0 ~ ": shifts toward paternal",
            "Thresholds: " * q[interaction] <= .(fdr_hline) ~ 
              " and |logFC| ≥ log(" * .(format(exp(min_abs_lfc_int), digits = 3)) * ")"
          )
        ),
        x = expression(beta[interaction]),
        y = bquote(-log[10](q[interaction]))
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray30", hjust = 0)
      ) +
      coord_cartesian(
        xlim = range(df$beta_int, na.rm = TRUE)
      )
    
    # y scaling
    y_break_at <- 10
    y_max <- max(df$.yqq[is.finite(df$.yqq)], na.rm = TRUE)
    need_break <- is.finite(y_max) &&
      y_max > (y_break_at + 5) &&
      sum(df$.yqq > y_break_at, na.rm = TRUE) >= 1 &&
      stats::quantile(df$.yqq, 0.95, na.rm = TRUE) <= y_break_at
    
    if (need_break) {
      top_end <- ceiling(y_max)
      bot_frac <- 0.75
      trans_fun <- function(y) {
        ifelse(y <= y_break_at,
               y * (bot_frac / y_break_at),
               bot_frac + (y - y_break_at) * ((1 - bot_frac) / (top_end - y_break_at)))
      }
      inv_fun <- function(x) {
        ifelse(x <= bot_frac,
               x * (y_break_at / bot_frac),
               y_break_at + (x - bot_frac) * ((top_end - y_break_at) / (1 - bot_frac)))
      }
      br_trans <- scales::trans_new("ybreak", trans_fun, inv_fun)
      
      b_bot <- scales::pretty_breaks(n = 5)(c(0, y_break_at))
      b_top <- scales::pretty_breaks(n = 3)(c(y_break_at, top_end))
      y_breaks <- sort(unique(c(b_bot[b_bot <= y_break_at], b_top[b_top >= y_break_at])))
      
      p <- p + scale_y_continuous(trans = br_trans, breaks = y_breaks, limits = c(0, top_end))
    } else {
      y_lim <- ceiling(max(10, y_max))
      p <- p + scale_y_continuous(limits = c(0, y_lim))
    }
    
    p
  }
  
  # ---------- build annotated table and plot for the chosen level ----------
  if (level == "tx") {
    a <- annotate_level(res_tx, tiers_tx, level = "tx")
    title_base <- "Isoform-level differential POE"
    if (!is.null(exposure_label) && nzchar(exposure_label)) title_base <- paste(title_base, exposure_label)
    p_diff <- make_diff_volcano(a, title = title_base,
                                fdr_hline = fdr_hline,
                                min_abs_lfc_int = min_abs_lfc_int,
                                label_selector = label_tier)
    return(list(diff_tx = p_diff))
  } else {
    a <- annotate_level(res_gene, tiers_gene, level = "gene")
    title_base <- "Gene-level differential POE"
    if (!is.null(exposure_label) && nzchar(exposure_label)) title_base <- paste(title_base, exposure_label)
    p_diff <- make_diff_volcano(a, title = title_base,
                                fdr_hline = fdr_hline,
                                min_abs_lfc_int = min_abs_lfc_int,
                                label_selector = label_tier)
    return(list(diff_gene = p_diff))
  }
}



# Forest plots
poe_forest_absbeta <- function(
    res,
    ensembl, tx2g = NULL,
    tx = FALSE,
    pfas_label = "PFAS",
    comparison_level = c("cohort","subset"),
    fdr_hline = 0.10,
    min_abs_lfc_int = log(1.5),
    label_tier = c("any","strict"),
    alpha = 0.05,
    top_n = Inf,
    plot_exposure_effect = TRUE  # if FALSE, only plot baseline (ref), not after
) {
  suppressPackageStartupMessages({
    library(dplyr); library(forcats); library(ggplot2); library(scales); library(tidyr)
  })
  comparison_level <- match.arg(comparison_level)
  label_tier <- match.arg(label_tier)
  
  `%||%` <- function(a,b) if (is.null(a)) b else a
  inv_logit <- function(x) 1/(1+exp(-x))
  clamp01 <- function(p, eps=1e-9) pmax(pmin(p, 1-eps), eps)
  
  # Determine which columns to use
  required_base <- c("q_int")
  if (!all(required_base %in% names(res))) {
    stop("Missing required column: q_int")
  }
  q_col_suffix <- ""

  
  # --- make sure diff flags exist ---
  if (!"pass_diff_loose" %in% names(res))  res$pass_diff_loose  <- FALSE
  if (!"pass_diff_strict" %in% names(res)) res$pass_diff_strict <- FALSE
  
  # --- significance mask (same as volcano) ---
  res$diff_poe <- (res$pass_diff_loose == TRUE) | (res$pass_diff_strict == TRUE)
  res$tier <- dplyr::case_when(
    res$pass_diff_strict ~ "strict",
    res$pass_diff_loose  ~ "loose",
    TRUE ~ "ns"
  )
  
  abs_beta <- abs(res$beta_int)
  lfc_ok   <- is.finite(abs_beta) & (abs_beta >= min_abs_lfc_int)
  pass_q   <- is.finite(res$q_int) & (res$q_int <= fdr_hline)
  pass_sig <- pass_q & lfc_ok
  label_mask <- if (label_tier == "strict") res$tier == "strict" else res$diff_poe
  sig_ids <- res$feature_id[ pass_sig & label_mask ]
  df <- res %>% filter(feature_id %in% sig_ids)
  if (!nrow(df)) stop("No features pass volcano significance; nothing to plot.")
  
  # ----- gene symbols -----
  if (tx) {
    df <- df %>%
      left_join(tx2g, by = c("feature_id" = "transcript_id")) %>%
      mutate(gene_id = sub("\\..*$","", gene_id))
  } else {
    df$gene_id <- sub("\\..*$","", df$feature_id)
  }
  df <- df %>% left_join(ensembl, by="gene_id") %>%
    mutate(label = ifelse(!is.na(symbol) & nzchar(symbol), symbol, feature_id))
  
  # ----- choose ref/after columns based on calibration and comparison level -----
  z <- stats::qnorm(1 - alpha/2)
  
  if (comparison_level == "cohort") {
    ref_beta_col   <- "beta_POE_at_ref"
    ref_se_col     <- "se_POE_at_ref"
    after_beta_col <- "beta_POE_after"
    after_lo_col   <- "beta_POE_after_ci_lo"
    after_hi_col   <- "beta_POE_after_ci_hi"
    legend_ref     <- paste0("Cohort μ ", pfas_label)
    legend_after   <- paste0("Cohort μ + 1 SD ", pfas_label)
  } else {
    ref_beta_col   <- "beta_POE_subset"
    ref_se_col     <- "se_beta_POE_subset"
    after_beta_col <- "beta_POE_after_subset"
    after_lo_col   <- "beta_POE_after_subset_ci_lo"
    after_hi_col   <- "beta_POE_after_subset_ci_hi"
    legend_ref     <- paste0("Subset μ ", pfas_label)
    legend_after   <- paste0("Subset μ + 1 SD ", pfas_label)
  }

  
  # Check if required columns exist
  required_cols <- c(ref_beta_col, ref_se_col)
  if (plot_exposure_effect) {
    required_cols <- c(required_cols, after_beta_col)
  }
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns for forest plot: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create long format data
  if (plot_exposure_effect) {
    long <- df %>%
      transmute(
        feature_id, label,
        beta_ref   = .data[[ref_beta_col]],
        lo_ref     = .data[[ref_beta_col]] - z * .data[[ref_se_col]],
        hi_ref     = .data[[ref_beta_col]] + z * .data[[ref_se_col]],
        beta_after = .data[[after_beta_col]],
        lo_after   = if (after_lo_col %in% names(df)) .data[[after_lo_col]] else NA_real_,
        hi_after   = if (after_hi_col %in% names(df)) .data[[after_hi_col]] else NA_real_,
        abs_beta_int = abs(beta_int),
        q_int = q_int
      ) %>%
      # If CI columns don't exist, compute them
      mutate(
        lo_after = if (is.na(lo_after[1])) beta_after - z * .data[[ref_se_col]] else lo_after,
        hi_after = if (is.na(hi_after[1])) beta_after + z * .data[[ref_se_col]] else hi_after
      ) %>%
      pivot_longer(
        c(beta_ref, lo_ref, hi_ref, beta_after, lo_after, hi_after),
        names_to = c(".value","group"),
        names_pattern = "^(beta|lo|hi)_(ref|after)$"
      ) %>%
      filter(is.finite(beta), is.finite(lo), is.finite(hi))
  } else {
    # Only plot baseline (ref)
    long <- df %>%
      transmute(
        feature_id, label,
        beta = .data[[ref_beta_col]],
        lo   = .data[[ref_beta_col]] - z * .data[[ref_se_col]],
        hi   = .data[[ref_beta_col]] + z * .data[[ref_se_col]],
        group = "ref",
        abs_beta_int = abs(beta_int),
        q_int = q_int
      ) %>%
      filter(is.finite(beta), is.finite(lo), is.finite(hi))
  }
  
  # No additional calibration needed - we're already using the appropriate columns
  long <- long %>% mutate(
    beta_hat = beta,
    lo_hat   = lo,
    hi_hat   = hi
  )
  
  # ----- order (|β_int| desc, q asc) -----
  ord <- long %>% distinct(feature_id, label, abs_beta_int, q_int) %>%
    arrange(desc(abs_beta_int), q_int, label) %>% head(top_n)
  keep_ids <- ord$feature_id
  plot_df <- long %>%
    filter(feature_id %in% keep_ids) %>%
    mutate(label = factor(label, levels = unique(rev(ord$label))))
  
  if (plot_exposure_effect) {
    plot_df$Group <- factor(plot_df$group, levels = c("ref","after"),
                            labels = c(legend_ref, legend_after))
  } else {
    plot_df$Group <- legend_ref
  }
  
  # ----- x-axis ticks (centered at beta_POE = 0) -----
  x_range <- range(c(plot_df$beta_hat, plot_df$lo_hat, plot_df$hi_hat), na.rm = TRUE)
  
  # standard breaks
  s_breaks <- pretty(x_range, n = 5)
  xlab <- expression(beta[POE])

  
  pd <- position_dodge(width = 0.6)
  
  p <- ggplot(plot_df, aes(x = beta_hat, y = label, color = Group)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +  # Always center at 0 (50:50)
    geom_errorbar(orientation = "y", aes(xmin = lo_hat, xmax = hi_hat), width = 0.25, position = pd) +
    geom_point(size = 2.6, position = pd) +
    scale_x_continuous(
      name = xlab,
      breaks = s_breaks,
      labels = function(x) sprintf("%.2f", x)
    ) +
    labs(
      title = if (plot_exposure_effect) {
        paste0("POE at Mean vs. Mean + 1 SD ", pfas_label)
      } else {
        paste0("POE at Mean ", pfas_label)
      },
      subtitle = paste0("Features with significant exposure x POE effects: q ≤", fdr_hline,
                        " & |logFC| ≥", format(exp(min_abs_lfc_int), digits = 3),
                        ""),
      y = NULL, color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0),
      plot.subtitle = element_text(hjust = 0)
    )
  
  # Remove legend if only plotting baseline
  if (!plot_exposure_effect) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}






#### Bar chart
poe_strict_barplot <- function(
    poe_res_tx,
    poe_res_gene,
    subset_type = c("mat","child"),
    title = "Number of features with strict POE",
    tx2g.nochr = NULL
) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(purrr); library(stringr); library(forcats)
  })
  
  subset_type <- match.arg(subset_type)
  
  # Helper to summarize one list
  summarize_list <- function(lst, feature_level) {
    keep <- grep(paste0("_", subset_type, "$"), names(lst), value = TRUE)
    if (length(keep) == 0) return(NULL)
    
    map_dfr(keep, function(nm) {
      df <- lst[[nm]]
      
      if (feature_level == "Genes with Isoforms" && !is.null(tx2g.nochr)) {
        # Map transcripts to genes and count unique genes
        df <- df %>%
          left_join(tx2g.nochr, by = c("feature_id" = "transcript_id")) %>%
          filter(pass_diff_strict %in% TRUE) %>%
          distinct(gene_id)
        
        n_strict <- nrow(df)
      } else {
        # Just count features (genes already)
        n_strict <- sum(df$pass_diff_strict %in% TRUE, na.rm = TRUE)
      }
      
      tibble(
        condition = nm,
        feature_level = feature_level,
        n_strict = n_strict
      )
    })
  }
  
  # Collect results
  df_all <- bind_rows(
    summarize_list(poe_res_tx, "Genes with Isoforms"),
    summarize_list(poe_res_gene, "Genes")
  )
  
  if (nrow(df_all) == 0) stop("No matching data for subset_type = ", subset_type)
  
  # Make condition labels nicer (remove suffix)
  df_all <- df_all %>%
    mutate(condition = str_replace(condition, paste0("_", subset_type, "$"), ""))
  
  # Order x-axis by Genes with Isoforms counts
  iso_counts <- df_all %>%
    filter(feature_level == "Genes with Isoforms") %>%
    select(condition, n_strict)
  
  df_all <- df_all %>%
    left_join(iso_counts, by = "condition", suffix = c("", "_iso")) %>%
    mutate(condition = forcats::fct_reorder(condition, n_strict_iso, .desc = TRUE))
  
  # Plot
  ggplot(df_all, aes(x = condition, y = n_strict, fill = feature_level)) +
    geom_col(position = "dodge") +
    labs(
      title = title,
      x = NULL,
      y = "Count",
      fill = NULL,
      subtitle = "Threshold on β = ± 0.405 and FDR < 10%"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0),
      legend.position = "bottom"
    )
}








########## DOSE-RESPONSE CURVES ##########
# =============================================================================
# Function to plot dose-response curves for a given condition
# =============================================================================
plot_dose_response <- function(
    res,
    filter_tier,
    ensembl,
    metadata,                 # REQUIRED
    tx2g = NULL,
    level = "gene",
    condition_name = "",
    centering = "cohort",     # "cohort" or "subset"
    cohort_pfas = NULL,
    output_pdf_path = NULL,
    plot_title = NULL,
    y_range_min = 0,
    y_range_max = 100,
    exclude = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2').")
  if (!requireNamespace("tidyr", quietly = TRUE))  stop("Please install.packages('tidyr').")
  level     <- match.arg(level,     c("gene","transcript"))
  centering <- match.arg(centering, c("cohort","subset"))
  if (level == "transcript" && is.null(tx2g)) stop("tx2g mapping is required when level='transcript'")
  if (!is.null(tx2g) && !is.data.frame(tx2g)) tx2g <- as.data.frame(tx2g)
  if (is.null(metadata)) stop("'metadata' is required.")
  
  # --- helpers ---
  dose_range <- -8:8
  dose_suf   <- ifelse(dose_range < 0, paste0("m", abs(dose_range)), as.character(dose_range))
  inv_logit  <- function(x) 1/(1+exp(-x))
  logit      <- function(p) log(p/(1-p))
  clamp01    <- function(p, eps=1e-9) pmax(pmin(p, 1-eps), eps)
  
  # resolve condition column name
  condition_col <- if (nzchar(condition_name)) {
    condition_name
  } else if ("condition" %in% names(res)) {
    unique(res$condition)[1]
  } else stop("Provide 'condition_name' or include a 'condition' column in 'res'.")
  
  # helper to pick column from a df
  pick_col <- function(df) {
    cands <- c(condition_col, paste0(condition_col, "_mat"), paste0(condition_col, "_child"))
    col <- cands[cands %in% names(df)][1]
    if (is.na(col)) stop("Could not find '", condition_col, "' (or *_mat/_child).")
    v <- df[[col]]; v <- v[is.finite(v)]
    if (!length(v)) stop("No finite values in '", col, "'.")
    list(values = v, colname = col)
  }
  
  # cohort & subset distributions (for x annotations only)
  coh_src <- if (!is.null(cohort_pfas)) pick_col(cohort_pfas) else pick_col(metadata)
  mu_cohort <- mean(coh_src$values); sd_cohort <- stats::sd(coh_src$values)
  
  sub_src <- pick_col(metadata)
  subset_min <- min(sub_src$values); subset_max <- max(sub_src$values)
  A <- mean(sub_src$values); B <- stats::sd(sub_src$values)
  if (!is.finite(B) || B <= 0) stop("Subset SD (computed from metadata) must be > 0.")
  
  # filter tier
  filter_col <- switch(
    filter_tier,
    "strict" = "pass_diff_strict",
    "loose"  = "pass_diff_loose",
    "any"    = "pass_diff_any",
    stop("Invalid filter_tier")
  )
  if (filter_tier == "any" && !"pass_diff_any" %in% names(res)) {
    res$pass_diff_any <- (res$pass_diff_strict %in% TRUE) | (res$pass_diff_loose %in% TRUE)
  }
  plot_data <- res[res[[filter_col]] %in% TRUE, ]
  if (nrow(plot_data) == 0) return(invisible(NULL))
  
  # --- build long_data depending on centering ---
  if (centering == "subset") {
    beta_cols  <- paste0("beta_POE_subset_at_dose_", dose_suf)
    ci_lo_cols <- paste0(beta_cols, "_ci_lo")
    ci_hi_cols <- paste0(beta_cols, "_ci_hi")
    if (!all(beta_cols %in% names(plot_data))) stop("Subset grid columns missing in 'res'.")
    
    long_data <- do.call(rbind, lapply(seq_along(dose_range), function(i) {
      data.frame(
        feature_id     = plot_data$feature_id,
        x              = dose_range[i],  # subset SD units
        beta_POE       = plot_data[[beta_cols[i]]],
        beta_POE_ci_lo = plot_data[[ci_lo_cols[i]]],
        beta_POE_ci_hi = plot_data[[ci_hi_cols[i]]]
      )
    }))
    
    x_range_min <- (subset_min - A) / B
    x_range_max <- (subset_max - A) / B
    x_axis_label <- "Dose (subset SD units)"
    x_axis_source <- sprintf("X-range: subset (metadata '%s') [%.2f, %.2f]", sub_src$colname, x_range_min, x_range_max)
    show_subset_dotted <- FALSE
    
    # convert cohort mean/sd into subset units for purple/red
    mu_line <- (mu_cohort - A) / B
    sd_line <- (mu_cohort + sd_cohort - A) / B
    
  } else {
    beta_cols  <- paste0("beta_POE_at_dose_", dose_suf)
    ci_lo_cols <- paste0(beta_cols, "_ci_lo")
    ci_hi_cols <- paste0(beta_cols, "_ci_hi")
    if (!all(beta_cols %in% names(plot_data))) stop("Cohort grid columns missing in 'res'.")
    
    long_data <- do.call(rbind, lapply(seq_along(dose_range), function(i) {
      data.frame(
        feature_id     = plot_data$feature_id,
        x              = dose_range[i],  # cohort units
        beta_POE       = plot_data[[beta_cols[i]]],
        beta_POE_ci_lo = plot_data[[ci_lo_cols[i]]],
        beta_POE_ci_hi = plot_data[[ci_hi_cols[i]]]
      )
    }))
    
    if (!is.null(cohort_pfas)) {
      x_range_min <- min(coh_src$values); x_range_max <- max(coh_src$values)
      x_axis_source <- sprintf("X-range: cohort_pfas '%s' [%.2f, %.2f]", coh_src$colname, x_range_min, x_range_max)
      show_subset_dotted <- TRUE
    } else {
      x_range_min <- min(sub_src$values); x_range_max <- max(sub_src$values)
      x_axis_source <- sprintf("X-range: metadata '%s' [%.2f, %.2f]", sub_src$colname, x_range_min, x_range_max)
      show_subset_dotted <- FALSE
    }
    x_axis_label <- "Dose (cohort units)"
    
    mu_line <- mu_cohort
    sd_line <- mu_cohort + sd_cohort
  }
  
  # trim to range
  long_data <- long_data[long_data$x >= x_range_min & long_data$x <= x_range_max, ]
  if (!nrow(long_data)) return(invisible(NULL))
  
  # display names
  long_data$display_name <- convert_to_symbols(long_data$feature_id, ensembl, tx2g, level)
  if (!is.null(exclude)) long_data <- long_data[!long_data$display_name %in% exclude, ]
  
  # --- OPTIONAL RECALIBRATION on combined dose response ---
  calib_note <- ""

  
  # proportions (paternal %) from betas
  to_percent <- function(b) 100 * inv_logit(b)
  long_data$paternal_prop <- to_percent(long_data$beta_POE)
  long_data$ci_lo_prop    <- to_percent(long_data$beta_POE_ci_lo)
  long_data$ci_hi_prop    <- to_percent(long_data$beta_POE_ci_hi)
  
  # title & subtitle
  if (is.null(plot_title)) {
    plot_title <- paste0("Dose-Response for Allelic Ratios (", condition_col, " – ", centering, ")")
  }
  if(centering=="cohort"){
    plot_subtitle <- paste0(
      "Ribbons = 95% CI. Purple = subset mean; Red = subset mean + 1 SD.\n",
      x_axis_source, calib_note
    )
  } else {
    plot_subtitle <- paste0(x_axis_source, calib_note)
  }
  if (show_subset_dotted) {
    plot_subtitle <- paste0(plot_subtitle, " Black dotted = subset range [",
                            round(subset_min,2), ", ", round(subset_max,2), "].")
  }
  
  # x-axis breaks (integers within range)
  breaks_x <- seq(ceiling(x_range_min), floor(x_range_max), by = 1)
  
  # plot
  p <- ggplot2::ggplot(long_data, ggplot2::aes(x = x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lo_prop, ymax = ci_hi_prop), alpha = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = paternal_prop), linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.8) +
    ggplot2::facet_wrap(~ display_name, scales = "free_y") +
    ggplot2::labs(title = plot_title, subtitle = plot_subtitle,
                  x = x_axis_label, y = "Paternal:Maternal Expression") +
    ggplot2::scale_y_continuous(limits = c(y_range_min, y_range_max),
                                breaks = seq(0, 100, by = 20),
                                labels = function(y) paste0(y, ":", 100 - y)) +
    ggplot2::scale_x_continuous(
      limits = c(x_range_min, x_range_max),
      breaks = breaks_x,
      labels = scales::label_number(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0,0))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (show_subset_dotted) {
    p <- p +
      ggplot2::geom_vline(xintercept = subset_min, linetype = "dotted", color = "black") +
      ggplot2::geom_vline(xintercept = subset_max, linetype = "dotted", color = "black")
  }
  if(centering=="cohort"){
    p <- p +
      ggplot2::geom_vline(xintercept = mu_line, color = "purple", linewidth = 0.7) +
      ggplot2::geom_vline(xintercept = sd_line, color = "red", linewidth = 0.7)
  }
  
  if (!is.null(output_pdf_path)) {
    ggplot2::ggsave(output_pdf_path, plot = p, device = "pdf", width = 8, height = 6)
  } else {
    print(p)
  }
  invisible(p)
}

