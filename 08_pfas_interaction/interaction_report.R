# =============================================================================
# Stage 08 | interaction_report.R
# PFAS / GDM x parent-of-origin interaction results: tiers, volcano, p-value,
# QQ and forest plots per condition and level; hit summaries; overlap with
# known imprinted genes; counts of strict hits per PFAS (maternal, child).
#
# In : MODEL_DIR/poe_by_condition_<level>_<condition>.rds, MODEL_DIR/model_metadata.rds,
#      TX_CLASSIF, PLACENTA_DB, IMPRINT_DB
# Out: MODEL_DIR/poe_diff_tiers.rds          (input to dose_response.R)
#      WORK_DIR/poe_interaction_all.tsv.gz    all features x conditions x levels
#      WORK_DIR/poe_interaction_hits.tsv      q_int < 0.10 and |beta_int| > log(1.5)
#      WORK_DIR/poe_interaction_hit_summary.tsv, poe_interaction_gene_tx_overlap.tsv
#      FIG_DIR/interaction/{gene,tx}_diffPOE_<cond>_BetaInt{,_pvals,_pQQ}.png,
#              {gene,tx}_forest_<cond>.png, pfas_bar_{mat,child}.png
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
load_poe_library()

mm <- readRDS(file.path(P("MODEL_DIR"), "model_metadata.rds"))
conds <- mm$conditions
tx2g <- as.data.frame(read_tx2g())
ensembl <- read_ensembl()
imprintdb <- fread(P("IMPRINT_DB"))
fdr <- 0.10
min_abs_beta_int <- log(1.5)
fdir <- fig("interaction"); dir.create(fdir, showWarnings = FALSE)
ifig <- function(name) file.path(fdir, name)

read_level <- function(lv) lapply(setNames(conds, conds), function(cn)
  readRDS(file.path(P("MODEL_DIR"), sprintf("poe_by_condition_%s_%s.rds", lv, cn))))
res <- list(gene = read_level("gene"), transcript = read_level("transcript"))
tiers <- lapply(res, function(l) lapply(l, classify_poe_nbmm_diff))
saveRDS(tiers, file.path(P("MODEL_DIR"), "poe_diff_tiers.rds"))

exposure_label <- function(cn) {
  if (mm$cinfo[[cn]]$type == "binary") return(paste("by", cn))
  src <- if (grepl("_mat$", cn)) "maternal" else if (grepl("_child$", cn)) "child" else ""
  paste("by 1 SD", src, sub("_(mat|child)$", "", cn))
}

for (lv in names(res)) for (cn in conds) {
  message("plots: ", lv, " ", cn)
  r <- res[[lv]][[cn]]; t <- tiers[[lv]][[cn]]
  pre <- if (lv == "gene") "gene" else "tx"
  if (!any(is.finite(r$p_int))) { message("  no fitted features; skipped"); next }
  vol <- if (lv == "gene") {
    make_poe_diff_volcanos(res_gene = r, tiers_gene = t, tx2g = tx2g, ensembl = ensembl,
                           exposure_label = exposure_label(cn), min_abs_lfc_int = min_abs_beta_int,
                           fdr_hline = fdr, label_tier = "any")$diff_gene
  } else {
    make_poe_diff_volcanos(res_tx = r, tiers_tx = t, tx2g = tx2g, ensembl = ensembl,
                           exposure_label = exposure_label(cn), min_abs_lfc_int = min_abs_beta_int,
                           fdr_hline = fdr, label_tier = "any")$diff_tx
  }
  ggsave(ifig(sprintf("%s_diffPOE_%s_BetaInt.png", pre, cn)), vol, width = 9, height = 6, dpi = 300, bg = "white")

  gp <- ggplot(t, aes(p_int)) + geom_histogram(bins = 100, boundary = 0, closed = "left") +
    labs(x = "Interaction p-value", y = "Count", title = paste(lv, "level POE x", cn)) +
    theme_minimal(base_size = 12)
  ggsave(ifig(sprintf("%s_diffPOE_%s_BetaInt_pvals.png", pre, cn)), gp, width = 6.5, height = 4.5, dpi = 300, bg = "white")
  qq <- make_nb_int_qq(t, mode = "p", which = "interaction") + ggtitle(paste(lv, "level POE x", cn))
  ggsave(ifig(sprintf("%s_diffPOE_%s_BetaInt_pQQ.png", pre, cn)), qq, width = 6.5, height = 4.5, dpi = 300, bg = "white")

  if (mm$cinfo[[cn]]$type == "continuous") {
    fo <- try(poe_forest_absbeta(res = r, min_abs_lfc_int = min_abs_beta_int, label_tier = "any",
                                 ensembl = ensembl, tx2g = tx2g, tx = lv == "transcript",
                                 comparison_level = "cohort", pfas_label = cn), silent = TRUE)
    if (inherits(fo, "ggplot"))
      ggsave(ifig(sprintf("%s_forest_%s.png", pre, cn)), fo, width = 7.5, height = 6.5, dpi = 300, bg = "white")
    else message("  forest plot skipped: ", if (inherits(fo, "try-error")) conditionMessage(attr(fo, "condition")) else "no hits")
  }
}

# ---- tables ---------------------------------------------------------------------------
all_res <- rbindlist(lapply(res, function(l) rbindlist(l, fill = TRUE)), fill = TRUE)
fwrite(all_res, out("poe_interaction_all.tsv.gz"), sep = "\t")

hits <- all_res[is.finite(q_int) & q_int < fdr & abs(beta_int) > min_abs_beta_int,
                .(level, condition, feature_id, beta_int, se_int, p_int, q_int,
                  beta_POE_before, beta_POE_after, delta_mag = abs(beta_POE_after) - abs(beta_POE_before))]
hits[, gene_id := ifelse(level == "transcript", tx2g$gene_id[match(feature_id, tx2g$transcript_id)], feature_id)]
hits[, symbol := ifelse(level == "transcript",
                        convert_to_symbols(feature_id, ensembl, tx2g, "transcript"),
                        convert_to_symbols(feature_id, ensembl, level = "gene"))]
hits[, imprinted := sub("\\..*$", "", gene_id) %in% imprintdb$ensembl_gene_id]
fwrite(hits, out("poe_interaction_hits.tsv"), sep = "\t")

qs <- function(v) list(n = sum(is.finite(v)), mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE),
                       min = suppressWarnings(min(v, na.rm = TRUE)),
                       q25 = quantile(v, 0.25, na.rm = TRUE, names = FALSE),
                       median = median(v, na.rm = TRUE),
                       q75 = quantile(v, 0.75, na.rm = TRUE, names = FALSE),
                       max = suppressWarnings(max(v, na.rm = TRUE)))
hit_summary <- rbind(hits[, c(statistic = "beta_int", qs(beta_int)), by = .(level, condition)],
                     hits[, c(statistic = "delta_abs_beta_POE", qs(delta_mag)), by = .(level, condition)])
fwrite(hit_summary[order(level, condition, statistic)], out("poe_interaction_hit_summary.tsv"), sep = "\t")
# genes with hits at gene level, through transcripts, and at both
ov <- rbindlist(lapply(conds, function(cn) {
  g <- unique(hits[condition == cn & level == "gene", gene_id])
  t <- unique(hits[condition == cn & level == "transcript", gene_id])
  data.table(condition = cn, genes_gene_level = length(g), genes_via_transcripts = length(t),
             genes_both = length(intersect(g, t)))
}))
fwrite(ov, out("poe_interaction_gene_tx_overlap.tsv"), sep = "\t")

for (st in c("mat", "child")) {
  if (!any(grepl(paste0("_", st, "$"), conds))) next
  bar <- poe_strict_barplot(res$transcript, res$gene, tx2g.nochr = tx2g, subset_type = st,
                            title = sprintf("Features with significant PFAS x POE interactions (%s)",
                                            if (st == "mat") "maternal" else "child"))
  ggsave(fig(sprintf("pfas_bar_%s.png", st)), bar + theme(plot.subtitle = element_blank()),
         width = 3.5, height = 5.75, dpi = 300)
}
message("interaction report done")
