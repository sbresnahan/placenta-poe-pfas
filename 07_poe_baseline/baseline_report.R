# =============================================================================
# Stage 07 | baseline_report.R
# Cohort-level POE results: tiers, volcano / forest / p-value / QQ /
# calibration plots, transcript-vs-gene overlap and enrichment of known
# imprinted genes among POE genes.
#
# In : MODEL_DIR/poe_nbmm_results_{gene,transcript}.rds, MODEL_INPUTS,
#      TX_CLASSIF, PLACENTA_DB, IMPRINT_DB
# Out: MODEL_DIR/poe_tiers_{gene,transcript}.rds
#      WORK_DIR/poe_baseline_significant_{gene,transcript}.tsv
#      WORK_DIR/poe_baseline_overlap.tsv, imprinting_enrichment.tsv,
#               imprinting_overlap_genes.tsv
#      FIG_DIR/beta_by_majority_allele.png, {isoPOE,genePOE}{,_Bpdist,_pQQ,_forest,_calibration}.png
# Significance: BH q < 0.10 and the 95% CI of beta outside +/- log(1.5)
# (60:40 allelic ratio).
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
load_poe_library()

res_nb_tx   <- readRDS(file.path(P("MODEL_DIR"), "poe_nbmm_results_transcript.rds"))
res_nb_gene <- readRDS(file.path(P("MODEL_DIR"), "poe_nbmm_results_gene.rds"))
inp <- readRDS(P("MODEL_INPUTS"))
poMat.gene <- inp$poMat.gene; poMat.tx <- inp$poMat.tx
tx2g    <- as.data.frame(read_tx2g())
ensembl <- read_ensembl()
q_thresh <- 0.10
ratio_thr <- 0.60

# ---- tiers -------------------------------------------------------------------------
tiers.tx   <- classify_poe_nbmm(res_nb_tx,   val_thresh = q_thresh, require_CI_low = TRUE)
tiers.gene <- classify_poe_nbmm(res_nb_gene, val_thresh = q_thresh, require_CI_low = TRUE)
saveRDS(tiers.tx,   file.path(P("MODEL_DIR"), "poe_tiers_transcript.rds"))
saveRDS(tiers.gene, file.path(P("MODEL_DIR"), "poe_tiers_gene.rds"))
tiers_all <- function(t) bind_rows(mutate(t$significant, sig_plot = TRUE),
                                   mutate(t$not_significant, sig_plot = FALSE))
tiers_all_tx <- tiers_all(tiers.tx); tiers_all_gene <- tiers_all(tiers.gene)

sig_tx <- tiers.tx$significant %>% left_join(tx2g, by = c("feature_id" = "transcript_id"))
sig_tx$symbol <- convert_to_symbols(sig_tx$feature_id, ensembl, tx2g, "transcript")
sig_gene <- tiers.gene$significant
sig_gene$symbol <- convert_to_symbols(sig_gene$feature_id, ensembl, level = "gene")
fwrite(sig_tx,   out("poe_baseline_significant_transcript.tsv"), sep = "\t")
fwrite(sig_gene, out("poe_baseline_significant_gene.tsv"), sep = "\t")

# ---- transcript vs gene level (q < 0.1, |beta| > log(1.5)) ---------------------------
simple_sig <- function(r) r$feature_id[is.finite(r$q) & r$q < q_thresh & abs(r$beta) > log(1.5)]
sig.tx.ids <- simple_sig(res_nb_tx)
sig.tx.g   <- unique(tx2g$gene_id[tx2g$transcript_id %in% sig.tx.ids])
sig.gene   <- simple_sig(res_nb_gene)
gene_only_tx <- setdiff(sig.tx.g, sig.gene)
fwrite(data.table(
  quantity = c("transcripts significant", "genes with a significant transcript",
               "genes significant at gene level", "genes significant at both levels",
               "genes significant only through transcripts",
               "transcripts of genes significant only through transcripts"),
  n = c(length(sig.tx.ids), length(sig.tx.g), length(sig.gene), length(intersect(sig.tx.g, sig.gene)),
        length(gene_only_tx),
        sum(sig.tx.ids %in% tx2g$transcript_id[tx2g$gene_id %in% gene_only_tx]))),
  out("poe_baseline_overlap.tsv"), sep = "\t")

# ---- ratio <-> |beta| ------------------------------------------------------------------
df <- data.frame(ratio = seq(0.50, 0.995, by = 0.001))
df$abs_beta <- ratio_to_abs_beta(df$ratio)
cut_likely <- ratio_to_abs_beta(ratio_thr)
folds <- c(1, 1.5, 2, 3, 5, 10, 20)
g3 <- ggplot(df, aes(ratio, abs_beta)) +
  geom_line(linewidth = 1, color = "#0279EE") +
  geom_vline(xintercept = ratio_thr, linetype = "dashed") +
  geom_hline(yintercept = cut_likely, linetype = "dotted") +
  scale_x_continuous(name = "Majority-allele ratio", limits = c(0.50, 0.995),
                     breaks = c(0.50, 0.60, 0.70, 0.80, 0.90, 0.95), expand = c(0, 0),
                     sec.axis = sec_axis(~ . / (1 - .), breaks = folds, labels = paste0(folds, "x"),
                                         name = "Fold over 50:50")) +
  labs(y = "|beta|", title = "Majority-allele ratio and |beta|") +
  theme_classic(base_size = 12)
ggsave(fig("beta_by_majority_allele.png"), g3, width = 8, height = 5.5, dpi = 300)

# ---- per-level plots -------------------------------------------------------------------
level_plots <- function(res, tiers, tiers_all, poMat, lv, pre) {
  vol <- if (lv == "transcript") plot_isoform_poe_volcano(tiers_all, tx2g, ensembl, ratio_thr, q_thresh)
         else plot_gene_poe_volcano(tiers_all, ensembl, ratio_thr, q_thresh)
  ggsave(fig(paste0(pre, ".png")), vol, width = 9.5, height = if (lv == "transcript") 4.2 else 7, dpi = 300)

  gp <- ggplot(res, aes(p)) + geom_histogram(bins = 100, boundary = 0, closed = "left") +
    labs(x = "p-value", y = "Count", title = paste0(tools::toTitleCase(lv), "-level POE p-values")) +
    theme_minimal(base_size = 12)
  ggsave(fig(paste0(pre, "_Bpdist.png")), gp, width = 6.5, height = 4.5, dpi = 300)
  ggsave(fig(paste0(pre, "_pQQ.png")), make_glm_qq(res, mode = "p") + ggtitle(paste(lv, "level")),
         width = 6.5, height = 4.5, dpi = 300)

  if (nrow(tiers$significant)) {
    fo <- poe_forest_nbmm(res = tiers$significant, ensembl = ensembl, tx2g = tx2g,
                          tx = lv == "transcript",
                          title = if (lv == "transcript") "Isoforms with significant POE" else "Genes with significant POE",
                          x_scale = if (lv == "transcript") "maternal_ratio" else "beta.ratio",
                          beta.thresh = log(1.5))
    ggsave(fig(paste0(pre, "_forest.png")), fo + theme(plot.subtitle = element_blank()),
           width = if (lv == "transcript") 4 else 5.75, height = 6.6, dpi = 300)
  }

  sel <- res[is.finite(res$q) & res$q < q_thresh, ]
  if (nrow(sel)) {
    cal <- plot_calibration_observed_vs_predicted(res = sel, poMat = as.data.frame(poMat), ensembl = ensembl,
                                                  tx2g = tx2g, level = lv, smooth_method = "lm")
    ggsave(fig(paste0(pre, "_calibration.png")), cal$plot, width = 6.5, height = 5.5, dpi = 300)
  }
}
level_plots(res_nb_tx,   tiers.tx,   tiers_all_tx,   poMat.tx,   "transcript", "isoPOE")
level_plots(res_nb_gene, tiers.gene, tiers_all_gene, poMat.gene, "gene",       "genePOE")

# ---- enrichment of known imprinted genes ----------------------------------------------
universe <- unique(sub("\\..*$", "", tx2g$gene_id))
imprintdb <- fread(P("IMPRINT_DB"))
imprint_genes <- unique(imprintdb$ensembl_gene_id)
tx_genes <- tx2g$gene_id[tx2g$transcript_id %in% tiers.tx$significant$feature_id]
skewed_genes <- unique(sub("\\..*$", "", c(tx_genes, tiers.gene$significant$feature_id)))

n11 <- length(intersect(skewed_genes, imprint_genes))
n10 <- length(setdiff(skewed_genes, imprint_genes))
n01 <- length(setdiff(imprint_genes, skewed_genes))
n00 <- length(setdiff(universe, union(skewed_genes, imprint_genes)))
cont_tab <- matrix(c(n11, n01, n10, n00), nrow = 2,
                   dimnames = list(POE = c("Yes", "No"), Imprinted = c("Yes", "No")))
ft <- fisher.test(cont_tab)
fwrite(data.table(poe_imprinted = n11, poe_not_imprinted = n10, imprinted_not_poe = n01, neither = n00,
                  odds_ratio = unname(ft$estimate), conf_lo = ft$conf.int[1], conf_hi = ft$conf.int[2],
                  p_value = ft$p.value),
       out("imprinting_enrichment.tsv"), sep = "\t")
fwrite(imprintdb[ensembl_gene_id %in% intersect(skewed_genes, imprint_genes)],
       out("imprinting_overlap_genes.tsv"), sep = "\t")
print(cont_tab); print(ft)
