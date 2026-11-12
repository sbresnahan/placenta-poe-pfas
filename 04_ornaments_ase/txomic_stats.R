# =============================================================================
# Stage 04 | txomic_stats.R
# Distribution of imputed variants per transcript (RNA-seq infants).
# In : WORK_DIR/txomic_variants_per_tx.tsv, TX_CLASSIF
# Out: WORK_DIR/txomic_variants_summary.tsv, FIG_DIR/txomic_variants_per_tx.png
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

tab <- fread(out("txomic_variants_per_tx.tsv"))
tab <- tab[, .(n_variants = sum(n_variants)), by = transcript]   # robust to unsorted input
n_tx <- nrow(fread(P("TX_CLASSIF"), select = 1))

summ <- data.table(
  transcripts_total = n_tx,
  transcripts_with_variant = nrow(tab),
  pct_with_variant = round(100 * nrow(tab) / n_tx, 2),
  transcripts_gt1_variant = tab[n_variants > 1, .N],
  pct_gt1_variant = round(100 * tab[n_variants > 1, .N] / n_tx, 2),
  median_variants = median(tab$n_variants),
  max_variants = max(tab$n_variants))
print(summ)
fwrite(summ, out("txomic_variants_summary.tsv"), sep = "\t")

g <- ggplot(tab, aes(n_variants)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  labs(x = "Variants per transcript", y = "Transcripts") +
  theme_minimal()
ggsave(fig("txomic_variants_per_tx.png"), g, width = 6, height = 5, dpi = 300, bg = "white")
