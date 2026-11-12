# =============================================================================
# Stage 03 | imputation_qc.R
# Summaries of the SEAD imputation: sites retained by the MAF filter, merged
# variant density per kb of each autosome, imputed MAF / R2 distributions.
# In : WORK_DIR/imputation_site_counts.tsv, WORK_DIR/imputation_maf.tsv, HG38_FASTA.fai
# Out: WORK_DIR/imputation_summary.tsv,
#      FIG_DIR/imputation_density_per_kb.png, FIG_DIR/imputation_maf_r2.png
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

counts <- fread(out("imputation_site_counts.tsv"), na.strings = "NA")
fai <- fread(paste0(P("HG38_FASTA"), ".fai"), header = FALSE, select = 1:2, col.names = c("seq", "len"))
counts[, len := fai$len[match(paste0("chr", chr), fai$seq)]]
counts[, `:=`(frac_kept = filtered / dose, per_kb = filtered / (len / 1000))]
fwrite(counts, out("imputation_summary.tsv"), sep = "\t")

tot <- counts[, .(dose = sum(dose, na.rm = TRUE), filtered = sum(filtered, na.rm = TRUE)), by = group]
print(tot)

g <- ggplot(counts, aes(factor(chr, levels = unique(chr)), per_kb, fill = group)) +
  geom_col(position = "dodge") +
  labs(x = "Chromosome", y = "Variants per kb (MAF > 0.01)", fill = "Group") +
  theme_minimal()
ggsave(fig("imputation_density_per_kb.png"), g, width = 9, height = 4, dpi = 300, bg = "white")

maf_file <- out("imputation_maf.tsv")
if (file.exists(maf_file) && length(readLines(maf_file, n = 2)) > 1) {
  maf <- fread(maf_file, na.strings = c(".", "NA"))
  maf_long <- melt(maf, id.vars = c("group", "chr"), measure.vars = c("MAF", "R2"))
  print(maf[, .(n = .N, frac_MAF_gt_0.01 = mean(MAF > 0.01, na.rm = TRUE),
                median_R2 = median(R2, na.rm = TRUE)), by = group])
  g <- ggplot(maf_long[!is.na(value)], aes(value, fill = group)) +
    geom_histogram(bins = 50, position = "identity", alpha = 0.5) +
    facet_wrap(~variable, scales = "free") +
    labs(x = NULL, y = "Sites", fill = "Group") + theme_minimal()
  ggsave(fig("imputation_maf_r2.png"), g, width = 8, height = 4, dpi = 300, bg = "white")
} else {
  message("no .info files found; MAF/R2 summary skipped")
}
