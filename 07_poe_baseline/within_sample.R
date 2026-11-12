# =============================================================================
# Stage 07 | within_sample.R
# Within-sample parent-of-origin tests: per-SNP exact binomial tests of
# maternal vs paternal counts, ACAT-combined per feature within each sample
# (count weights), BH across features within the sample. Also summarizes
# per-feature / per-sample read support.
#
# In : MODEL_INPUTS (stage 07 model_inputs.R), TX_CLASSIF
# Out: MODEL_DIR/within_sample_tests.rds
#      WORK_DIR/within_sample_biased_genes_per_sample.tsv
#      WORK_DIR/feature_count_summary.tsv
#      FIG_DIR/significant_POE_totals.png, significant_POE_props.png,
#              per_sample_POE_pdists.png
# Biased = POE (combined q < 0.05) and summary maternal ratio < 0.1 or > 0.9.
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
load_poe_library()

inp  <- readRDS(P("MODEL_INPUTS"))
tx2g <- read_tx2g()
drop <- excluded_samples()
keep_cols <- function(m) m[, !sub("_(maternal|paternal)$", "", colnames(m)) %in% drop, drop = FALSE]
poMat.tx   <- keep_cols(inp$poMat.tx)
poMat.gene <- keep_cols(inp$poMat.gene)

res.tx   <- poe_tx_tests_hierarchical(as.data.frame(poMat.tx), rounding = "round", q_threshold = 0.05,
                                      acat_weights = "counts", n_workers = NCPU, verbose = FALSE)
res.gene <- poe_tx_tests_hierarchical(as.data.frame(poMat.gene), rounding = "round", q_threshold = 0.05,
                                      acat_weights = "counts", n_workers = NCPU, verbose = FALSE)
saveRDS(list(transcript = res.tx, gene = res.gene), file.path(P("MODEL_DIR"), "within_sample_tests.rds"))

# ---- biased genes per sample: transcript-derived, gene-level, or both ------------
hi <- 0.9; lo <- 0.1
sig <- function(s) dplyr::filter(s, POE, !is.na(ratio_summary), ratio_summary < lo | ratio_summary > hi)
sig_tx_genes <- sig(res.tx$summary) %>%
  select(sample_id, transcript_id = feature_id) %>%
  inner_join(as.data.frame(tx2g), by = "transcript_id") %>%
  distinct(sample_id, gene_id)
sig_gene_genes <- sig(res.gene$summary) %>%
  transmute(sample_id, gene_id = feature_id) %>%
  distinct(sample_id, gene_id)

counts <- bind_rows(
  anti_join(sig_tx_genes, sig_gene_genes, by = c("sample_id", "gene_id")) %>% mutate(category = "transcript"),
  anti_join(sig_gene_genes, sig_tx_genes, by = c("sample_id", "gene_id")) %>% mutate(category = "gene"),
  inner_join(sig_tx_genes, sig_gene_genes, by = c("sample_id", "gene_id")) %>% mutate(category = "both")
) %>% count(sample_id, category, name = "N") %>%
  complete(sample_id = union(sig_tx_genes$sample_id, sig_gene_genes$sample_id),
           category = c("transcript", "gene", "both"), fill = list(N = 0L))
fwrite(counts, out("within_sample_biased_genes_per_sample.tsv"), sep = "\t")

counts_ordered <- counts %>% group_by(sample_id) %>% mutate(total = sum(N)) %>% ungroup() %>%
  mutate(sample_id = reorder(sample_id, -total),
         category = factor(category, levels = c("transcript", "gene", "both")))
median_total <- median(distinct(counts_ordered, sample_id, total)$total)

g_tot <- ggplot(counts_ordered, aes(sample_id, N, fill = category)) +
  geom_col(width = 1, color = NA) +
  geom_hline(yintercept = median_total, linetype = "dashed") +
  scale_fill_viridis_d(end = 0.85) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Sample", y = "Genes", fill = NULL, title = "Genes with biased haplotype expression") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(fig("significant_POE_totals.png"), g_tot, width = 8.2, height = 3.5, dpi = 300)

g_prop <- ggplot(counts_ordered, aes(sample_id, N, fill = category)) +
  geom_col(width = 1, color = NA, position = "fill") +
  scale_fill_viridis_d(end = 0.85) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(x = "Sample", y = "Proportion of biased genes", fill = NULL, title = "Level of detection") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(fig("significant_POE_props.png"), g_prop, width = 7.5, height = 3.5, dpi = 300)

p_lists <- list("Transcript, per SNP (binomial)" = res.tx$tests$p_value,
                "Transcript, per feature (ACAT)" = res.tx$summary$combined_p,
                "Gene, per SNP (binomial)"       = res.gene$tests$p_value,
                "Gene, per feature (ACAT)"       = res.gene$summary$combined_p)
p_long <- imap_dfr(p_lists, ~ tibble(p = .x, panel = .y)) %>%
  filter(is.finite(p)) %>%
  mutate(p = pmax(p, .Machine$double.xmin), panel = factor(panel, levels = names(p_lists)))
g_p <- ggplot(p_long, aes(p)) +
  geom_histogram(bins = 50, boundary = 0, closed = "left") +
  facet_wrap(~ panel, nrow = 2, scales = "free_y") +
  labs(x = "p-value", y = "Count", title = "Within-sample p-value distributions") +
  theme_minimal(base_size = 12)
ggsave(fig("per_sample_POE_pdists.png"), g_p, width = 8, height = 5.5, dpi = 300)

# ---- read support per feature and per sample ---------------------------------
support <- function(m, level) {
  fid <- sub(":[^:]*$", "", rownames(m))
  feat <- rowsum(m, fid, reorder = FALSE)                       # feature x (sample_allele)
  samp <- sub("_(maternal|paternal)$", "", colnames(feat))
  per_sample <- t(rowsum(t(feat), samp, reorder = FALSE))       # feature x sample
  ft <- rowSums(feat)
  sm <- apply(per_sample, 2, function(v) if (any(v > 0)) mean(v[v > 0]) else NA_real_)
  q <- function(v) as.list(stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, names = FALSE))
  rbind(data.table(level = level, quantity = "total reads per feature (features with reads)",
                   n = sum(ft > 0), mean = mean(ft[ft > 0]))[, c("min", "q25", "median", "q75", "max") := q(ft[ft > 0])],
        data.table(level = level, quantity = "mean reads per expressed feature, per sample",
                   n = sum(is.finite(sm)), mean = mean(sm, na.rm = TRUE))[, c("min", "q25", "median", "q75", "max") := q(sm)])
}
fwrite(rbind(support(poMat.gene, "gene"), support(poMat.tx, "transcript")),
       out("feature_count_summary.tsv"), sep = "\t")
message("within-sample tests done")
