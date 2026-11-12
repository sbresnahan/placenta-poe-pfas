# =============================================================================
# Stage 07 | combine_baseline.R
# Gather batch results, BH-adjust across features within each level.
# In : BATCH_DIR/results/baseline_<level>_NNN.rds (all batches required)
# Out: MODEL_DIR/poe_nbmm_results_{gene,transcript}.rds, WORK_DIR/poe_nbmm_results_*.tsv
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages(library(data.table))
source(file.path(repo, "R", "poe_nbmm.R"))
mm <- readRDS(file.path(P("MODEL_DIR"), "model_metadata.rds"))
nb <- poe_n_batches()
for (lv in names(nb)) {
  files <- vapply(seq_len(nb[[lv]]), function(k) poe_batch_result("baseline", lv, k), "")
  if (any(miss <- !file.exists(files)))
    stop(sprintf("[%s] %d of %d batch results missing, e.g. %s; rerun 07_fit",
                 lv, sum(miss), length(files), basename(files[miss][1])))
  res <- poe_finish_baseline(rbindlist(lapply(files, readRDS), fill = TRUE))
  attr(res, "size_factors") <- mm$size_factors[[lv]]
  saveRDS(res, file.path(P("MODEL_DIR"), sprintf("poe_nbmm_results_%s.rds", lv)))
  fwrite(res, out(sprintf("poe_nbmm_results_%s.tsv", lv)), sep = "\t")
  message(sprintf("[%s] %d features, %d fitted, %d with q < 0.1", lv, nrow(res),
                  sum(is.finite(res$p)), sum(res$q < 0.1, na.rm = TRUE)))
}
