# =============================================================================
# Stage 08 | combine_interaction.R
# Gather batch results; per level x condition, BH-adjust the interaction
# p-values across features and add derived columns.
# In : BATCH_DIR/results/interaction_<level>_NNN.rds (all batches required)
# Out: MODEL_DIR/poe_by_condition_<level>_<condition>.rds,
#      WORK_DIR/poe_interaction_summary.tsv
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages(library(data.table))
source(file.path(repo, "R", "poe_nbmm.R"))
mm <- readRDS(file.path(P("MODEL_DIR"), "model_metadata.rds"))
nb <- poe_n_batches()
summ <- list()
for (lv in names(nb)) {
  files <- vapply(seq_len(nb[[lv]]), function(k) poe_batch_result("interaction", lv, k), "")
  if (any(miss <- !file.exists(files)))
    stop(sprintf("[%s] %d of %d batch results missing, e.g. %s; rerun 08_fit",
                 lv, sum(miss), length(files), basename(files[miss][1])))
  all <- rbindlist(lapply(files, readRDS), fill = TRUE)
  for (cn in mm$conditions) {
    res <- poe_finish_interaction(all[condition == cn], mm$cinfo[[cn]], lv, cn, mm$opt)
    saveRDS(res, file.path(P("MODEL_DIR"), sprintf("poe_by_condition_%s_%s.rds", lv, cn)))
    summ[[length(summ) + 1]] <- data.table(level = lv, condition = cn, features = nrow(res),
      fitted = sum(res$boot_status == "ok", na.rm = TRUE),
      q_int_lt_0.10 = sum(res$q_int < 0.10, na.rm = TRUE),
      q_int_le_0.05 = sum(res$pass_diff_strict, na.rm = TRUE))
  }
}
summ <- rbindlist(summ)
fwrite(summ, out("poe_interaction_summary.tsv"), sep = "\t")
print(summ)
