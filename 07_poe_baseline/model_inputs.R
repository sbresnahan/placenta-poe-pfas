# =============================================================================
# Stage 07 | model_inputs.R
# Count matrices + sample metadata -> model inputs and per-batch long tables
# for the baseline (stage 07) and interaction (stage 08) NB-GLMMs.
#
# In : POE_GENE_COUNTS, POE_TX_COUNTS (stage 06), RNASEQ_COVARS_CSV,
#      PFAS_MAT_CSV, PFAS_CHILD_CSV, optional EXCLUDE_SAMPLES_FILE
# Out: MODEL_INPUTS                 list(poMat.gene, poMat.tx, metadata)
#      MODEL_DIR/model_metadata.rds modelled samples, covariates, coded
#                                   conditions, condition info, size factors
#      BATCH_DIR/<level>_batch_NNN.rds  long tables, one per array element
#      WORK_DIR/model_inputs_summary.tsv
#
# Sample metadata
#   RNASEQ_COVARS_CSV columns SubjectID, sex, mother_ethnicity (-> ancestry),
#   GA, gdm (-> GDM); PFAS tables are joined by SubjectID (first 9 columns:
#   SubjectID + 8 PFAS with _mat / _child suffixes). Sample IDs follow the
#   genotype VCF convention <SubjectID>_B<SubjectID> with "-" -> "_".
#   PFAS values are used on the scale supplied in the CSVs.
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages(library(data.table))
source(file.path(repo, "R", "utils.R")); source(file.path(repo, "R", "poe_nbmm.R"))
opt <- POE_FIT_DEFAULTS
opt$n_boot <- as.integer(Sys.getenv("N_BOOT", opt$n_boot))
conds <- words("POE_CONDITIONS")

# ---- metadata ----------------------------------------------------------------
cov <- fread(P("RNASEQ_COVARS_CSV"), na.strings = c("", "NA"))
need <- c("SubjectID", "sex", "mother_ethnicity", "GA", "gdm")
if (length(miss <- setdiff(need, names(cov)))) stop("RNASEQ_COVARS_CSV lacks: ", paste(miss, collapse = ", "))
meta <- cov[, ..need]
for (f in c("PFAS_MAT_CSV", "PFAS_CHILD_CSV")) {
  pf <- fread(P(f), na.strings = c("", "NA"))[, 1:9]
  if (names(pf)[1] != "SubjectID") stop(f, ": first column must be SubjectID")
  meta <- merge(meta, pf, by = "SubjectID", all.x = TRUE, sort = FALSE)
}
setnames(meta, c("SubjectID", "sex", "mother_ethnicity", "GA", "gdm"),
         c("sample_id", "sex", "ancestry", "GA", "GDM"))
meta[, sample_id := gsub("-", "_", sample_id)]
meta[, sample_id := paste0(sample_id, "_B", sample_id)]
if (anyDuplicated(meta$sample_id)) stop("duplicated SubjectID in RNASEQ_COVARS_CSV")
if (length(miss <- setdiff(conds, names(meta))))
  stop("POE_CONDITIONS not found in the metadata: ", paste(miss, collapse = ", "))

# Condition coding uses the whole RNA-seq cohort (step = SD for continuous
# exposures; 0/1 for binary, reference = first level)
coded <- lapply(setNames(conds, conds), function(cn) poe_code_condition(meta[[cn]]))

# ---- count matrices ------------------------------------------------------------
poMat.gene <- read_poe_counts(P("POE_GENE_COUNTS"), opt$min_snps)
poMat.tx   <- read_poe_counts(P("POE_TX_COUNTS"),   opt$min_snps)

drop <- excluded_samples()
mat_samples <- unique(sub("_(maternal|paternal)$", "", colnames(poMat.gene)))
model_samples <- meta[sample_id %in% mat_samples & !sample_id %in% drop &
                        !is.na(sex) & !is.na(ancestry), sample_id]
message(sprintf("samples: %d in counts, %d with metadata, %d excluded by file, %d modelled",
                length(mat_samples), sum(mat_samples %in% meta$sample_id),
                sum(mat_samples %in% drop), length(model_samples)))
if (length(model_samples) < opt$min_samples) stop("too few modelled samples")
if (length(no_meta <- setdiff(mat_samples, meta$sample_id)))
  message("count-matrix samples without metadata (dropped): ", length(no_meta))

meta_m <- meta[match(model_samples, sample_id)]
cinfo <- lapply(setNames(conds, conds), function(cn) {
  cc <- coded[[cn]]
  x <- cc$x[match(model_samples, meta$sample_id)]
  list(type = cc$type, step = cc$step, levels = cc$levels,
       A = mean(x, na.rm = TRUE), B = stats::sd(x, na.rm = TRUE),
       n_samples = sum(!is.na(x)))
})
for (cn in conds) meta_m[, (paste0(".c_", cn)) := coded[[cn]]$x[match(model_samples, meta$sample_id)]]

# ---- per-level size factors, feature filter and batches -------------------------
dir.create(file.path(P("BATCH_DIR"), "results"), recursive = TRUE, showWarnings = FALSE)
nb <- poe_n_batches()
sfs <- list(); summ <- list()
for (lv in c("gene", "transcript")) {
  mat <- if (lv == "gene") poMat.gene else poMat.tx
  sp <- poe_split_matrix(mat)
  keep_s <- intersect(model_samples, sp$samples)
  M <- sp$M[, keep_s, drop = FALSE]; Pm <- sp$P[, keep_s, drop = FALSE]
  sf <- poe_size_factors(M, Pm); sfs[[lv]] <- sf
  keep_f <- poe_keep_features(M, Pm, sp$feature_id, opt$min_snps, opt$min_samples, opt$min_total_counts)
  rows_per_f <- table(sp$feature_id[sp$feature_id %in% keep_f])
  bat <- poe_assign_batches(setNames(as.numeric(rows_per_f), names(rows_per_f)), nb[[lv]])
  row_batch <- bat[sp$feature_id]
  for (k in seq_len(nb[[lv]])) {
    rows <- which(row_batch == k)
    saveRDS(poe_long(M, Pm, sp$feature_id, sp$snppos, sf, rows), poe_batch_input(lv, k), compress = FALSE)
  }
  summ[[lv]] <- data.table(level = lv, rows = nrow(mat), features = length(unique(sp$feature_id)),
                           features_modelled = length(keep_f), samples = length(keep_s),
                           batches = nb[[lv]])
  message(sprintf("[%s] %d features pass filters -> %d batches", lv, length(keep_f), nb[[lv]]))
}

dir.create(P("MODEL_DIR"), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(poMat.gene = poMat.gene, poMat.tx = poMat.tx, metadata = as.data.frame(meta)),
        P("MODEL_INPUTS"))
saveRDS(list(meta = meta_m, conditions = conds, cinfo = cinfo, size_factors = sfs,
             opt = opt, n_batches = nb),
        file.path(P("MODEL_DIR"), "model_metadata.rds"))

cond_tab <- rbindlist(lapply(conds, function(cn) data.table(
  condition = cn, type = cinfo[[cn]]$type, step = cinfo[[cn]]$step,
  modelled_mean = cinfo[[cn]]$A, modelled_sd = cinfo[[cn]]$B,
  n_samples = cinfo[[cn]]$n_samples,
  levels = paste(cinfo[[cn]]$levels, collapse = "/"))))
fwrite(rbindlist(summ), out("model_inputs_summary.tsv"), sep = "\t")
fwrite(cond_tab, out("model_conditions.tsv"), sep = "\t")
print(rbindlist(summ)); print(cond_tab)
