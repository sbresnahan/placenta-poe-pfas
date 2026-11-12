# =============================================================================
# Stage 08 | fit_interaction.R    one array element = one feature batch,
# all conditions in POE_CONDITIONS
# Interaction NB-GLMM per feature x condition X:
#   count ~ allele * X + ancestry + sex + (1|SNPpos) + (1|sample_id), offset log(SF)
# Samples with missing X are dropped for that condition.
# In : BATCH_DIR/<level>_batch_NNN.rds, MODEL_DIR/model_metadata.rds
# Out: BATCH_DIR/results/interaction_<level>_NNN.rds  (column `condition`)
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({ library(data.table); library(parallel) })
source(file.path(repo, "R", "poe_nbmm.R"))

b  <- poe_batch_of(as.integer(Sys.getenv("IDX", "1")))
mm <- readRDS(file.path(P("MODEL_DIR"), "model_metadata.rds"))
conds <- mm$conditions
ccols <- paste0(".c_", conds)
dt <- readRDS(poe_batch_input(b$level, b$k))
dt <- merge(dt, mm$meta[, c("sample_id", "sex", "ancestry", ccols), with = FALSE], by = "sample_id")
feats <- split(dt, by = "feature_id", keep.by = TRUE)
message(sprintf("[%s batch %d] %d features x %d conditions on %d cores",
                b$level, b$k, length(feats), length(conds), NCPU))

fit_feature <- function(sub) {
  set.seed(poe_seed(sub$feature_id[1]))
  rbindlist(lapply(conds, function(cn) {
    s <- sub[!is.na(get(paste0(".c_", cn)))]
    s[, .cond := get(paste0(".c_", cn))]
    s[, (ccols) := NULL]
    r <- fit_poe_interaction_one(s, mm$cinfo[[cn]], mm$opt)
    r$feature_id <- sub$feature_id[1]          # also defined when s is empty
    r$condition <- cn
    r
  }), fill = TRUE)
}

t0 <- Sys.time()
res <- mclapply(feats, fit_feature, mc.cores = NCPU, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop(sum(bad), " features failed; first error: ", res[[which(bad)[1]]])
res <- rbindlist(res, fill = TRUE)
res[, level := b$level]
saveRDS(res, poe_batch_result("interaction", b$level, b$k))
message(sprintf("done: %d fits (%d ok) in %.1f min", nrow(res), sum(res$boot_status == "ok", na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
