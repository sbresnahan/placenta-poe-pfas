# =============================================================================
# Stage 07 | fit_baseline.R      one array element = one feature batch
# Baseline NB-GLMM per feature:
#   count ~ allele + sex + ancestry + (1|SNPpos) + (1|sample_id), offset log(SF)
# In : BATCH_DIR/<level>_batch_NNN.rds, MODEL_DIR/model_metadata.rds
# Out: BATCH_DIR/results/baseline_<level>_NNN.rds
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({ library(data.table); library(parallel) })
source(file.path(repo, "R", "poe_nbmm.R"))

b  <- poe_batch_of(as.integer(Sys.getenv("IDX", "1")))
mm <- readRDS(file.path(P("MODEL_DIR"), "model_metadata.rds"))
dt <- readRDS(poe_batch_input(b$level, b$k))
dt <- merge(dt, mm$meta[, .(sample_id, sex, ancestry)], by = "sample_id")
feats <- split(dt, by = "feature_id", keep.by = TRUE)
message(sprintf("[%s batch %d] %d features on %d cores", b$level, b$k, length(feats), NCPU))

t0 <- Sys.time()
res <- mclapply(feats, fit_poe_baseline_one, opt = mm$opt,
                mc.cores = NCPU, mc.preschedule = FALSE)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop(sum(bad), " features failed; first error: ", res[[which(bad)[1]]])
res <- rbindlist(res, fill = TRUE)
res[, level := b$level]
saveRDS(res, poe_batch_result("baseline", b$level, b$k))
message(sprintf("done: %d fits (%d with estimates) in %.1f min", nrow(res), sum(is.finite(res$beta)),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
