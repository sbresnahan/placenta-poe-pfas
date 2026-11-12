# =============================================================================
# Stage 08 | dose_response.R
# Dose-response curves of the parental allelic ratio for interaction hits
# (tier "any") of the conditions in DOSE_RESPONSE_CONDITIONS.
#   cohort-centred: x = exposure on the CSV scale, annotated with the whole
#                   GUSTO cohort mean / SD (PFAS_*_COHORT_XLSX)
#   subset-centred: x = SDs from the mean of the modelled RNA-seq samples
#
# In : MODEL_DIR/poe_diff_tiers.rds, MODEL_DIR/model_metadata.rds,
#      PFAS_MAT_COHORT_XLSX, PFAS_CORD_COHORT_XLSX, TX_CLASSIF, PLACENTA_DB
# Out: FIG_DIR/dose_response/<cond>_Dose_Response_{Gene_cohort,Gene_subset,Transcript_subset}.png
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
load_poe_library()
suppressPackageStartupMessages(library(readxl))

mm <- readRDS(file.path(P("MODEL_DIR"), "model_metadata.rds"))
tiers <- readRDS(file.path(P("MODEL_DIR"), "poe_diff_tiers.rds"))
tx2g <- as.data.frame(read_tx2g()); ensembl <- read_ensembl()
meta <- as.data.frame(mm$meta)
dose_conds <- words("DOSE_RESPONSE_CONDITIONS")
exclude_symbols <- strsplit(Sys.getenv("DOSE_RESPONSE_EXCLUDE"), "[[:space:]]+")[[1]]
fdir <- fig("dose_response"); dir.create(fdir, showWarnings = FALSE)

read_cohort <- function(f, suffix) {
  x <- as.data.frame(read_excel(P(f)))[, 1:9]
  names(x)[2:9] <- paste0(names(x)[2:9], suffix)
  x$SubjectID <- gsub("-", "_", x$SubjectID)
  x
}
cohort <- merge(read_cohort("PFAS_CORD_COHORT_XLSX", "_child"),
                read_cohort("PFAS_MAT_COHORT_XLSX", "_mat"), by = "SubjectID", all = TRUE)

for (cn in dose_conds) {
  if (!cn %in% mm$conditions) stop(cn, " is not in POE_CONDITIONS")
  if (mm$cinfo[[cn]]$type != "continuous") { message(cn, ": binary condition, skipped"); next }
  label <- sprintf("%s %s", if (grepl("_mat$", cn)) "Maternal" else "Child", sub("_(mat|child)$", "", cn))
  runs <- list(
    Gene_cohort       = list(lv = "gene", centering = "cohort", cohort = cohort, ex = exclude_symbols),
    Gene_subset       = list(lv = "gene", centering = "subset", cohort = NULL, ex = exclude_symbols),
    Transcript_subset = list(lv = "transcript", centering = "subset", cohort = NULL, ex = NULL))
  for (nm in names(runs)) {
    r <- runs[[nm]]
    p <- plot_dose_response(res = tiers[[r$lv]][[cn]], filter_tier = "any", condition_name = cn,
                            tx2g = tx2g, ensembl = ensembl, level = r$lv, metadata = meta,
                            plot_title = sprintf("%s-level dose-response of allelic ratios (%s)",
                                                 tools::toTitleCase(r$lv), label),
                            centering = r$centering, cohort_pfas = r$cohort,
                            exclude = if (length(r$ex) && any(nzchar(r$ex))) r$ex else NULL)
    if (is.null(p)) { message(cn, " ", nm, ": no features to plot"); next }
    ggsave(file.path(fdir, sprintf("%s_Dose_Response_%s.png", cn, nm)), p,
           width = 10, height = 6, dpi = 300)
  }
}
message("dose-response plots done")
