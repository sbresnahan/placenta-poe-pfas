# =============================================================================
# reports/render_reports.R -- re-knit the figure-based reports into docs/reports/.
#
#   Rscript reports/render_reports.R            # both reports
#   Rscript reports/render_reports.R methods    # Placenta_POE_PFAS_methods only
#   Rscript reports/render_reports.R all_PFAS   # all_PFAS_report only
#
# Figures are read from reports/figures/{methods,interaction}/ (stage them with
# reports/stage_figures.sh). To read stage 08 output directly instead, set
#   REPORT_INTERACTION_DIR=/path/to/FIG_DIR/interaction
# Needs: rmarkdown, knitr, cowplot, magick, ggplot2 (and pandoc).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
which <- if (length(args)) args else c("methods", "all_PFAS")

f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
rep_dir <- normalizePath(if (length(f)) dirname(f) else "reports")
out_dir <- file.path(dirname(rep_dir), "docs", "reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("rmarkdown", "knitr", "cowplot", "magick", "ggplot2"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("install R package '", pkg, "'", call. = FALSE)

if ("methods" %in% which)
  rmarkdown::render(file.path(rep_dir, "Placenta_POE_PFAS_methods.Rmd"),
                    output_format = "html_document", output_dir = out_dir,
                    params = list(fig_dir = file.path(rep_dir, "figures", "methods")),
                    envir = new.env())

if ("all_PFAS" %in% which) {
  idir <- Sys.getenv("REPORT_INTERACTION_DIR", file.path(rep_dir, "figures", "interaction"))
  rmarkdown::render(file.path(rep_dir, "all_PFAS_report.Rmd"),
                    output_format = "html_document", output_dir = out_dir,
                    params = list(fig_dir = normalizePath(idir, mustWork = FALSE)),
                    envir = new.env())
}
message("rendered into ", out_dir)
