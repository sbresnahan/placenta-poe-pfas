# =============================================================================
# R/setup.R -- sourced first by every R script in the pipeline
#
#   repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
#
# * imports every variable in config/paths.sh into the R environment, so R and
#   the job scripts share one path configuration (read values with P("NAME"))
# * adds R_LIB_DIR to the library path
# * NCPU = cores granted by LSF (1 when run interactively, unless set)
# =============================================================================
if (!exists("repo") || !nzchar(repo)) repo <- Sys.getenv("REPO", ".")
repo <- normalizePath(repo)
Sys.setenv(REPO = repo)

paths_init <- function(repo) {
  f <- file.path(repo, "config", "paths.sh")
  lines <- readLines(f)
  vars <- unique(unlist(regmatches(lines, gregexpr("(?<=\\$\\{)[A-Z0-9_]+(?=:=)", lines, perl = TRUE))))
  cmd <- sprintf("source '%s'; %s", f,
                 paste(sprintf("printf '%%s\\n' \"${%s}\"", vars), collapse = "; "))
  vals <- system2("bash", c("-c", shQuote(cmd)), stdout = TRUE)
  stopifnot(length(vals) == length(vars))
  do.call(Sys.setenv, as.list(setNames(vals, vars)))
  invisible(vars)
}
paths_init(repo)

P <- function(name) {
  v <- Sys.getenv(name)
  if (!nzchar(v)) stop("path variable ", name, " is empty; set it in config/paths.sh")
  v
}
words <- function(name) strsplit(trimws(Sys.getenv(name)), "[[:space:]]+")[[1]]

if (dir.exists(Sys.getenv("R_LIB_DIR"))) .libPaths(c(Sys.getenv("R_LIB_DIR"), .libPaths()))
NCPU <- max(1L, as.integer(Sys.getenv("NCPU", "1")))
for (d in c("FIG_DIR", "WORK_DIR")) dir.create(P(d), recursive = TRUE, showWarnings = FALSE)

if (!interactive()) grDevices::pdf(NULL)              # no stray Rplots.pdf from printed plots

fig <- function(name) file.path(P("FIG_DIR"), name)   # figure output path
out <- function(name) file.path(P("WORK_DIR"), name)  # table output path

# Attach the packages and functions used by the POE model drivers (stages 07-09)
load_poe_library <- function() {
  pk <- c("data.table", "ggplot2", "dplyr", "tidyr", "scales", "purrr", "tibble", "ggrepel", "grid")
  for (p in pk) suppressPackageStartupMessages(library(p, character.only = TRUE))
  for (f in c("utils", "poe_within_sample", "poe_nbmm", "poe_calibration", "poe_classify", "poe_plots"))
    source(file.path(repo, "R", paste0(f, ".R")), local = globalenv())
  invisible(TRUE)
}
