# =============================================================================
# Stage 10 | aggregate_PFAS_sim_results.R
# Type-I error of each mediation method across the simulated studies, and
# the DAG of the planned mediation analysis.
# In : SIM_DIR/run_NNN.RData (SIM_N_RUNS runs; all required)
# Out: WORK_DIR/sim_type1_error.tsv, FIG_DIR/sim_type1_error.png (DAG.R: FIG_DIR/PFAS_POE_DAG.png)
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

n_runs <- as.integer(P("SIM_N_RUNS"))
files <- file.path(P("SIM_DIR"), sprintf("run_%03d.RData", seq_len(n_runs)))
if (any(miss <- !file.exists(files)))
  stop(sum(miss), " of ", n_runs, " simulation runs missing, e.g. ", basename(files[miss][1]), "; rerun 10_sim")
t1 <- rbindlist(lapply(files, function(f) { e <- new.env(); load(f, envir = e); as.data.table(e$combined_type1_errors) }),
                fill = TRUE)
t1[, label := paste0(Method, ": ", path)]
fwrite(t1, out("sim_type1_error.tsv"), sep = "\t")

g <- ggplot(t1, aes(label, Type1_Error_Rate, color = Method)) +
  geom_boxplot(fill = NA, outliers = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_color_viridis_d(end = 0.85) +
  coord_flip() +
  labs(x = NULL, y = "Type I error rate", title = sprintf("Type I error across %d simulated studies", n_runs)) +
  theme_minimal()
ggsave(fig("sim_type1_error.png"), g, width = 8.5, height = 4.5, dpi = 300)
message("aggregated ", n_runs, " runs")
