# =============================================================================
# Stage 09 | power_analysis.R
# Power of the cohort-level POE test as a function of sample size, simulated
# from the observed standard errors of the baseline NB-GLMMs.
#   * SEs are resampled from the fitted features and scaled by sqrt(n0 / N)
#     (n0 = largest number of samples contributing to a feature)
#   * a fraction pi1 of the m tested features carries a true effect beta
#   * discoveries: BH q < 0.10 and |beta_hat| >= log(1.5)
#   * power = true discoveries / m1, averaged over 1000 simulated studies
#
# In : MODEL_DIR/poe_nbmm_results_{transcript,gene}.rds
# Out: WORK_DIR/power_table_<level>.tsv
#      FIG_DIR/power_curves_<level>_pi1_<pct>.png
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(scales); library(parallel) })

q_bh          <- 0.10
beta_cut      <- log(1.5)
effect_ratios <- c(1.5, 1.6, 1.7, 1.8)
pi1_values    <- c(0.005, 0.01, 0.02)
n_sims        <- 1000

simulate_power_poe <- function(
    se_vec, m = length(se_vec), pi1 = 0.01,
    beta_true_grid = log(1.5),
    q = 0.10, abs_beta_min = log(1.5),
    n_sims = 1000, se_scale = 1.0, two_sided = TRUE, seed = 1,
    use_effect_cutoff = TRUE
){
  stopifnot(is.numeric(se_vec), length(se_vec) > 100)
  se_vec <- se_vec[is.finite(se_vec) & se_vec > 0]
  set.seed(seed)
  draw_se <- function(k) sample(se_vec, size = k, replace = TRUE) * se_scale
  
  one_beta <- function(beta_true){
    pow <- fdr <- disc <- numeric(n_sims)
    m1 <- round(m * pi1); m0 <- m - m1
    for (s in seq_len(n_sims)) {
      se0 <- draw_se(m0); se1 <- draw_se(m1)
      bh0 <- rnorm(m0, 0,         se0)         # beta_hat under H0
      bh1 <- rnorm(m1, beta_true, se1)         # beta_hat under H1
      z0  <- bh0/se0; z1 <- bh1/se1
      p0  <- if (two_sided) 2*pnorm(-abs(z0)) else 1 - pnorm(z0)
      p1  <- if (two_sided) 2*pnorm(-abs(z1)) else 1 - pnorm(z1)
      p_all <- c(p0, p1)
      
      # BH step-up
      ord <- order(p_all); ps <- p_all[ord]; thr <- (seq_along(ps)/length(ps)) * q
      k   <- max(which(ps <= thr), na.rm = TRUE); if (!is.finite(k)) k <- 0L
      call_bh <- rep(FALSE, length(p_all)); if (k > 0) call_bh[ord[seq_len(k)]] <- TRUE
      
      # Effect-size gate
      if (use_effect_cutoff) {
        call_mag <- c(abs(bh0) >= abs_beta_min, abs(bh1) >= abs_beta_min)
        call <- call_bh & call_mag
      } else {
        call <- call_bh
      }
      
      tp <- sum(call[(m0+1):(m0+m1)])
      fp <- sum(call[1:m0])
      R  <- sum(call)
      pow[s] <- ifelse(m1 > 0, tp / m1, NA_real_)
      fdr[s] <- ifelse(R > 0, fp / R, 0)
      disc[s] <- R
    }
    data.frame(beta_true = beta_true,
               power = mean(pow, na.rm = TRUE),
               mean_FDR = mean(fdr, na.rm = TRUE),
               mean_discoveries = mean(disc, na.rm = TRUE))
  }
  
  do.call(rbind, lapply(beta_true_grid, one_beta))
}




for (lv in c("transcript", "gene")) {
  res <- readRDS(file.path(P("MODEL_DIR"), sprintf("poe_nbmm_results_%s.rds", lv)))
  stopifnot(all(c("se", "p", "total_counts", "n_samples", "beta") %in% names(res)))
  n0     <- max(res$n_samples, na.rm = TRUE)
  m_tot  <- sum(is.finite(res$p))                    # features tested
  se_vec <- res$se[is.finite(res$se) & res$se > 0 & is.finite(res$total_counts)]
  N_grid <- seq(max(30, floor(0.5 * n0)), max(7 * n0, n0 + 5), by = 100)
  grid <- CJ(pi1 = pi1_values, N = N_grid)
  message(sprintf("[%s] n0 = %d, m = %d, %d scenarios on %d cores", lv, n0, m_tot, nrow(grid), NCPU))
  pw <- rbindlist(mclapply(seq_len(nrow(grid)), function(i) {
    ans <- simulate_power_poe(se_vec = se_vec, m = m_tot, pi1 = grid$pi1[i],
                              beta_true_grid = log(effect_ratios), q = q_bh,
                              abs_beta_min = beta_cut, n_sims = n_sims,
                              se_scale = sqrt(n0 / grid$N[i]), seed = i)
    data.table(level = lv, pi1 = grid$pi1[i], N = grid$N[i], ans)
  }, mc.cores = NCPU, mc.preschedule = FALSE))
  pw[, `:=`(effect_ratio = exp(beta_true), m = m_tot, n0 = n0)]
  fwrite(pw[order(pi1, effect_ratio, N)], out(sprintf("power_table_%s.tsv", lv)), sep = "\t")

  for (p1 in pi1_values) {
    g <- ggplot(pw[pi1 == p1], aes(N, power, color = factor(round(effect_ratio, 2)))) +
      geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
      geom_hline(yintercept = 0.80, linetype = 3) +
      scale_color_viridis_d(end = 0.85) +
      scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, 1)) +
      labs(title = sprintf("POE power, %s level (pi1 = %s, m = %d)", lv, percent(p1, accuracy = 0.1), m_tot),
           x = "Sample size (N)", y = "Power", color = "Effect ratio") +
      theme_minimal(base_size = 12)
    ggsave(fig(sprintf("power_curves_%s_pi1_%s.png", lv, sub("\\.", "p", format(100 * p1)))), g,
           width = 8, height = 4.8, dpi = 300)
  }
}
message("power analysis done")
