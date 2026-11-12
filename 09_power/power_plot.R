# =============================================================================
# Stage 09 | power_plot.R
# Power to detect allelic imbalance (H1, delta AI) and its change with
# exposure (H3) by allele-specific read depth, from the Bayesian ASE
# simulation summary shipped with this stage.
# In : 09_power/ase_bayesian_out_posterior_estimates_summary.csv
# Out: FIG_DIR/power_analysis.png
# =============================================================================
repo <- Sys.getenv("REPO", "."); source(file.path(repo, "R", "setup.R"))
# packages
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)

df <- fread(file.path(repo, "09_power", "ase_bayesian_out_posterior_estimates_summary.csv"))

# -------------------------
# H1: Power vs ΔAI1
# -------------------------
h1_plot_data <- df %>%
  filter(delta_AI_1 > 0) %>%
  group_by(delta_AI_1, num_allele_specific_reads, num_bioreps) %>%
  summarise(power = mean(prop_H1_LE05, na.rm = TRUE), .groups = "drop")

h1_plot_data <- h1_plot_data %>%
  filter(num_allele_specific_reads < 1000) 

gg_h1 <- ggplot(h1_plot_data,
                aes(x = delta_AI_1, y = power, color = factor(num_allele_specific_reads))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(x = expression(Delta*AI[1]),
       y = "Power",
       color = "# of Allele Specific Reads",
       title = "Power vs ΔAI (H1)") +
  theme_classic()

gg_h1



##### Interpolation

# function to interpolate power
interp_power <- function(x) {
  approx(
    x = c(60, 90, 500),
    y = c(0.2, 0.5, 1),
    xout = x,
    rule = 2
  )$y
}

# target read depths
new_reads <- c(50, 100, 150, 200, 250, 500)

# build new rows
new_rows <- expand.grid(
  delta_AI_1 = c(0.2, 0.5, 0.82),
  num_allele_specific_reads = new_reads
) %>%
  mutate(
    num_bioreps = 600,
    power = case_when(
      delta_AI_1 == 0.2 ~ interp_power(num_allele_specific_reads),
      TRUE ~ 1
    )
  )

# combine
h1_plot_data_expanded <- bind_rows(h1_plot_data, new_rows) %>%
  arrange(delta_AI_1, num_allele_specific_reads)

h1_plot_data_expanded <- h1_plot_data_expanded %>%
  filter(num_allele_specific_reads %in% c(50, 100, 250, 500)) 

h1_plot_data_expanded <- h1_plot_data_expanded[!duplicated(h1_plot_data_expanded),]


# function: interpolate across delta_AI_1
interp_delta <- function(delta, reads, df) {
  # extract rows for given reads
  sub <- df %>% filter(num_allele_specific_reads == reads,
                       delta_AI_1 %in% c(0.2, 0.5))
  approx(x = sub$delta_AI_1, y = sub$power,
         xout = delta, rule = 2)$y
}

# target deltas
new_deltas <- c(0.4, 0.6, 0.8)

# build new rows
new_rows <- expand.grid(
  delta_AI_1 = new_deltas,
  num_allele_specific_reads = unique(h1_plot_data_expanded$num_allele_specific_reads)
) %>%
  mutate(
    num_bioreps = 600,
    power = case_when(
      delta_AI_1 <= 0.5 ~ mapply(interp_delta, delta_AI_1, num_allele_specific_reads,
                                 MoreArgs = list(df = h1_plot_data_expanded)),
      TRUE ~ 1
    )
  )

# combine
h1_all <- bind_rows(h1_plot_data_expanded, new_rows) %>%
  arrange(delta_AI_1, num_allele_specific_reads)

h1_all <- h1_all %>%
  filter(delta_AI_1 %in% c(0.2,0.4,0.6,0.8)) 

h1_all <- h1_all[!duplicated(h1_all),]


## Plot with delta AI x axis
gg_h1 <- ggplot(h1_all,
                aes(x = delta_AI_1, y = power,
                    color = factor(num_allele_specific_reads))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1), 
                     expand = c(0.01, 0.01),
                     labels = label_percent(accuracy = 1)) +
  scale_x_continuous(
    breaks = c(0.2, 0.4, 0.6, 0.8),
    labels = c("0.2", "0.4", "0.6", "0.8"),
    sec.axis = dup_axis(
      breaks = c(0.2, 0.4, 0.6, 0.8),
      labels = c("60:40", "70:30", "80:20", "90:10"),
      name = "Majority allele ratio"
    )
  ) +
  labs(x = expression(Delta*AI[1]),
       y = "Power",
       color = "Reads") +
  theme_classic() + theme(legend.position="bottom",
                          plot.title=element_blank()) +
  geom_hline(yintercept = 0.8, linetype = "dashed")

gg_h1

## Plot with Beta POE x axis
# helper to convert ΔAI -> Beta_POE
delta_to_beta <- function(delta) {
  theta <- 0.5 + delta/2
  log(theta / (1 - theta))
}

# set ΔAI and corresponding Beta_POE
delta_breaks <- c(0.2, 0.4, 0.6)
beta_breaks  <- delta_to_beta(delta_breaks)
beta_labels  <- round(beta_breaks, 3)

h1_all <- h1_all %>% 
  filter(delta_AI_1%in%c(0.2,0.4,0.6))

gg_h1 <- ggplot(h1_all,
                aes(x = delta_AI_1, y = power,
                    color = factor(num_allele_specific_reads))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1), 
                     expand = c(0.01, 0.01),
                     labels = label_percent(accuracy = 1)) +
  scale_x_continuous(
    breaks = delta_breaks,
    labels = beta_labels,
    sec.axis = dup_axis(
      breaks = delta_breaks,
      labels = c("60:40", "70:30", "80:20"),
      name = expression(bold("Majority allele ratio"))
    )
  ) +
  labs(x = expression(bold(beta[POE])),
       y = "Power",
       color = "Reads") +
  theme_classic(base_size = 14, base_family = "") +
  theme(
    text = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(face = "plain", size = 10),
    axis.text.y = element_text(face = "plain", size = 10),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(face = "bold", size = 8),
    legend.position = "bottom",
    legend.margin = margin(t = -5, unit = "pt"),        # tighten inside legend box
    legend.box.margin = margin(-5, 0, -5, 0, unit="pt"), # tighten legend vs plot
    plot.title = element_blank(),
    axis.line = element_blank(),  # remove default axis lines
    panel.border = element_rect(
      colour = "black",
      fill   = NA,
      linewidth = 0.5
    )) +
  geom_hline(yintercept = 0.8, linetype = "dashed")

gg_h1





# -------------------------
# H3: Power vs ΔAI3
# -------------------------
h3_plot_data <- df %>%
  # keep non-null effect sizes
  filter(delta_AI_3 > 0) %>%
  group_by(delta_AI_3, num_allele_specific_reads, num_bioreps) %>%
  summarise(power = mean(prop_H3_LE05, na.rm = TRUE), .groups = "drop")

gg_h3 <- ggplot(h3_plot_data,
                aes(x = delta_AI_3, y = power, color = factor(num_allele_specific_reads))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  labs(x = expression(Delta*AI[3]),
       y = "Power",
       color = "# of Allele Specific Reads",
       title = "Power vs ΔAI (H3)") +
  theme_classic() + theme(legend.position="bottom")

gg_h3


##### Interpolation
### --- Conversions ---
h3_plot_data <- h3_plot_data %>%
  mutate(
    theta = 0.5 + delta_AI_3/2,
    fc = theta / (1 - theta),
    beta_POE = log(fc)
  )

# Step 1: interpolate across reads at observed FCs
interp_reads <- function(reads, fc_val, df) {
  sub <- df %>% filter(abs(fc - fc_val) < 1e-6)
  approx(
    x = sub$num_allele_specific_reads,
    y = sub$power,
    xout = reads,
    rule = 2   # clamp
  )$y
}

new_reads <- c(50, 100, 150, 250, 500)
observed_fc <- unique(h3_plot_data$fc)

new_rows <- expand.grid(
  fc = observed_fc,
  num_allele_specific_reads = new_reads
) %>%
  mutate(
    num_bioreps = 600,
    power = mapply(interp_reads,
                   reads = num_allele_specific_reads,
                   fc_val = fc,
                   MoreArgs = list(df = h3_plot_data))
  )

# Combine and recalc all derived cols
h3_plot_data_expanded <- bind_rows(h3_plot_data, new_rows) %>%
  mutate(
    theta = fc / (1 + fc),
    delta_AI_3 = 2 * (theta - 0.5),
    beta_POE = log(fc)
  ) %>%
  arrange(fc, num_allele_specific_reads)


# Step 2: interpolate across FC for each read depth
interp_fc <- function(fc_target, read_val, df) {
  sub <- df %>% 
    filter(num_allele_specific_reads == read_val) %>%
    distinct(fc, .keep_all = TRUE)   # remove duplicates
  approx(
    x = sub$fc,
    y = sub$power,
    xout = fc_target,
    rule = 2
  )$y
}

target_fc <- c(1.5,2.5,4)

new_rows_fc <- expand.grid(
  fc = target_fc,
  num_allele_specific_reads = new_reads   # use new_reads, not unique(df)
) %>%
  mutate(
    num_bioreps = 600,
    power = mapply(interp_fc,
                   fc_target = fc,
                   read_val = num_allele_specific_reads,
                   MoreArgs = list(df = h3_plot_data_expanded))
  )

h3_all <- bind_rows(h3_plot_data_expanded, new_rows_fc) %>%
  arrange(fc, num_allele_specific_reads) %>%
  filter(fc %in% target_fc) %>%
  distinct()

### --- Add back delta + beta_POE for plotting ---
h3_all <- h3_all %>%
  mutate(
    theta = fc / (1 + fc),
    delta_AI_3 = 2 * (theta - 0.5),
    beta_POE = log(fc)
  )

h3_all <- h3_all %>% 
  filter(num_allele_specific_reads%in%c(50,100,250,500)) %>%
  filter(fc%in%c(1.5,2.5,4))

### --- Plot ---
gg_h3 <- ggplot(h3_all,
                aes(x = beta_POE, y = power,
                    color = factor(num_allele_specific_reads))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1),
                     expand = c(0.01, 0.01),
                     labels = label_percent(accuracy = 1)) +
  scale_x_continuous(
    breaks = log(target_fc),
    labels = round(log(target_fc), 3),
    sec.axis = dup_axis(
      breaks = log(target_fc),
      labels = target_fc,
      name = expression(bold("Fold change"))
    )
  ) +
  labs(
    x = expression(bold(beta[int])),
    y = "Power",
    color = "Reads"
  ) +
  theme_classic(base_size = 14, base_family = "") +
  theme(
    text = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(face = "plain", size = 10),
    axis.text.y = element_text(face = "plain", size = 10),
    legend.title = element_text(face = "bold", size = 9),
    legend.text  = element_text(face = "bold", size = 8),
    legend.position = "bottom",
    legend.margin = margin(t = -5, unit = "pt"),
    legend.box.margin = margin(-5, 0, -5, 0, unit = "pt"),
    plot.title = element_blank(),
    axis.line = element_blank(),  # remove default axis lines
    panel.border = element_rect(
      colour = "black",
      fill   = NA,
      linewidth = 0.5
    )) +
  geom_hline(yintercept = 0.8, linetype = "dashed")

gg_h3


# With delta AI x axis
gg_h3 <- ggplot(h3_all,
                aes(x = delta_AI_3, y = power,
                    color = factor(num_allele_specific_reads))) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0, 1), 
                     expand = c(0.01, 0.01),
                     labels = label_percent(accuracy = 1)) +
  scale_x_continuous(
    breaks = c(0.2, 0.4, 0.6, 0.8),
    labels = c("0.2", "0.4", "0.6", "0.8"),
    sec.axis = dup_axis(
      breaks = c(0.2, 0.4, 0.6, 0.8),
      labels = c("1.5", "2.3", "4.0", "9.0"),
      name = "Fold change"
    )
  ) +
  labs(x = expression(Delta*AI[3]),
       y = "Power",
       color = "Reads") +
  theme_classic() +
  theme(legend.position="bottom",plot.title=element_blank()) +
  geom_hline(yintercept = 0.8, linetype = "dashed")

gg_h3



### Combined

library(patchwork)

combined_plot <- gg_h1 + 
  (gg_h3 + theme(
    axis.title.y = element_text(color = NA), # invisible text, keeps space
    axis.text.y  = element_blank()           # remove tick labels
  )) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

combined_plot

### Optional title
combined_plot <- combined_plot +
  plot_annotation(
    title = "Power analysis",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 10,
        hjust = 0.5,
        margin = margin(b = 0)   # tighten bottom spacing
      )
    )
  )

combined_plot

ggsave(fig("power_analysis.png"), combined_plot, width = 3.5, height = 4.5, units = "in", dpi = 300)
