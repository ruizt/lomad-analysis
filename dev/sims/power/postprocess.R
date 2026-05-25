## postprocess.R — summarise and plot power study results
##
## Reads the assembled results.rds produced by collect.R and generates
## summary tables and power curve plots.
##
## Usage (from the repo root):
##   Rscript dev/sims/power/postprocess.R

library(dplyr)
library(ggplot2)

OUT_DIR <- "dev/sims/power/results"
IMG_DIR <- file.path(OUT_DIR, "_img")
alpha   <- 0.05

# ---- Load --------------------------------------------------------------------

results <- readRDS(file.path(OUT_DIR, "results.rds"))

# ---- Summary -----------------------------------------------------------------

results_summary <- results |>
  group_by(struct, d, n, snr, phi, method) |>
  summarise(
    S         = n(),
    detection = mean(detected, na.rm = TRUE),
    ci_lo     = qbeta(0.025, sum(detected), S - sum(detected) + 1),
    ci_hi     = qbeta(0.975, sum(detected) + 1, S - sum(detected)),
    .groups   = "drop"
  )

# saveRDS(results_summary, file.path(OUT_DIR, "results_summary.rds"))

# ---- Draft figure: power curves faceted by snr x phi -------------------------

# dir.create(IMG_DIR, showWarnings = FALSE, recursive = TRUE)

df <- filter(results_summary, n == 600)

ggplot(df, 
       aes(d, detection, 
           colour = struct,
           linetype = method, 
           group = interaction(struct, method, n))) +
  geom_hline(yintercept = alpha, 
             linetype = "dashed", 
             colour = "grey60", 
             linewidth = 0.5) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = struct),
              alpha = 0.2, color = NULL) +
  # # geom_point(size = 2.5) +
  # scale_x_continuous(breaks = 0:2) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c("smooth" = "#0072B2",
                                 "cross"  = "#D55E00",
                                 "rate"   = "#009E73")) +
  # scale_linetype_manual(values = c("estimated" = "solid",
  #                                  "oracle"    = "dashed")) +
  facet_grid(snr ~ phi, labeller = labeller(
    phi = \(x) paste0("phi == ", x),
    snr = \(x) paste0("SNR == ", x),
    .default = label_parsed
  )) +
  labs(x        = "Separation (d)",
       y        = "Power",
       colour   = "Structure",
       fill     = "Structure") +
  guides(linetype = guide_none()) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, color = 'darkgray'))

# ---- Nested facet layout (all sample sizes in one plot) ----------------------

library(ggh4x)

results_summary |>
  mutate(struct = str_to_sentence(struct)) |>
ggplot(aes(d, detection, colour = struct,
           linetype = method,
           group = interaction(struct, method))) +
  geom_hline(yintercept = alpha, linetype = "dashed",
             colour = "grey60", linewidth = 0.5) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = struct),
              alpha = 0.2, colour = NA) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c("Smooth" = "#0072B2",
                                 "Cross"  = "#D55E00",
                                 "Rate"   = "#009E73")) +
  scale_fill_manual(values = c("Smooth" = "#0072B2",
                               "Cross"  = "#D55E00",
                               "Rate"   = "#009E73")) +
  facet_nested(phi ~ snr + n, labeller = labeller(
    phi = \(x) paste0("phi == ", x),
    snr = \(x) paste0("SNR == ", x),
    n   = \(x) paste0("T == ", x),
    .default = label_parsed
  )) +
  labs(x = "Separation (d)", y = "Power",
       colour = "Structure", fill = "Structure") +
  guides(linetype = guide_none()) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        axis.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, color = "darkgray"))


ggsave(file.path(IMG_DIR, "power_curves.png"), 
       width = 9, height = 4, dpi = 450)
