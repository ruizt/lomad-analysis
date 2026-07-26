## postprocess.R — summarise and plot power study results
##
## Reads the assembled results.rds produced by collect.R and generates
## summary tables and power curve plots.
##
## Usage (from the repo root):
##   Rscript sims/power/postprocess.R

library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

OUT_DIR <- "sims/power/results"
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


# ---- Localization helpers ---------------------------------------------------

RAW_DIR <- file.path(OUT_DIR, "_raw")

#' Convert one replicate's series data to a tidy data frame.
#' Ground truth is sep_t = |x1_t - x2_t| at valid indices.
series_to_df <- function(rep_series) {
  vi <- rep_series$vi
  tibble(
    t        = vi,
    sep      = rep_series$sep[vi],
    w        = if (!is.null(rep_series$w)) rep_series$w[vi] else NA_real_,
    rejected = rep_series$rejected[vi]
  )
}

#' Sensitivity curve: P(rejected | ground_truth > c) as a function of c.
compute_sens_curve <- function(pool, col, c_grid) {
  vals <- pool[[col]]
  rej  <- pool$rejected
  vapply(c_grid, function(c) {
    above <- vals > c
    if (sum(above, na.rm = TRUE) == 0L) return(NA_real_)
    mean(rej[above], na.rm = TRUE)
  }, double(1))
}

#' Specificity curve: P(not rejected | ground_truth <= c) as a function of c.
compute_spec_curve <- function(pool, col, c_grid) {
  vals <- pool[[col]]
  rej  <- pool$rejected
  vapply(c_grid, function(c) {
    below <- vals <= c
    if (sum(below, na.rm = TRUE) == 0L) return(NA_real_)
    mean(!rej[below], na.rm = TRUE)
  }, double(1))
}

#' Summarise ground-truth separation within each rejection class.
truth_given_class <- function(df) {
  df |>
    filter(!is.na(rejected)) |>
    group_by(rejected) |>
    summarise(
      n_windows  = n(),
      mean_sep   = mean(sep, na.rm = TRUE),
      median_sep = median(sep, na.rm = TRUE),
      q75_sep    = quantile(sep, 0.75, na.rm = TRUE),
      q90_sep    = quantile(sep, 0.90, na.rm = TRUE),
      max_sep    = max(sep, na.rm = TRUE),
      .groups    = "drop"
    )
}

#' Parse a series filename into simulation parameters.
parse_series_filename <- function(fname) {
  bn <- basename(fname) |> str_remove("-series\\.rds$")
  oracle <- str_detect(bn, "-oracle")
  bn <- str_remove(bn, "-oracle")
  parts <- str_match(
    bn,
    "^(\\w+)_d([0-9-]+)_n(\\d+)_snr([0-9-]+)_phi([0-9-]+)$"
  )
  tibble(
    file   = fname,
    struct = parts[, 2],
    d      = as.numeric(str_replace(parts[, 3], "-", ".")),
    n      = as.integer(parts[, 4]),
    snr    = as.numeric(str_replace(parts[, 5], "-", ".")),
    phi    = as.numeric(str_replace(parts[, 6], "-", ".")),
    oracle = oracle
  )
}

#' Pool all replicates in a series object into one data frame.
pool_reps <- function(series_obj) {
  imap_dfr(series_obj, function(rep, seed_name) {
    if (is.null(rep)) return(NULL)
    series_to_df(rep) |> mutate(seed = seed_name)
  })
}

# ---- Localization: load series files ----------------------------------------

series_files <- list.files(RAW_DIR, pattern = "-series\\.rds$", full.names = TRUE)

file_meta <- map_dfr(series_files, parse_series_filename) |>
  filter(d > 0, !oracle)

# ---- Localization: pilot on two settings ------------------------------------

pick_file <- function(meta, struct_, d_, n_, snr_, phi_) {
  meta |>
    filter(struct == struct_, d == d_, n == n_, snr == snr_, phi == phi_) |>
    pull(file)
}

fav_file   <- pick_file(file_meta, "smooth", 0.5, 600, 1.5, 0.3)
unfav_file <- pick_file(file_meta, "smooth", 0.5, 600, 0.5, 0.3)

fav_series   <- readRDS(fav_file)
unfav_series <- readRDS(unfav_file)

pool_fav   <- pool_reps(fav_series)
pool_unfav <- pool_reps(unfav_series)

# ---- Localization: truth given class ----------------------------------------

cat("\n--- Favorable setting: truth given class ---\n")
print(truth_given_class(pool_fav))

cat("\n--- Unfavorable setting: truth given class ---\n")
print(truth_given_class(pool_unfav))

# ---- Localization: sensitivity curves ---------------------------------------

c_grid <- seq(0, max(c(pool_fav$sep, pool_unfav$sep), na.rm = TRUE),
              length.out = 200)

sens_df <- bind_rows(
  tibble(c = c_grid,
         p_reject = compute_sens_curve(pool_fav, "sep", c_grid),
         setting  = "Favorable (snr=1.5)"),
  tibble(c = c_grid,
         p_reject = compute_sens_curve(pool_unfav, "sep", c_grid),
         setting  = "Unfavorable (snr=0.5)")
)

ggplot(sens_df, aes(x = c, y = p_reject, colour = setting)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = alpha, linetype = "dashed", colour = "grey60") +
  labs(x = expression("Threshold " * c),
       y = expression(P(rejected ~ "|" ~ sep[t] > c)),
       title = "Sensitivity curve — pointwise trend separation",
       colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# ---- Localization: specificity curves ---------------------------------------

spec_df <- bind_rows(
  tibble(c = c_grid,
         p_notrej = compute_spec_curve(pool_fav, "sep", c_grid),
         setting  = "Favorable (snr=1.5)"),
  tibble(c = c_grid,
         p_notrej = compute_spec_curve(pool_unfav, "sep", c_grid),
         setting  = "Unfavorable (snr=0.5)")
)

ggplot(spec_df, aes(x = c, y = p_notrej, colour = setting)) +
  geom_line(linewidth = 0.8) +
  labs(x = expression("Threshold " * c),
       y = expression(P(not~rejected ~ "|" ~ sep[t] <= c)),
       title = "Specificity curve — pointwise trend separation",
       colour = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

# ---- Localization: individual replicate visualisation -----------------------

rep_stats <- pool_fav |>
  group_by(seed) |>
  summarise(prop_rej = mean(rejected, na.rm = TRUE), .groups = "drop")

mixed_seed <- rep_stats |>
  arrange(abs(prop_rej - 0.5)) |>
  slice(1) |>
  pull(seed)

mixed_df <- pool_fav |> filter(seed == mixed_seed)

ggplot(mixed_df, aes(x = t, y = sep)) +
  geom_line(linewidth = 0.4) +
  geom_point(data = \(x) filter(x, rejected), size = 1.2, colour = "#D55E00") +
  labs(x = "Time index",
       y = expression("|" * x[1*t] - x[2*t] * "|"),
       title = "Localization: trend separation with rejected windows",
       subtitle = paste("seed", mixed_seed)) +
  theme_minimal(base_size = 13)
