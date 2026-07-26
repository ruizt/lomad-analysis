
##prototype localization metrics using stored series output from the power study
##
##Compare continuous ground truth w series to inferential classifications
##from p-value/rejection series.

library(tidyverse)

# ---------------- data-----------
RAW_DIR <- "sims/power/results/_simulations-power/_raw"
alpha <- 0.05

# length(series_files) is 360

series_files <- list.files(
  RAW_DIR,
  pattern = "-series\\.rds$",
  full.names = TRUE
)

# ---- Convert one replicate to data frame -------------------------------------

# the test cannot evaluate all time points because rolling windows so we need 
#enough observations on both sides. So we only analyze these incies

series_to_df <- function(rep_series) {
  
  w_full <- rep_series$w
  vi <- rep_series$vi
  p_raw_full <- rep_series$p_raw
  p_adj_full <- rep_series$p_adj
  rejected_full <- rep_series$rejected
  
  tibble(
    idx = seq_along(vi),
    #original time index in simulated series
    t = vi,
    # ground truth coupling value at valid time points
    w = w_full[vi],
    # larger values mean more decoupling 
    decouple = abs(1 - w_full[vi]),
    # raw unadjusted p-values at valid time points
    p_raw = p_raw_full[vi],
    p_adj = p_adj_full[vi],
    rejected = rejected_full[vi]
  )
}


# ---- Ground truth conditional on classification ------------------------------

# among windows classified as reject or not rejected, how strong was 
# the true decoupling? - precision/FDR direction

# Instead of immediately calling something a true/false discovery,
# we summarize the continuous ground truth within each inferred class.
truth_given_class <- function(df) {
  df |>
    filter(!is.na(rejected)) |>
    group_by(rejected) |>
    summarise(
      n_windows = n(),
      mean_decouple = mean(decouple, na.rm = TRUE),
      median_decouple = median(decouple, na.rm = TRUE),
      q75_decouple = quantile(decouple, 0.75, na.rm = TRUE),
      q90_decouple = quantile(decouple, 0.90, na.rm = TRUE),
      max_decouple = max(decouple, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- Classification conditional on ground truth ------------------------------
# If true decoupling is larger than cutoff c, how often does the method reject
# sensitiity/specificcity direction

#for each cutoff c: true_decouples = decouple > c
# then compute p(rejected | decouple > c)
# p(not rejected | decouple <= c)

class_given_truth <- function(df, c_grid = seq(0, 1, by = 0.05)) {
  
  df <- df |>
    filter(!is.na(rejected))
  
  map_dfr(c_grid, function(cutoff) {
    
    true_decoupled <- df$decouple > cutoff
    
    tibble(
      c = cutoff,
      # number of windows above and below the truth cutoff
      n_true = sum(true_decoupled, na.rm = TRUE),
      n_not_true = sum(!true_decoupled, na.rm = TRUE),
      # among windows with strong true decoupling, how many were flagged?
      p_flag_given_true = ifelse(
        sum(true_decoupled, na.rm = TRUE) > 0,
        mean(df$rejected[true_decoupled], na.rm = TRUE),
        NA_real_
      ),
      # among windows without strong true decouplng how many were not flagged?
      p_no_flag_given_not_true = ifelse(
        sum(!true_decoupled, na.rm = TRUE) > 0,
        mean(!df$rejected[!true_decoupled], na.rm = TRUE),
        NA_real_
      ),
      # overall proportion of valid windows flagged
      flagged_rate = mean(df$rejected, na.rm = TRUE)
    )
  })
}

# ---------- count rejections ------------
# one series file contains 500 simulation replicates
# this function summarizes each replicate byhow many valid windows rejected 
#
#find a replicate with a mixed pattern:some rejected and some unrejected windows
count_rejections <- function(series_obj) {
  imap_dfr(series_obj, function(rep_series, seed_name) {
    
    df <- series_to_df(rep_series)
    
    tibble(
      seed = seed_name,
      n_valid = nrow(df),
      n_rejected = sum(df$rejected, na.rm = TRUE),
      prop_rejected = mean(df$rejected, na.rm = TRUE),
      mean_decouple = mean(df$decouple, na.rm = TRUE),
      max_decouple = max(df$decouple, na.rm = TRUE)
    )
  })
}

# ---------- summarize rejection behavior for one file -------------
# function reads one -series.rds file and summarizes the rejection of
# proportions across 500 replicates

summarize_file_rejections <- function(file) {
  
  obj <- readRDS(file)
  
  counts <- count_rejections(obj)
  
  tibble(
    file = file,
    basename = basename(file),
    n_reps = nrow(counts),
    mean_prop_rejected = mean(counts$prop_rejected, na.rm = TRUE),
    median_prop_rejected = median(counts$prop_rejected, na.rm = TRUE),
    min_prop_rejected = min(counts$prop_rejected, na.rm = TRUE),
    max_prop_rejected = max(counts$prop_rejected, na.rm = TRUE)
  )
}


# ---------- Choose Candidate Files ---------------
# start with one favorable broad setting 
# structure = smooth n =600, SNR = 1.5, phi = 0.3, method = estimated
# 
# within that broad setting we compare d values. 
# want a d value where rejection is mized because localization can't be studied
#if every window is rejected or no windows rejected
#
candidate_files <- series_files[
  str_detect(basename(series_files), "smooth") &
    str_detect(basename(series_files), "n600") &
    str_detect(basename(series_files), "snr1-5") &
    str_detect(basename(series_files), "phi0-3") &
    !str_detect(basename(series_files), "oracle")
]

basename(candidate_files)

candidate_summary <- map_dfr(candidate_files, summarize_file_rejections)

candidate_summary |>
  arrange(mean_prop_rejected)

#--------- Analyze one mixed-rejection setting ------------
mixed_file <- candidate_summary |>
  filter(str_detect(basename, "d0-5")) |>
  pull(file)

basename(mixed_file)

mixed <- readRDS(mixed_file)
# count rejection proportions for all 500 replications in this file
mixed_counts <- count_rejections(mixed)

# show replicates closest to 50% rejected
# these are so useful because contain both rejection and unrejected windows 
mixed_counts |>
  arrange(abs(prop_rejected - 0.5)) |>
  print(n = 20)

# choose the replicate whose rejection proportion is closest to 0.5. 
mixed_seed <- mixed_counts |>
  arrange(abs(prop_rejected - 0.5)) |>
  slice(1) |>
  pull(seed)

mixed_seed

# convert that replicate into a clean time-indexed data frame
mixed_df <- series_to_df(mixed[[mixed_seed]])

table(mixed_df$rejected, useNA = "ifany")
summary(mixed_df$decouple)
summary(mixed_df$p_adj)

# --------- localization summaries for the mixed replicate ----------
# ground trut conditional on classification comparing true fecoupling strength
#for rejected vs unrejected windows
truth_given_class(mixed_df)

# for a grid of cutoffs c, compute p(rejected | |1-w|> c)
mixed_curve <- class_given_truth(mixed_df)
print(mixed_curve, n = 21)

# -----------plot 1: true decoupling over time with rejected windows marked-----
mixed_df |>
  ggplot(aes(x = t, y = decouple)) +
  geom_line(linewidth = 0.4) +
  geom_point(
    data = \(x) filter(x, rejected),
    size = 1.2
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Localization example: mixed rejection pattern",
    subtitle = paste(basename(mixed_file), "| seed", mixed_seed),
    x = "Time index",
    y = expression(abs(1 - w[t]))
  )

#------- plot 2: clasification conditioal on ground truth----------------------
mixed_curve |>
  filter(!is.na(p_flag_given_true)) |>
  ggplot(aes(x = c, y = p_flag_given_true)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Classification conditional on ground truth",
    subtitle = expression(P(rejected ~ "|" ~ abs(1 - w[t]) > c)),
    x = "Ground-truth decoupling threshold c",
    y = "Probability rejected"
  )

#------------- plot 3: specificity-like localization curve----------
mixed_curve |>
  filter(!is.na(p_no_flag_given_not_true)) |>
  ggplot(aes(x = c, y = p_no_flag_given_not_true)) +
  geom_line() +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Specificity-like localization curve",
    subtitle = expression(P(not~rejected ~ "|" ~ abs(1 - w[t]) <= c)),
    x = "Ground-truth decoupling threshold c",
    y = "Probability not rejected"
  )