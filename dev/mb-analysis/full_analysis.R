library(tidyverse)
devtools::load_all()
source('dev/mb-analysis/utils.R')   # presmooth_tidal()

# output directory for figures
img_out <- '_img'
fs::dir_create(img_out)

# ggplot theming
ggthm <- theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = 'black', linewidth = 0.1))

# --- Data -------------------------------------------------------------------

block_data <- read_csv('_mb-data/ph_o2_blocks.csv')

block_data |>
  count(block_id) |>
  mutate(days = n / 24) |>
  arrange(desc(n))

# --- Presmooth --------------------------------------------------------------
# Remove tidal periodicity (~25 h cycle at hourly resolution) and downsample
# to 6-hour resolution. Applied independently within each block.

blocks_presm <- block_data |>
  group_split(block_id) |>
  lapply(presmooth_tidal, cols = c('o2', 'ph'), cycle_len = 25, step = '6 hours')

names(blocks_presm) <- vapply(blocks_presm, \(b) as.character(b$block_id[[1]]), character(1))

# --- Fit null model across all blocks ---------------------------------------

blocks <- lapply(blocks_presm, \(b) list(x1 = b$o2, x2 = b$ph))

fit <- lomad_fit_blocks(blocks,
                        rho0   = 0,
                        q      = 14 * 4,
                        h      = 14 * 2 * 4,
                        max_pq = 2)

fit$null_model
fit$observed
fit$expected_asymptotic

# --- Visualize individual blocks --------------------------------------------
# plot_lomad_fit() requires x1/x2 and fit$R/I/valid_idx to have the same
# length, so we build a per-block view of the fit object for each plot.

block_fit_view <- function(fit, k) {
  blk <- fit$blocks[[k]]
  fit$trend_hat <- blk$trend_hat
  fit$ma1       <- blk$ma1
  fit$ma2       <- blk$ma2
  fit$R         <- blk$R
  fit$I         <- blk$I
  fit$valid_idx <- blk$valid_idx   # local indices, correct for this block's length
  fit
}

pdf(paste(img_out, 'all-blocks-fit.pdf', sep = '/'), width = 5, height = 4)
for (nm in names(fit$blocks)) {
  plot_lomad_fit(blocks[[nm]]$x1, blocks[[nm]]$x2,
                 fit  = block_fit_view(fit, nm),
                 alpha = 0.4)
  title(main = nm, line = 0.5)
}
dev.off()

# --- Inference --------------------------------------------------------------

# Asymptotic test (rough approximation)
test_analytic <- lomad_test_analytic(fit)
test_analytic$p_values

# Markov-chain bootstrap (fast)
test_mc <- lomad_test_mc(fit, B = 1000, seed = 4721, 
                         verbose = TRUE)
test_mc$p_values

# Full parametric bootstrap (use parallel cores to speed up)
test_boot <- lomad_test_boot(fit, B = 500, seed = 4721,
                             ncores = parallel::detectCores() - 1,
                             verbose = TRUE)
test_boot$p_values
test_boot$expected
test_boot$observed
