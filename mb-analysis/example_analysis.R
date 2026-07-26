library(tidyverse)
library(lomad)

# output directory for figures
img_out <- '_img'
fs::dir_create(img_out)

# ggplot theming
ggthm <- theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = 'black', linewidth = 0.1))

# read in data (blocked and standardized)
block_data <- read_csv('_mb-data/ph_o2_blocks.csv')

# longest blocks
block_data |>
  count(block_id, location) |>
  mutate(days = n/24) |>
  slice_max(n = 5, order_by = n)

# select one (BS1 block 2)
example_block <- block_data |>
  filter(block_id == 2, location == "BS1") |>
  arrange(datetime)

# plot raw series
example_block |>
  pivot_longer(cols = c(o2, ph)) |>
  ggplot(aes(x = datetime, y = value, 
             color = name, group = block_id)) +
  geom_path(alpha = 0.4, linewidth = 0.2) +
  facet_grid(name ~ block_id, scales = 'free_x', space = 'free_x') +
  ggthm +
  labs(x = NULL, y = 'z score', title = 'Bay South') +
  guides(color = guide_none())

paste(img_out, 'bs1-20.png', sep = '/') |>
  ggsave(dpi = 400, width = 6, height = 4)

# --- Presmooth ---------------------------------------------------------------
# Remove tidal periodicity (~25h cycle at hourly) and downsample to 6h

source('mb-analysis/utils.R')

example_block_presm <- presmooth_tidal(example_block,
                                        cols = c('o2', 'ph'),
                                        cycle_len = 25,
                                        step = '6 hours')

# plot presmoothed
example_block_presm |>
  pivot_longer(cols = c(o2, ph)) |>
  ggplot(aes(x = datetime, y = value, 
             color = name, group = block_id)) +
  geom_path(alpha = 0.4, linewidth = 0.4) +
  facet_grid(name ~ block_id, scales = 'free_x', space = 'free_x') +
  ggthm +
  labs(x = NULL, y = 'z score', title = 'Bay South (presmoothed)') +
  guides(color = guide_none())

paste(img_out, 'bs1-20-presm.png', sep = '/') |>
  ggsave(dpi = 400, width = 6, height = 4)

# --- Fit and test (CLT method with ARMA noise) --------------------------------

x1 <- example_block_presm$o2
x2 <- example_block_presm$ph
dt <- example_block_presm$datetime

fit <- lomad_fit(x1, x2, h = 3, s = 50, noise_method = "arma")
tst <- lomad_test(fit, alpha = 0.05)

cat(sprintf("Rejected %d / %d time points at alpha_eff = %.4f\n",
            sum(tst$rejected, na.rm = TRUE),
            sum(!is.na(tst$rejected)),
            tst$alpha_eff))

# plot fit and detected decoupling periods
paste(img_out, 'bs1-20-fit.png', sep = '/') |>
  png(width = 5, height = 4, units = 'in', res = 400)
lomad_plot(fit, tst, dates = dt, alpha = 0.4)
dev.off()
