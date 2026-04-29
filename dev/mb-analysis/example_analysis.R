library(tidyverse)
library(data.table)
devtools::load_all()

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
  count(block_id) |>
  mutate(days = n/24) |>
  slice_max(n = 5, order_by = n)

# select one
example_block <- block_data |>
  filter(block_id == 2)

# plot
example_block |>
  arrange(datetime) |>
  pivot_longer(cols = c(o2, ph)) |>
  ggplot(aes(x = datetime, y = value, 
             color = name, group = block_id)) +
  geom_path(alpha = 0.4, linewidth = 0.2) +
  facet_grid(name ~ block_id, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  ggthm +
  labs(x = NULL, y = 'z score', title = 'Bay South') +
  guides(color = guide_none())

paste(img_out, 'bs1-20.png', sep = '/') |>
  ggsave(dpi = 400, width = 6, height = 4)

# (pre)smooth to remove tidal fluctuations
w_smooth_pre <- 25
example_block_presm <- example_block |>
  mutate(across(.cols = c(o2, ph), 
                .fns = list(smooth = ~frollmean(.x, 
                                                n = w_smooth_pre, 
                                                fill = NA, 
                                                align = 'center')),
                .names = "{.col}")) |>
  drop_na(o2, ph) |>
  arrange(datetime) |>
  group_by(time.rounded = floor_date(datetime, unit = '6 hours')) |>
  slice_min(datetime)

# plot
example_block_presm |>
  pivot_longer(cols = c(o2, ph)) |>
  ggplot(aes(x = datetime, y = value, 
             color = name, group = block_id)) +
  geom_path(alpha = 0.4, linewidth = 0.4) +
  facet_grid(name ~ block_id, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  ggthm +
  labs(x = NULL, y = 'z score', title = 'Bay South') +
  guides(color = guide_none())

paste(img_out, 'bs1-20-presm.png', sep = '/') |>
  ggsave(dpi = 400, width = 6, height = 4)

# fit null model
x1 <- example_block_presm$o2
x2 <- example_block_presm$ph
dt <- example_block_presm$datetime
fit <- lomad_fit(x1, x2,
                 method = "state",
                 rho0   = 0,
                 q      = 14*4,
                 h      = 14*2*4,
                 max_pq = 2)

# inspect
fit$null_model

# plot fit and detected decoupling periods
paste(img_out, 'bs1-20-fit.png', sep = '/') |>
  png(width = 5, height = 4, units = 'in', res = 400)
lomad_plot(x1, x2, fit, dates = dt, alpha = 0.4)
dev.off()

## END ------------

# observed decoupling statistics and (rough) expectations
fit$observed
fit$expected_asymptotic

# asymptotic test (very rough... approximation of approximation)
test1 <- lomad_test(fit, method = "analytic")
test1$p_values

# state process bootstrap (quick)
test2 <- lomad_test(fit, method = "mc", B = 1000, seed = 123)
test2$p_values

# full parametric bootstrap (takes a bit)
test3 <- lomad_test(fit, method = "boot", B = 500, seed = 123,
                    ncores = parallel::detectCores() - 1,
                    verbose = TRUE)
test3$p_values
test3$expected
test3$observed
